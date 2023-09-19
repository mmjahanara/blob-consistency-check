use std::env::set_var;
use std::env::var;

use ff::{Field, PrimeField};
use halo2_base::QuantumCell::Constant;
use halo2_base::halo2_proofs::{
    dev::MockProver, 
    halo2curves::bls12_381::Scalar, 
    halo2curves::bn256::{Fr, Fq},
};
use halo2_base::safe_types::GateInstructions;
use halo2_base::utils::ScalarField;
use halo2_base::AssignedValue;
use halo2_ecc::bigint::ProperCrtUint;
use halo2_ecc::fields::fp::FpChip;
use halo2_ecc::fields::FieldChip;

use halo2_base::gates::builder::{GateThreadBuilder, RangeWithInstanceCircuitBuilder};
use halo2_base::gates::{RangeChip, GateChip};
use halo2_base::Context;

use poseidon::PoseidonChip;
use serde::{Deserialize, Serialize};

const T: usize = 3;
const RATE: usize = 2;
const R_F: usize = 8;
const R_P: usize = 57;
const BLOB_WIDTH: usize = 4;

const K: usize = 10;

const LIMB_BITS: usize = 88;
const NUM_LIMBS: usize = 3;

const LOOKUP_BITS: usize = 9;

fn main() {
    let input = CircuitInput {
        batch_commit: "77".to_string(),
        blob: [
            "1".to_string(),
            "2".to_string(),
            "3".to_string(),
            "4".to_string(),
        ],
    };

    set_var("LOOKUP_BITS", LOOKUP_BITS.to_string());

    let mut builder = GateThreadBuilder::<Fr>::mock();
    let ctx = builder.main(0);
    let mut make_public: Vec<AssignedValue<Fr>> = vec![];
    blob_consistency_check::<Fr>(ctx, input, &mut make_public);

    builder.config(K, Some(10));
    let circuit = RangeWithInstanceCircuitBuilder::mock(builder, make_public.clone());
        
    let public: Vec<Fr> = make_public.iter().map(|x|*x.value()).collect();
    MockProver::run(K as u32, &circuit, vec![public])
        .unwrap()
        .assert_satisfied();
}

fn blob_consistency_check<F: ScalarField>(
    ctx: &mut Context<F>,
    input: CircuitInput,
    make_public: &mut Vec<AssignedValue<F>>,
) {
    // === STEP 1 ===
    // load batch_commit and blob then compute the challenge point

    let batch_commit =
        F::from_str_vartime(&input.batch_commit).expect("failed to deserialize batch_commit");
    let batch_commit = ctx.load_witness(batch_commit);
    make_public.push(batch_commit);

    let blob: Vec<_> = input
        .blob
        .iter()
        .map(|x| F::from_str_vartime(x).expect("failed to deserialize blob"))
        .collect();
    let blob = blob
        .iter()
        .map(|x| ctx.load_witness(*x))
        .collect::<Vec<_>>();

    let range = RangeChip::<F>::default(LOOKUP_BITS);
    let gate = &range.gate;

    let mut poseidon = PoseidonChip::<F, T, RATE>::new(ctx, R_F, R_P).unwrap();
    poseidon.update(&[batch_commit]);
    poseidon.update(blob.as_slice());

    let challenge_point = poseidon.squeeze(ctx, gate).unwrap();

    // === STEP 2 === 
    // load challenge_point and blob in fp_chip and compute the evaluation

    // TODO: in the following, change Fq to bls12_381::Scalar 

    let fq_chip = FpChip::<F, Fq>::new(&range, LIMB_BITS, NUM_LIMBS);
    let one_fq = fq_chip.load_constant(ctx, Fq::one());

    // load challenge_point and blob to fq
    let challenge_point_fq = cross_field_load_private(ctx, &fq_chip, &challenge_point);

    let blob_fq: Vec<_> = input.blob.iter().map(|x| Fq::from_str_vartime(x).expect("failed to deserialize blob")).collect();
    let blob_fq = blob_fq.iter().map(|x| fq_chip.load_private(ctx, *x)).collect::<Vec<_>>();
    for i in 0..BLOB_WIDTH {
        cross_field_constrain_equal(ctx, gate, &blob[i], &blob_fq[i]);
    }
    
    // TODO: set `blob_width_th_root_of_unity` correctly
    let blob_width_th_root_of_unity= Fq::from(2);
    let roots_of_unity: Vec<_> = (0..BLOB_WIDTH)
        .map(|i| blob_width_th_root_of_unity.pow(&[i as u64, 0,0,0]))
        .collect();

    let roots_of_unity = roots_of_unity.iter()
        .map(|x| fq_chip.load_constant(ctx, *x))
        .collect::<Vec<_>>();
    
    // TODO: apply `bit_reversal_permutation` to `roots_of_unity`

    let mut result = fq_chip.load_constant(ctx, Fq::zero());
    let mut cp_is_not_root_of_unity = fq_chip.load_constant(ctx, Fq::one());
    let mut barycentric_evaluation = fq_chip.load_constant(ctx, Fq::zero());
    for i in 0..BLOB_WIDTH as usize {
        let numinator_i = fq_chip.mul(ctx, roots_of_unity[i].clone(), blob_fq[i].clone());
        
        let denominator_i_no_carry = fq_chip.sub_no_carry(ctx, challenge_point_fq.clone(), roots_of_unity[i].clone());
        let denominator_i = fq_chip.carry_mod(ctx, denominator_i_no_carry);
        
        // avoid division by zero
        // safe_denominator_i = denominator_i   (denominator_i != 0)
        // safe_denominator_i = 1               (denominator_i == 0)
        let is_zero_denominator_i = fq_is_zero(ctx, gate, &denominator_i);
        let is_zero_denominator_i = cross_field_load_private(ctx, &fq_chip, &is_zero_denominator_i);
        let safe_denominator_i = fq_chip.add_no_carry(ctx, denominator_i, is_zero_denominator_i.clone());
        let safe_denominator_i = fq_chip.carry_mod(ctx, safe_denominator_i);

        // update `cp_is_not_root_of_unity`
        cp_is_not_root_of_unity = fq_chip.mul(ctx, cp_is_not_root_of_unity, is_zero_denominator_i.clone());

        // update `result`, select blob_fq[i] if denominator_i == 0
        let select_blob_fq = fq_chip.mul(ctx, blob_fq[i].clone(), is_zero_denominator_i);
        let tmp_result = fq_chip.add_no_carry(ctx, result, select_blob_fq);
        result = fq_chip.carry_mod(ctx, tmp_result);

        let term_i = fq_chip.divide(ctx, numinator_i, safe_denominator_i);

        let evaluation_not_proper = fq_chip.add_no_carry(ctx, barycentric_evaluation, term_i);
        barycentric_evaluation = fq_chip.carry_mod(ctx, evaluation_not_proper);
    }
    // evaluation = evaluation * (challenge_point**BLOB_WIDTH - 1) / BLOB_WIDTH
    let cp_to_the_width = fq_pow(ctx, &fq_chip, &challenge_point_fq, BLOB_WIDTH as u32);
    let cp_to_the_width_minus_one = fq_chip.sub_no_carry(ctx, cp_to_the_width, one_fq);
    let cp_to_the_width_minus_one = fq_chip.carry_mod(ctx, cp_to_the_width_minus_one);
    let width_fq = fq_chip.load_constant(ctx, Fq::from(BLOB_WIDTH as u64));
    let factor = fq_chip.divide(ctx, cp_to_the_width_minus_one, width_fq);
    barycentric_evaluation = fq_chip.mul(ctx, barycentric_evaluation, factor);

    let select_evaluation = fq_chip.mul(ctx, barycentric_evaluation, cp_is_not_root_of_unity);
    let tmp_result = fq_chip.add_no_carry(ctx, result, select_evaluation);
    result = fq_chip.carry_mod(ctx, tmp_result);
    make_public.extend(result.limbs());    
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CircuitInput {
    pub batch_commit: String,
    pub blob: [String; BLOB_WIDTH],
}

fn cross_field_load_private<F: ScalarField>(
    ctx: &mut Context<F>, 
    fq_chip: &FpChip<F, Fq>, 
    x: &AssignedValue<F>, 
) -> ProperCrtUint<F> {
    let x_fq = Fq::from_bytes_le(x.value().to_bytes_le().as_slice());
    let x_fq = fq_chip.load_private(ctx, x_fq); 
    cross_field_constrain_equal(ctx, &fq_chip.range().gate, x, &x_fq);
    x_fq
}

fn cross_field_constrain_equal<F: ScalarField>(
    ctx: &mut Context<F>, 
    gate: &GateChip<F>, 
    x: &AssignedValue<F>, 
    x_fq: &ProperCrtUint<F>
) {
    //TODO: add overflow protection constraints
    //Bound the last limb of x_fq to be less than 2^LIMB_BITS-10)?

    let x_fq_limbs = x_fq.limbs();
    let mut sum = ctx.load_zero();
    let mut mul = ctx.load_constant(F::from(1));
    let limb_multiplier = ctx.load_constant(F::from_u128(2u128.pow(LIMB_BITS as u32)));
    for limb in x_fq_limbs {
        sum = gate.mul_add(ctx, limb.clone(), mul, sum);
        mul = gate.mul(ctx, limb_multiplier, mul);
    }
    ctx.constrain_equal(&sum, &x);
}

fn fq_is_zero<F: ScalarField>(
    ctx: &mut Context<F>, 
    gate: &GateChip<F>, 
    x_fq: &ProperCrtUint<F>, 
) -> AssignedValue<F> {
    let zero = ctx.load_zero();
    let x_fq_limbs = x_fq.limbs();
    let mut partial_and = ctx.load_constant(F::from(1));
    for limb in x_fq_limbs {
        let is_zero_limb = gate.is_equal(ctx, limb.clone(), zero);
        partial_and = gate.and(ctx, is_zero_limb, Constant(F::from(1)));
    }
    partial_and
}

fn fq_pow<F: ScalarField>(
    ctx: &mut Context<F>,
    fq_chip: &FpChip<F, Fq>,
    x: &ProperCrtUint<F>,
    pow: u32
)-> ProperCrtUint<F>{
    if pow == 0 {
        return fq_chip.load_constant(ctx, Fq::one());
    }
    else if pow == 1 {
        return x.clone();
    }
    
    let mut result = fq_pow(ctx, fq_chip, x, pow/2);
    result = fq_chip.mul(ctx, result.clone(), result);
    if pow % 2 == 1 {
        result = fq_chip.mul(ctx, result, x.clone());
    }
    result
}

fn bit_reversal_permutation<T: Clone + Copy>(seq: Vec<T>) -> Vec<T> {
    // return a permutation of seq, where the indices are bit-reversed
    // e.g. bit_reversal_permutation([0, 1, 2, 3]) = [0, 2, 1, 3]
    let n = seq.len();
    let log_n = (n as f64).log2() as usize;
    let mut result: Vec<T> = vec![seq[0]; n];
    for i in 0..n {
        let mut j = i;
        let mut k = 0;
        for _ in 0..log_n {
            k = (k << 1) | (j & 1);
            j >>= 1;
        }
        result[i] = seq[k];
    }
    result
}

#[test]
fn test_bit_reversal() {
    let seq = vec![0, 1, 2, 3];
    let expected = vec![0, 2, 1, 3];
    assert_eq!(bit_reversal_permutation(seq), expected);
}
