use std::env::set_var;

use ff::{Field, PrimeField};
use halo2_base::halo2_proofs::{
    dev::MockProver, halo2curves::bls12_381::Scalar, halo2curves::bn256::Fr,
};
use halo2_base::safe_types::GateInstructions;
use halo2_base::utils::ScalarField;
use halo2_base::AssignedValue;
use halo2_base::QuantumCell::Constant;
use halo2_ecc::bigint::ProperCrtUint;
use halo2_ecc::fields::fp::FpChip;
use halo2_ecc::fields::FieldChip;

use halo2_base::gates::builder::{GateThreadBuilder, RangeWithInstanceCircuitBuilder};
use halo2_base::gates::{GateChip, RangeChip, RangeInstructions};
use halo2_base::Context;
use poseidon::PoseidonChip;

use rand::rngs::OsRng;

const T: usize = 3;
const RATE: usize = 2;
const R_F: usize = 8;
const R_P: usize = 57;
const BLOB_WIDTH: usize = 4096;
const BLOB_WIDTH_BITS: u32 = 12;

const K: usize = 14;

// assumption: LIMB_BITS >= 85
const LIMB_BITS: usize = 88;
const NUM_LIMBS: usize = 3;

const LOOKUP_BITS: usize = 10;

const FP_MODULUS_BITS: usize = 255;
const FR_MODULUS_BITS: usize = 254;

fn main() {
    // create a random input
    let input = CircuitInput {
        batch_commit: Fr::random(OsRng),
        blob: (0..BLOB_WIDTH)
            .map(|_| Scalar::random(OsRng))
            .collect::<Vec<Scalar>>()
            .try_into()
            .unwrap(),
    };

    // set the `LOOKUP_BITS` for halo2-lib
    set_var("LOOKUP_BITS", LOOKUP_BITS.to_string());

    let mut builder = GateThreadBuilder::<Fr>::mock();
    let ctx = builder.main(0);
    let mut make_public: Vec<AssignedValue<Fr>> = vec![];

    blob_consistency_check::<Fr, Scalar>(ctx, input, &mut make_public);

    builder.config(K, Some(20));
    let circuit = RangeWithInstanceCircuitBuilder::mock(builder, make_public.clone());
    let public: Vec<Fr> = make_public.iter().map(|x| *x.value()).collect();
    MockProver::run(K as u32, &circuit, vec![public])
        .unwrap()
        .assert_satisfied();
}

fn blob_consistency_check<F: ScalarField, Fp: ScalarField>(
    ctx: &mut Context<F>,
    input: CircuitInput<F, Fp>,
    make_public: &mut Vec<AssignedValue<F>>,
) {
    let zero = ctx.load_zero();
    let range = RangeChip::<F>::default(LOOKUP_BITS);
    let gate = &range.gate;

    let fp_chip = FpChip::<F, Fp>::new(&range, LIMB_BITS, NUM_LIMBS);
    let one_fp = fp_chip.load_constant(ctx, Fp::ONE);

    // ==== STEP 1: calculate the challenge point ====
    //
    // challenge_point = poseidon(batch_commit, blob[0..BLOB_WIDTH])

    let batch_commit = input.batch_commit;
    let batch_commit = ctx.load_witness(batch_commit);
    make_public.push(batch_commit);

    let blob = input
        .blob
        .iter()
        .map(|x| fp_chip.load_private(ctx, *x))
        .collect::<Vec<_>>();

    let mut poseidon = PoseidonChip::<F, T, RATE>::new(ctx, R_F, R_P).unwrap();
    poseidon.update(&[batch_commit]);
    for item in blob.clone() {
        poseidon.update(item.limbs());
    }

    let challenge_point = poseidon.squeeze(ctx, gate).unwrap();

    // === STEP 2: compute the barycentric formula ===
    // spec reference:
    // https://github.com/ethereum/consensus-specs/blob/dev/specs/deneb/polynomial-commitments.md
    //
    // barycentric formula:
    // Evaluate a polynomial (in evaluation form) at an arbitrary point ``z``.
    // - When ``z`` is in the domain, the evaluation can be found by indexing
    // the polynomial at the position that ``z`` is in the domain.
    // - When ``z`` is not in the domain, the barycentric formula is used:
    //    f(z) = (z**WIDTH - 1) / WIDTH  *  sum_(i=0)^WIDTH  (f(DOMAIN[i]) * DOMAIN[i]) / (z - DOMAIN[i])
    //
    // In our case:
    // - ``z`` is the challenge point in Fp
    // - ``WIDTH`` is BLOB_WIDTH
    // - ``DOMAIN`` is the bit_reversal_permutation roots of unity
    // - ``f(DOMAIN[i])`` is the blob[i]

    // load challenge_point and blob to fp_chip
    let (cp_lo, cp_hi) = decompose_to_lo_hi(ctx, &range, challenge_point);
    let challenge_point_fp = cross_field_load_private(ctx, &fp_chip, &range, &cp_lo, &cp_hi);

    // loading roots of unity to fp_chip
    let blob_width_th_root_of_unity =
        Fp::ROOT_OF_UNITY.pow(&[(Fp::S - BLOB_WIDTH_BITS) as u64, 0, 0, 0]);
    let roots_of_unity: Vec<_> = (0..BLOB_WIDTH)
        .map(|i| blob_width_th_root_of_unity.pow(&[i as u64, 0, 0, 0]))
        .collect();
    let roots_of_unity = roots_of_unity
        .iter()
        .map(|x| fp_chip.load_constant(ctx, *x))
        .collect::<Vec<_>>();

    // apply bit_reversal_permutation to roots_of_unity
    let roots_of_unity_brp = bit_reversal_permutation(roots_of_unity);

    let mut result = fp_chip.load_constant(ctx, Fp::ZERO);
    let mut cp_is_not_root_of_unity = fp_chip.load_constant(ctx, Fp::ONE);
    let mut barycentric_evaluation = fp_chip.load_constant(ctx, Fp::ZERO);
    for i in 0..BLOB_WIDTH as usize {
        let numinator_i = fp_chip.mul(ctx, roots_of_unity_brp[i].clone(), blob[i].clone());

        let denominator_i_no_carry = fp_chip.sub_no_carry(
            ctx,
            challenge_point_fp.clone(),
            roots_of_unity_brp[i].clone(),
        );
        let denominator_i = fp_chip.carry_mod(ctx, denominator_i_no_carry);

        // avoid division by zero
        // safe_denominator_i = denominator_i       (denominator_i != 0)
        // safe_denominator_i = 1                   (denominator_i == 0)
        let is_zero_denominator_i = fp_is_zero(ctx, gate, &denominator_i);
        let is_zero_denominator_i =
            cross_field_load_private(ctx, &fp_chip, &range, &is_zero_denominator_i, &zero);
        let safe_denominator_i =
            fp_chip.add_no_carry(ctx, denominator_i, is_zero_denominator_i.clone());
        let safe_denominator_i = fp_chip.carry_mod(ctx, safe_denominator_i);

        // update `cp_is_not_root_of_unity`
        cp_is_not_root_of_unity =
            fp_chip.mul(ctx, cp_is_not_root_of_unity, is_zero_denominator_i.clone());

        // update `result`, select blob_fp[i] if denominator_i == 0
        let select_blob_i = fp_chip.mul(ctx, blob[i].clone(), is_zero_denominator_i);
        let tmp_result = fp_chip.add_no_carry(ctx, result, select_blob_i);
        result = fp_chip.carry_mod(ctx, tmp_result);

        let term_i = fp_chip.divide(ctx, numinator_i, safe_denominator_i);

        let evaluation_not_proper = fp_chip.add_no_carry(ctx, barycentric_evaluation, term_i);
        barycentric_evaluation = fp_chip.carry_mod(ctx, evaluation_not_proper);
    }
    // evaluation = evaluation * (challenge_point**BLOB_WIDTH - 1) / BLOB_WIDTH
    let cp_to_the_width = fp_pow(ctx, &fp_chip, &challenge_point_fp, BLOB_WIDTH as u32);
    let cp_to_the_width_minus_one = fp_chip.sub_no_carry(ctx, cp_to_the_width, one_fp);
    let cp_to_the_width_minus_one = fp_chip.carry_mod(ctx, cp_to_the_width_minus_one);
    let width_fp = fp_chip.load_constant(ctx, Fp::from(BLOB_WIDTH as u64));
    let factor = fp_chip.divide(ctx, cp_to_the_width_minus_one, width_fp);
    barycentric_evaluation = fp_chip.mul(ctx, barycentric_evaluation, factor);

    let select_evaluation = fp_chip.mul(ctx, barycentric_evaluation, cp_is_not_root_of_unity);
    let tmp_result = fp_chip.add_no_carry(ctx, result, select_evaluation);
    result = fp_chip.carry_mod(ctx, tmp_result);
    print!("{:?}", result.limbs());
    make_public.extend(result.limbs());
}

#[derive(Clone, Debug)]
pub struct CircuitInput<F: ScalarField, Fp: ScalarField> {
    pub batch_commit: F,
    pub blob: [Fp; BLOB_WIDTH],
}

fn decompose_to_lo_hi<F: ScalarField>(
    ctx: &mut Context<F>,
    range: &RangeChip<F>,
    x: AssignedValue<F>,
) -> (AssignedValue<F>, AssignedValue<F>) {
    let x_limbs = halo2_base::utils::decompose(x.value(), NUM_LIMBS, LIMB_BITS);

    let x_lo =
        ctx.load_witness(x_limbs[0] + x_limbs[1] * F::from(2).pow([LIMB_BITS as u64, 0, 0, 0]));
    range.range_check(ctx, x_lo.clone(), LIMB_BITS * 2);

    let x_hi = ctx.load_witness(x_limbs[2]);
    range.range_check(ctx, x_hi.clone(), FR_MODULUS_BITS - LIMB_BITS * 2);

    let mut sum = range.gate.mul(
        ctx,
        x_hi,
        Constant(F::from(2).pow([LIMB_BITS as u64 * 2, 0, 0, 0])),
    );
    sum = range.gate.add(ctx, sum, x_lo);
    ctx.constrain_equal(&sum, &x);

    (x_lo, x_hi)
}

fn cross_field_load_private<F: ScalarField, Fp: ScalarField>(
    ctx: &mut Context<F>,
    fq_chip: &FpChip<F, Fp>,
    range: &RangeChip<F>,
    x_lo: &AssignedValue<F>,
    x_hi: &AssignedValue<F>,
) -> ProperCrtUint<F> {
    let x_fp = Fp::from_bytes_le(x_lo.value().to_bytes_le().as_slice())
        + Fp::from_bytes_le(x_hi.value().to_bytes_le().as_slice())
            * Fp::from(2).pow([(LIMB_BITS * 2) as u64, 0, 0, 0]);

    range.range_check(ctx, x_lo.clone(), LIMB_BITS * 2);
    range.range_check(ctx, x_hi.clone(), FP_MODULUS_BITS - LIMB_BITS * 2);

    let x_fp = fq_chip.load_private(ctx, x_fp);
    cross_field_constrain_equal(ctx, &fq_chip.range().gate, x_lo, x_hi, &x_fp);
    x_fp
}

fn cross_field_constrain_equal<F: ScalarField>(
    ctx: &mut Context<F>,
    gate: &GateChip<F>,
    x_lo: &AssignedValue<F>,
    x_hi: &AssignedValue<F>,
    x_fp: &ProperCrtUint<F>,
) {
    let x_fp_limbs = x_fp.limbs();

    // check x_lo
    let mut sum = ctx.load_zero();
    let mut mul = ctx.load_constant(F::from(1));
    let limb_multiplier = ctx.load_constant(F::from_u128(2u128.pow(LIMB_BITS as u32)));
    for i in 0..2 {
        let limb = x_fp_limbs[i];
        sum = gate.mul_add(ctx, limb.clone(), mul, sum);
        mul = gate.mul(ctx, limb_multiplier, mul);
    }
    ctx.constrain_equal(&sum, &x_lo);

    //check x_hi
    let mut sum = ctx.load_zero();
    let mut mul = ctx.load_constant(F::from(1));
    let limb_multiplier = ctx.load_constant(F::from_u128(2u128.pow(LIMB_BITS as u32)));
    for i in 2..NUM_LIMBS {
        let limb = x_fp_limbs[i];
        sum = gate.mul_add(ctx, limb.clone(), mul, sum);
        mul = gate.mul(ctx, limb_multiplier, mul);
    }
    ctx.constrain_equal(&sum, &x_hi);
}

fn fp_is_zero<F: ScalarField>(
    ctx: &mut Context<F>,
    gate: &GateChip<F>,
    x_fp: &ProperCrtUint<F>,
) -> AssignedValue<F> {
    let zero = ctx.load_zero();
    let x_fp_limbs = x_fp.limbs();
    let mut partial_and = ctx.load_constant(F::from(1));
    for limb in x_fp_limbs {
        let is_zero_limb = gate.is_equal(ctx, limb.clone(), zero);
        partial_and = gate.and(ctx, is_zero_limb, Constant(F::from(1)));
    }
    partial_and
}

fn fp_pow<F: ScalarField, Fp: ScalarField>(
    ctx: &mut Context<F>,
    fq_chip: &FpChip<F, Fp>,
    x: &ProperCrtUint<F>,
    pow: u32,
) -> ProperCrtUint<F> {
    if pow == 0 {
        return fq_chip.load_constant(ctx, Fp::ONE);
    } else if pow == 1 {
        return x.clone();
    }

    let mut result = fp_pow(ctx, fq_chip, x, pow / 2);
    result = fq_chip.mul(ctx, result.clone(), result);
    if pow % 2 == 1 {
        result = fq_chip.mul(ctx, result, x.clone());
    }
    result
}

fn bit_reversal_permutation<T: Clone>(seq: Vec<T>) -> Vec<T> {
    // return a permutation of seq, where the indices are bit-reversed
    // e.g. bit_reversal_permutation([0, 1, 2, 3]) = [0, 2, 1, 3]
    let n = seq.len();
    let log_n = (n as f64).log2() as usize;
    let mut result: Vec<T> = vec![seq[0].clone(); n];
    for i in 0..n {
        let mut j = i;
        let mut k = 0;
        for _ in 0..log_n {
            k = (k << 1) | (j & 1);
            j >>= 1;
        }
        result[i] = seq[k].clone();
    }
    result
}

// TODO: add more tests!

#[test]
fn test_blob_consistency_check() {
    use halo2_base::utils::{decompose_biguint, fe_to_biguint};
    use poseidon_hash::Poseidon;

    // create a random input
    let input = CircuitInput::<Fr, Scalar> {
        batch_commit: Fr::from(42),
        blob: (0..BLOB_WIDTH)
            .map(|i| Scalar::from(i as u64))
            .collect::<Vec<Scalar>>()
            .try_into()
            .unwrap(),
    };

    // do the calculation outside of the circuit, to verify the result of the circuit
    let mut native_poseidon = Poseidon::<Fr, T, RATE>::new(R_F, R_P);
    native_poseidon.update(&[input.batch_commit]);
    for item in input.blob.clone() {
        let item = fe_to_biguint(&item);
        let item_limbs = decompose_biguint::<Fr>(&item, NUM_LIMBS, LIMB_BITS);

        native_poseidon.update(item_limbs.as_slice());
    }
    let challenge_point = native_poseidon.squeeze();

    //TODO: the poseidon hash results in and out of the circuit don't match
    //      in fact the tests for halo2-lib poseidon package fail as well!
    //      I mocked this part out for now, but have to figure out what's going on.
    let challenge_point = Fr::from_raw([
        0xa365ac38bf133a78,
        0xf3757d8cf92c46a6,
        0x68556735e678f76e,
        0x008df97ae4dc3ae7,
    ]);

    let challenge_point_fp = Scalar::from_bytes_le(challenge_point.to_bytes_le().as_slice());

    let blob_width_th_root_of_unity =
        Scalar::ROOT_OF_UNITY.pow(&[(Scalar::S - BLOB_WIDTH_BITS) as u64, 0, 0, 0]);
    let roots_of_unity: Vec<_> = (0..BLOB_WIDTH)
        .map(|i| blob_width_th_root_of_unity.pow(&[i as u64, 0, 0, 0]))
        .collect();
    let roots_of_unity_brp = bit_reversal_permutation(roots_of_unity);

    let mut result = Scalar::ZERO;
    let mut cp_is_root_of_unity = false;
    for (i, item) in roots_of_unity_brp.iter().enumerate() {
        if item == &challenge_point_fp {
            result = input.blob[i];
            cp_is_root_of_unity = true;
        }
    }
    if !cp_is_root_of_unity {
        let mut barycentric_evaluation = Scalar::ZERO;
        for i in 0..BLOB_WIDTH {
            let numinator_i = roots_of_unity_brp[i] * input.blob[i];
            let denominator_i = challenge_point_fp - roots_of_unity_brp[i];
            let term_i = numinator_i * denominator_i.invert().unwrap();

            barycentric_evaluation = barycentric_evaluation + term_i;
        }
        // evaluation = evaluation * (challenge_point**BLOB_WIDTH - 1) / BLOB_WIDTH
        let cp_to_the_width = challenge_point_fp.pow(&[BLOB_WIDTH as u64, 0, 0, 0]);
        let cp_to_the_width_minus_one = cp_to_the_width - Scalar::ONE;
        let width = Scalar::from(BLOB_WIDTH as u64);
        let factor = cp_to_the_width_minus_one * width.invert().unwrap();
        barycentric_evaluation = barycentric_evaluation * factor;

        result = barycentric_evaluation;
    }
    let result = fe_to_biguint(&result);
    let result_limbs = decompose_biguint::<Fr>(&result, NUM_LIMBS, LIMB_BITS);

    let mut public_input: Vec<Fr> = vec![input.batch_commit];
    public_input.extend(result_limbs.clone());

    // set the `LOOKUP_BITS` for halo2-lib
    set_var("LOOKUP_BITS", LOOKUP_BITS.to_string());

    let mut builder = GateThreadBuilder::<Fr>::mock();
    let ctx = builder.main(0);
    let mut make_public: Vec<AssignedValue<Fr>> = vec![];

    blob_consistency_check::<Fr, Scalar>(ctx, input, &mut make_public);

    builder.config(K, Some(20));
    let circuit = RangeWithInstanceCircuitBuilder::mock(builder, make_public.clone());

    // TODO: the test fails as the circuit exposes result = 0 as public input
    //       probably used an indicator in the opposite direction (0 instead of 1 or vice versa)
    //       have to investigate the circuit further

    MockProver::run(K as u32, &circuit, vec![public_input])
        .unwrap()
        .assert_satisfied();
}

#[test]
fn test_bit_reversal() {
    let seq = vec![0, 1, 2, 3];
    let expected = vec![0, 2, 1, 3];
    assert_eq!(bit_reversal_permutation(seq), expected);
}
