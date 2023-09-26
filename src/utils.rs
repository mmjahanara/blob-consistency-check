use halo2_base::QuantumCell::Constant;
use halo2_base::safe_types::{RangeInstructions, GateInstructions};
use halo2_base::{
    utils::ScalarField,
    Context, 
    safe_types::RangeChip, 
    AssignedValue, gates::GateChip,
};
use halo2_ecc::bigint::ProperCrtUint;
use halo2_ecc::fields::FieldChip;
use halo2_ecc::fields::fp::FpChip;


// assumption: LIMB_BITS >= 85
pub const LIMB_BITS: usize = 88;
pub const NUM_LIMBS: usize = 3;

// update this when FP is changed, e.g. 255 for BLS12-381 Scalar Field
const FP_MODULUS_BITS: usize = 254;
const FR_MODULUS_BITS: usize = 254;

// power of the largest power of two root of unity in Fp
// For BLS12-381, S = 32
//
// For BN254::Fq, S = 1, however, we need it to be higher than BLOB_WIDTH_BITS
// so we just set it to S = 32 for the test purposes. 
pub const FP_S: u32 = 32;


pub fn decompose_to_lo_hi<F: ScalarField>(
    ctx: &mut Context<F>,
    range: &RangeChip<F>,
    x: AssignedValue<F>,
) -> (AssignedValue<F>, AssignedValue<F>) {
    let x_limbs = halo2_base::utils::decompose(x.value(), NUM_LIMBS, LIMB_BITS);

    let x_lo =
        ctx.load_witness(x_limbs[0] + x_limbs[1] * F::from(2).pow(&[LIMB_BITS as u64, 0, 0, 0]));
    range.range_check(ctx, x_lo.clone(), LIMB_BITS * 2);

    let x_hi = ctx.load_witness(x_limbs[2]);
    range.range_check(ctx, x_hi.clone(), FR_MODULUS_BITS - LIMB_BITS * 2);

    let mut sum = range.gate.mul(
        ctx,
        x_hi,
        Constant(F::from(2).pow(&[LIMB_BITS as u64 * 2, 0, 0, 0])),
    );
    sum = range.gate.add(ctx, sum, x_lo);
    ctx.constrain_equal(&sum, &x);

    (x_lo, x_hi)
}

/*
given two AssignedValue<F>, x_lo and x_hi, in the native field F,
returns a ProperCrtUint<F> in the target field Fp (which is bigger)
*/
pub fn cross_field_load_private<F: ScalarField, Fp: ScalarField>(
    ctx: &mut Context<F>,
    fq_chip: &FpChip<F, Fp>,
    range: &RangeChip<F>,
    x_lo: &AssignedValue<F>,
    x_hi: &AssignedValue<F>,
) -> ProperCrtUint<F> {
    let x_fp = Fp::from_bytes_le(x_lo.value().to_bytes_le().as_slice())
        + Fp::from_bytes_le(x_hi.value().to_bytes_le().as_slice())
            * Fp::from(2).pow(&[(LIMB_BITS * 2) as u64, 0, 0, 0]);

    range.range_check(ctx, x_lo.clone(), LIMB_BITS * 2);
    range.range_check(ctx, x_hi.clone(), FP_MODULUS_BITS - LIMB_BITS * 2);

    let x_fp = fq_chip.load_private(ctx, x_fp);
    cross_field_constrain_equal(ctx, &fq_chip.range().gate, x_lo, x_hi, &x_fp);
    x_fp
}

/*
given x_fp, a ProperCrtUint<Fp> in the target field Fp,
and its decomposition x_lo and x_hi in the native field F,
constrains x_lo and x_hi to be equal to the decomposition of x_fp
*/
pub fn cross_field_constrain_equal<F: ScalarField>(
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

/*
given x_fp, a ProperCrtUint<Fp> in the target field Fp,
returns an AssignedValue 1 if x_fp is zero, and 0 otherwise.
*/
pub fn fp_is_zero<F: ScalarField>(
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

/*
raises x in Fp to the power of pow,
notice that pow is a constant
*/
pub fn fp_pow<F: ScalarField, Fp: ScalarField>(
    ctx: &mut Context<F>,
    fp_chip: &FpChip<F, Fp>,
    x: &ProperCrtUint<F>,
    pow: u32,
) -> ProperCrtUint<F> {
    if pow == 0 {
        return fp_chip.load_constant(ctx, Fp::one());
    } else if pow == 1 {
        return x.clone();
    }

    let mut result = fp_pow(ctx, fp_chip, x, pow / 2);
    result = fp_chip.mul(ctx, result.clone(), result);
    if pow % 2 == 1 {
        result = fp_chip.mul(ctx, result, x.clone());
    }
    result
}


/*
returns a clone of the input vector with indices bit-reversed
*/
pub fn bit_reversal_permutation<T: Clone>(seq: Vec<T>) -> Vec<T> {
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

#[test]
fn test_bit_reversal() {
    let seq = vec![0, 1, 2, 3];
    let expected = vec![0, 2, 1, 3];
    assert_eq!(bit_reversal_permutation(seq), expected);
}
