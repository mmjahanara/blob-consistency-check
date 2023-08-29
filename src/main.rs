use std::env::set_var;

use halo2_ecc::fields::fp::FpChip;
use halo2_ecc::fields::{FieldChip, PrimeField};
use halo2_base::halo2_proofs::{
    dev::MockProver,
    halo2curves::bls12_381::Scalar,
    halo2curves::bn256::Fr
};

use halo2_base::gates::builder::{GateThreadBuilder, RangeCircuitBuilder};
use halo2_base::gates::RangeChip;
use halo2_base::Context;

const K: usize = 10;

fn fp_chip_test(
    k: usize,
    lookup_bits: usize,
    limb_bits: usize,
    num_limbs: usize,
    f: impl Fn(&mut Context<Fr>, &FpChip<Fr, Scalar>),
) {
    set_var("LOOKUP_BITS", lookup_bits.to_string());
    let range = RangeChip::<Fr>::default(lookup_bits);
    let chip = FpChip::<Fr, Scalar>::new(&range, limb_bits, num_limbs);

    let mut builder = GateThreadBuilder::mock();
    f(builder.main(0), &chip);

    builder.config(k, Some(10));
    let circuit = RangeCircuitBuilder::mock(builder);
    MockProver::run(k as u32, &circuit, vec![]).unwrap().assert_satisfied();
}


fn main() {
    let limb_bits = 88;
    let num_limbs = 3;
    fp_chip_test(K, K - 1, limb_bits, num_limbs, |ctx, chip| {
        let _a = Scalar::from(2);
        let _b = Scalar::from(7);

        let [a, b] = [_a, _b].map(|x| chip.load_private(ctx, x));
        let c = chip.mul(ctx, a.clone(), b.clone());

        // print log value of c
        print!("{} * {} = {}\n", a.value(), b.value(), c.value());

    });

}
