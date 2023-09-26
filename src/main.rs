use blob_consistency_check::*;
use halo2_base::halo2_proofs::halo2curves::bn256::{Fq, Fr};

use std::env::set_var;
use halo2_base::halo2_proofs::dev::MockProver;
use halo2_base::{
    gates::builder::{GateThreadBuilder, RangeWithInstanceCircuitBuilder},
    AssignedValue,
};

fn main() {
    // The BlobField should be changed to BLS12-381 Scalar Field.
    type CircuitField = Fr;
    type BlobField = Fq;

    println!("== Generating a random batch commitment and blob");
    let input = CircuitInput::<CircuitField, BlobField>::random();

    // set the `LOOKUP_BITS` for halo2-lib
    set_var("LOOKUP_BITS", LOOKUP_BITS.to_string());

    println!("== Instantiating the blob consistency check circuit");
    let mut builder = GateThreadBuilder::<CircuitField>::mock();
    let ctx = builder.main(0);
    let mut make_public: Vec<AssignedValue<CircuitField>> = vec![];

    blob_consistency_check_gadget::<CircuitField, BlobField>(ctx, input, &mut make_public);

    builder.config(K, Some(20));
    let circuit = RangeWithInstanceCircuitBuilder::mock(builder, make_public.clone());
    let public: Vec<CircuitField> = make_public.iter().map(|x| *x.value()).collect();

    println!("== Running the mock prover");
    MockProver::run(K as u32, &circuit, vec![public])
        .unwrap()
        .assert_satisfied();

    println!("== Mock prover is satisfied!");
}
