use std::env::set_var;

use blob_consistency_check::*;

use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use halo2_base::{
    gates::builder::{GateThreadBuilder, RangeWithInstanceCircuitBuilder},
    AssignedValue, halo2_proofs::{poly::kzg::{commitment::{ParamsKZG, KZGCommitmentScheme}, multiopen::ProverSHPLONK}, plonk::{keygen_vk, keygen_pk, create_proof as halo2_create_proof}, transcript::{Blake2bWrite, Challenge255, TranscriptWriterBuffer}, halo2curves::bn256::G1Affine},
};
use halo2_base::halo2_proofs::halo2curves::bn256::{Fq, Fr, Bn256};
use rand::rngs::OsRng;

const BENCH_SAMPLES: usize = 10;

fn create_proof(c: &mut Criterion) {
    // The BlobField should be changed to BLS12-381 Scalar Field.
    type CircuitField = Fr;
    type BlobField = Fq;

    println!("== Generating params and keys");

    // generate a default input
    let default_input = CircuitInput::default();

    // set the `LOOKUP_BITS` for halo2-lib
    set_var("LOOKUP_BITS", LOOKUP_BITS.to_string());

    // create circuit for keygen
    let mut builder = GateThreadBuilder::<CircuitField>::keygen();
    let ctx = builder.main(0);
    let mut make_public: Vec<AssignedValue<CircuitField>> = vec![];
    blob_consistency_check_gadget::<CircuitField, BlobField>(ctx, default_input, &mut make_public);
    builder.config(K, Some(20));
    let circuit = RangeWithInstanceCircuitBuilder::keygen(builder, make_public.clone());

    let params = ParamsKZG::<Bn256>::setup(K as u32, OsRng);
    let vk = keygen_vk(&params, &circuit).expect("vk should not fail");
    let pk = keygen_pk(&params, vk, &circuit).expect("pk should not fail");

    let break_points = circuit.break_points();

    println!("== Generating a random circuit input");
    let random_input = CircuitInput::<CircuitField, BlobField>::random();

    println!("== Preparing circuit for proof generation");
    let mut builder = GateThreadBuilder::<CircuitField>::prover();
    let ctx = builder.main(0);
    let mut make_public: Vec<AssignedValue<CircuitField>> = vec![];
    blob_consistency_check_gadget::<CircuitField, BlobField>(ctx, random_input, &mut make_public);
    let public: Vec<CircuitField> = make_public.iter().map(|x| *x.value()).collect();
    builder.config(K, Some(20));        
    let circuit = RangeWithInstanceCircuitBuilder::prover(builder, make_public.clone(), break_points.clone());

    let mut group = c.benchmark_group("benches");
    group.sample_size(BENCH_SAMPLES);
    group.bench_with_input(
        BenchmarkId::new("create_proof", K as u32),
        &(&params, &pk),
        |bencher, &(params, pk)| {
            bencher.iter(|| {
                let mut transcript = Blake2bWrite::<_, _, Challenge255<_>>::init(vec![]);
                halo2_create_proof::<
                    KZGCommitmentScheme<Bn256>,
                    ProverSHPLONK<'_, Bn256>,
                    Challenge255<_>,
                    _,
                    Blake2bWrite<Vec<u8>, G1Affine, _>,
                    _,
                >(params, pk, &[circuit.clone()], &[&[public.as_slice()]], OsRng, &mut transcript)
                .expect("prover should not fail");
            })
        },
    );
    group.finish();
}

criterion_group!(benches, create_proof);
criterion_main!(benches);