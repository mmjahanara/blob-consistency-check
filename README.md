# blob consistency check gadget

I have implemented a proof of concept for an in-circuit *blob consistency check gadget*. It can be used by zk-rollups that intend to use blob storage as a data availability solution. **This is just a research artifact, rather than production-grade code, and has to be treated as such**. 

On an M1 Macbook Pro (10 CPU cores, 16 GB RAM) proof generation takes 138.97 seconds. The circuit has 28,083,027 advice cells and 3,393,116 lookup advice cells.

**remark #1:** *Currently `halo2curves` does not have an implementation of `BLS12-381`, for that reason in the code we use `BN254::Fq`, another comparable non-native field, as the blob field .*

We provide more context in the following.

## What changes in the Rollup Process post EIP-4844?

### commit

**current design:** the rollup node submits the transaction batch in an L1 transaction as calldata. The rollup contract subsequently **calculates the commitment** to this data, which will serve as **public input** to the batch verification circuit.

**new design:** the commitment transaction will be a blob-carrying transaction. The same transaction batch will be the blob and instead of calculating the commitment in the contract, we just copy the versioned hash carried by the transaction to storage.

**remark #2:** each transaction batch includes a certain number of L1 messages, currently `commitBatch` function takes these transactions into account when calculating the commitment to the batch. In the new design we still need to calculate the hash of these transactions, let's denote it by `C_L1` and the final commitment to the whole batch will be `C_batch = Keccak(C_L1 | blob_versioned_hash)`
*Savings: (1) blob vs calldata, (2) almost no more on-chain commitment (Keccak) calculation.*

### zkEVM Public Input circuit

**current design:** we don’t provide the whole transaction batch data as public input to the circuits. The KZG commitment to public inputs should be calculated by the verifier, and it is prohibitively costly to do that for the whole transaction batch in the verifier contract. Instead, we just feed in the hash of that data, and provide the data as an advice column (witness). 

The **purpose of the PI circuit is to ensure the witness matches the public commitment**, basically opening the commitment for all the other sub-circuits. The PI circuit calculates the hash of raw public input using the Keccak table. This is done by calculating the RLC of data and one lookup, note that we still have to pay for this in the Keccak circuit.

**new design:** we use random point evaluation to check the advice column containing raw public input and the commitment to the blob are consistent. The challenge here is that we use BN254 for KZG commitment of circuits, but BLS12-381 for blobs. If we had the same curve for both commitments, we could have just checked the equality of commitments in the verifier contract. However this is not the case, so we do the following:

- **off-chain - (Fiat-Shamir):** calculate `z = poseidon(blob + C_batch)` and `y = p(z)`
  
- **in-contract - blob commitment random evaluation:** use point evaluation precompile to assert `point_evaluation_precompile(blob_versioned_hash, z, y)`.
  
- **in-circuit - consistency check:** 
*public input:* `C_batch`, `z`, and `y`  
*advice:* `blob`, `C_L1`, and `blob_versioned_hash`
0. constrain `C_batch = Keccak(C_L1 + blob_versioned_hash)`
1. constrain `z = poseidon(blob + C_batch)`. 
2. constrain `y = p(z)`. circuit uses the barycentric formula which costs 4096 * 2 non-native field multiplications and divisions to evaluate the Lagrange polynomial of `blob` denoted by `p(X)` on `z` over the scalar field of BLS12-381.

*Extra gas costs: invoking the point evaluation precompile*

*Extra prover cost: calculating commitment to `blob` in circuit to obtain `z`, and evaluating the Lagrange polynomial over a non-native field with 4096 * 2 operations.*

**remark #3:** *This PoC currently skips item (0), and takes `C_batch` directly as a public input.*
