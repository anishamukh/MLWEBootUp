# HEAANBOOT -> ModHE Sage port

This directory contains a Sage adaptation of the 2018 HEAANBOOT implementation to the ModHE/MLWE-RNS proof-of-concept used in this project.

## Files

- `bootstrap_chain.sage` — extends the ModHE RNS chain with bootstrap primes and recomputes `Q`, `p`, `B`, `q_hat`, etc. Run this **before key generation**.
- `bootstrap_heaanboot_modhe.sage` — bootstrap implementation.
- `bootstrap_example_modhe.sage` — end-to-end driver analogous to `TestScheme::testBootstrap()`.
- `boost_rankup_bootstrap_example.sage` — rank-1 q0 client encryption -> rank-up -> full bootstrap, matching the current client-Boost research path.

## C++ -> Sage mapping

| HEAANBOOT C++ | Sage port |
|---|---|
| `Context::addBootContext` | `heaan_add_boot_context` |
| `Scheme::addBootKey` | `heaan_boot_keygen` |
| `normalizeAndEqual` + setting large `logq` | `heaan_mod_raise_q0` |
| SubSum loop in `bootstrapAndEqual` | `heaan_subsum` |
| `coeffToSlotAndEqual` | `heaan_coeff_to_slot` |
| `exp2piAndEqual` | `heaan_exp2pi` |
| `evalExpAndEqual` | `heaan_eval_exp` |
| `slotToCoeffAndEqual` | `heaan_slot_to_coeff` |
| `bootstrapAndEqual` | `heaan_bootstrap` |
| `leftRotateFast` | `hb_left_rotate_fast` |
| `conjugate` | `hb_conjugate` |
| `reScaleBy(logp)` | one ModHE RNS rescale (`hb_rescale_one`) |
| `divByPo2(bits)` | centered CRT divide-and-reshare (`hb_div_by_po2`) |

## MLWE adaptation

A ModHE ciphertext has the form `(c0, c1-vector)` and decrypts as
`c0 + <c1,s>`.  Applying a ring automorphism to every ciphertext component changes the decryption secret from `s` to `phi(s)`.  Therefore rotation/conjugation are implemented as:

1. apply the automorphism component-wise;
2. key-switch `phi(s) -> s` using the existing generic `rankup_keygen` / `rank_up` machinery at the **same rank**.

This is the module analogue of HEAANBOOT's RLWE rotation and conjugation keys.

## RNS adaptation of the HEAAN bit modulus

The original repository uses a coefficient modulus `2^logq` and bit shifts. ModHE uses a product of odd RNS primes, so a literal representation-level port is impossible. The Sage code preserves the same arithmetic circuit with these replacements:

- `modDownTo` -> drop high RNS shares without scaling;
- `reScaleBy(logp)` -> consume one RNS prime;
- small non-RNS-aligned shifts -> centered CRT reconstruction, rounded division, and RNS re-sharing;
- C++ `cipher.logp = logq + 4` -> temporary scale `16*q0` because in the original `q = 2^logq`.

The bootstrap primes are generated close to `16*q0`. This keeps the scale after each RNS rescale close to the HEAANBOOT scale while allowing the original ModHE evaluation primes to remain below them and reappear after bootstrap.

## Usage

From the directory containing the original ModHE Sage files:

```text
sage bootstrap_example_modhe.sage
```

The default `LOGI=4` follows the HEAANBOOT implementation comment for secret Hamming weight 64. `LOGT` is parameter-dependent.

## Important validation note

The files are a research/proof-of-concept translation. The extended coefficient modulus changes the MLWE security instance, so the resulting `(POLDEG, MLWE rank, Q, error)` parameters must be re-estimated before being treated as secure. The CRT-based emulation of non-RNS-aligned bit shifts is intentionally slow and intended for correctness testing, not performance benchmarking.
