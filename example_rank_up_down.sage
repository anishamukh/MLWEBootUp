load("constants.sage")
load("auxiliary_functions.sage")
load("KeyGen_rank_up_down.sage")
load("Encrypt-Decrypt_rank_up_down.sage")
load("HMult.sage")
load("Relin_rank_up_down.sage")
load("Flexibility.sage")

# -----------------------------------------------------------------------------
# Basic rank-2 ModHE/MLWE-CKKS multiplication + relinearisation
# -----------------------------------------------------------------------------
# This example deliberately creates a rank-2 key even if MLWE in constants.sage
# is changed.  The updated routines read the rank from each key/ciphertext.

key_r2 = public_private_keygen(rank=2)
evk_r2 = relin_keygen(key_r2)

m_plain1 = 0.4
m_plain2 = 0.6

ct1_r2 = encryption(m_plain1, key_r2)
ct2_r2 = encryption(m_plain2, key_r2)

print("m_plain1 rank-2 dec=", decryption(ct1_r2, key_r2, 0, 0))
print("m_plain2 rank-2 dec=", decryption(ct2_r2, key_r2, 0, 0))

ct_prod_nonlin = homomorphic_mult(ct1_r2, ct2_r2, L-1, key_r2)
ct_prod_r2 = relinearisation(ct_prod_nonlin, evk_r2)

print("rank-2 HE mult dec=", decryption(ct_prod_r2, key_r2, L-1, L-2))
print("original product=", m_plain1*m_plain2)

# -----------------------------------------------------------------------------
# Flexibility demo: rank 2 -> rank 4 -> rank 2
# -----------------------------------------------------------------------------
# key_r4 shares the first two secret polynomials with key_r2.  This makes the
# rank-down step a genuine ModHE prefix rank reduction back to the original
# rank-2 secret.

key_r4 = extend_key_with_extra_secret(key_r2, 4)

# Rank-up: key-switch the rank-2 ciphertext under the rank-4 secret.
rup_key_2_to_4 = rankup_keygen(key_r2[0], key_r4[0])
ct1_r4 = rank_up(ct1_r2, rup_key_2_to_4)
print("after rank-up 2->4 dec=", decryption(ct1_r4, key_r4, 0, 0))

# Optional: perform a rank-4 homomorphic multiplication after rank-up.
# This validates that the lifted ciphertext is compatible with the rest of the
# generalized rank code.
evk_r4 = relin_keygen(key_r4)
ct2_r4 = rank_up(ct2_r2, rup_key_2_to_4)
ct_prod_r4_nonlin = homomorphic_mult(ct1_r4, ct2_r4, L-1, key_r4)
ct_prod_r4 = relinearisation(ct_prod_r4_nonlin, evk_r4)
print("rank-4 HE mult after rank-up dec=", decryption(ct_prod_r4, key_r4, L-1, L-2))

# Rank-down: reduce rank 4 back to the first two secret components.
rdown_key_4_to_2 = rankred_keygen(key_r4, 2)
ct1_r2_again = rank_red(ct1_r4, rdown_key_4_to_2)
print("after rank-down 4->2 dec=", decryption(ct1_r2_again, key_r2, 0, 0))

ct_prod_r2_again = rank_red(ct_prod_r4, rdown_key_4_to_2)
print("rank-4 product rank-down 4->2 dec=", decryption(ct_prod_r2_again, key_r2, L-1, L-2))
