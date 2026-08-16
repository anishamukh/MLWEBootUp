load("constants.sage")
load("auxiliary_functions.sage")
load("KeyGen.sage")
load("Encrypt-Decrypt.sage")
load("HMult.sage")
load("Relin.sage")

# This example now works for arbitrary MLWE rank set in constants.sage.
# For the highest parameter set, keep POLDEG=2^14 and L=7.

key = public_private_keygen()
evk_key = relin_keygen(key)

m_plain1 = 0.4
m_plain2 = 0.6

ct1 = encryption(m_plain1, key)
ct2 = encryption(m_plain2, key)

print("m_plain1=", decryption(ct1, key, 0, 0))
print("m_plain2=", decryption(ct2, key, 0, 0))

ct_prod_nonlin = homomorphic_mult(ct1, ct2, L-1, key)
ct_prod = relinearisation(ct_prod_nonlin, evk_key)

print("HE_dec=", decryption(ct_prod, key, L-1, L-2))
print("original_dec=", m_plain1*m_plain2)
