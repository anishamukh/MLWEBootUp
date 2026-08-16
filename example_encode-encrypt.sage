# -----------------------------------------------------------------------------
# example_ckks_codec_integrated.sage
#
# Minimal test for integrating true base-ring CKKS encoding/decoding into the
# revamped MLWE-CKKS / ModHE Sage prototype.
#
# For fast testing, temporarily set POLDEG small in constants.sage, e.g. 8, 16,
# or 32.  The direct encoder solves a dense POLDEG x POLDEG complex linear system
# and is not intended for POLDEG=2^14 performance.
# -----------------------------------------------------------------------------

load("constants.sage")
load("auxiliary_functions.sage")
load("KeyGen.sage")
load("Encrypt-Decrypt.sage")
load("encode.sage")
load("encode-encrypt.sage")

'''
print("POLDEG=", POLDEG, "slots=", ckks_num_slots(), "MLWE=", MLWE)

# 1. Pure codec test.
slots = [CCK(0.25), CCK(-0.5), CCK(0.125) + CCK(0.0625)*CCK(I)]
pt_poly = ckks_encode(slots, p_enc)
dec_slots = ckks_decode(pt_poly, p_enc, len(slots))
print("\n[codec only]")
for i in range(len(slots)):
    print(i, "in=", slots[i], "out=", dec_slots[i], "err=", abs(dec_slots[i]-slots[i]))

# 2. Full ModHE encrypt/decrypt/decode test.
key = public_private_keygen()
ct = encryption_ckks(slots, key, p_enc)
out_slots = decryption_ckks(ct, key, 0, 0, p_enc, len(slots))
print("\n[ModHE encrypt/decrypt + CKKS decode]")
for i in range(len(slots)):
    print(i, "in=", slots[i], "out=", out_slots[i], "err=", abs(out_slots[i]-slots[i]))

# 3. Compatibility with old constant-message style.
ct_const = encryption_ckks([0.4], key, p_enc)
print("\nconstant slot 0 ~=", decryption_ckks(ct_const, key, 0, 0, p_enc, 1)[0])
print("constant coefficient helper ~=", decryption_constant_compat(ct_const, key, 0, 0, p_enc))
'''

print("POLDEG=", POLDEG, "slots=", ckks_num_slots(), "MLWE=", MLWE)
print("CKKS codec method: fft")

# -----------------------------------------------------------------------------
# 1. Pure codec test: slots -> encode -> decode
# -----------------------------------------------------------------------------
slots = [CCK(0.25), CCK(-0.5), CCK(0.125) + CCK(0.0625)*CCK(I)]
pt_poly = ckks_encode(slots, p_enc, method="fft")
dec_slots = ckks_decode(pt_poly, p_enc, len(slots), method="fft")

print("\n[codec only: FFT/iFFT]")
for i in range(len(slots)):
    print(i, "in=", slots[i], "out=", dec_slots[i], "err=", abs(dec_slots[i]-slots[i]))

# Optional: compare against the slow O(N^2) matrix-vector reference for small N.
if POLDEG <= CKKS_MATRIX_REF_MAX_N:
    print("\n[FFT codec vs matrix-vector reference]")
    ckks_codec_compare_fft_matrix(slots, p_enc)
else:
    print("\n[FFT codec vs matrix-vector reference]")
    print("skipped because POLDEG is larger than", CKKS_MATRIX_REF_MAX_N)

# -----------------------------------------------------------------------------
# 2. Full ModHE encrypt/decrypt/decode test
# -----------------------------------------------------------------------------
key = public_private_keygen()
ct = encryption_ckks(slots, key, p_enc, method="fft")
out_slots = decryption_ckks(ct, key, 0, 0, p_enc, len(slots), method="fft")

print("\n[ModHE encrypt/decrypt + CKKS FFT decode]")
for i in range(len(slots)):
    print(i, "in=", slots[i], "out=", out_slots[i], "err=", abs(out_slots[i]-slots[i]))

# -----------------------------------------------------------------------------
# 3. Compatibility with old constant-message style
# -----------------------------------------------------------------------------
ct_const = encryption_ckks([0.4], key, p_enc, method="fft")
print("\nconstant slot 0 ~=", decryption_ckks(ct_const, key, 0, 0, p_enc, 1, method="fft")[0])
print("constant coefficient helper ~=", decryption_constant_compat(ct_const, key, 0, 0, p_enc))