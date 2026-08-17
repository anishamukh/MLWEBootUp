# -----------------------------------------------------------------------------
# bootstrap_example_modhe.sage
#
#   configure_heaanboot_chain() MUST run before key generation.
# -----------------------------------------------------------------------------

load("constants.sage")
load("bootstrap_heaanboot_chain.sage")

# HEAANBOOT parameters.
# Repository default: logI = 4 (for h=64).  logT is experiment/parameter-set
# dependent.  Start with logT=0 for sanity chexk.
LOGT = 5
LOGI = 4

# Append exactly the RNS primes consumed by the packed bootstrap.  The existing
# q0,q1,... chain remains as the post-bootstrap evaluation chain.
boot_cfg = configure_heaanboot_chain(logT=LOGT, logI=LOGI,
                                     extra_guard_primes=0,
                                     keep_existing_eval_primes=True,
                                     verbose=True)

# Load the ModHE routines only after the chain has been finalized.
load("auxiliary_functions.sage")
load("KeyGen_rank_up_down.sage")
load("Encrypt-Decrypt_rank_up_down.sage")
load("HMult.sage")
load("Relin_rank_up_down.sage")
load("Flexibility.sage")
load("encode.sage")
load("encode-encrypt.sage")
load("bootstrap_heaanboot_modhe.sage")

# Use the rank selected in constants.sage.  The bootstrap code itself is rank
# generic: rotations/conjugation use same-rank MLWE key switching.
key = public_private_keygen(rank=MLWE)
relin_key = relin_keygen(key)

# Full packing by default, as in the standard packed HEAANBOOT path.
SLOTS = POLDEG // 2
if SLOTS < 1 or (SLOTS & (SLOTS-1)) != 0:
    raise ValueError("POLDEG/2 must be a positive power of two")

# Context::addBootContext + Scheme::addBootKey analogues.
boot_context = heaan_add_boot_context(SLOTS, nominal_logp=boot_cfg["boot_logp"], verbose=True)
boot_keys = heaan_boot_keygen(key, SLOTS,
                              exact_cpp_keyset=False,
                              verbose=True)

print("\nbootstrapping keys=", boot_keys)

heaan_bootstrap_depth_report(logT=LOGT, logI=LOGI, packed=True)

# Keep test values modest; the 2018 bootstrap assumes the encoded plaintext is
# small compared with q0 so that the scaled sine approximates modular reduction.
base_vals = [
    0.035 + 0.010*I,
   -0.020 + 0.015*I,
    0.010 - 0.025*I,
   -0.030 - 0.005*I,
]
vals = [base_vals[i % len(base_vals)] for i in range(SLOTS)]

print("\n================================================")
print("TEST 1: HEAAN primitive rotations/conjugation")
print("================================================")

# IMPORTANT:
# Use an ordinary full-chain encryption here.
# Do NOT use ModRaise, since we want to test the primitive operations
# independently of the bootstrap.
raw_full = encryption_ckks(vals, key, scale=p_enc, method="fft")
ct_full = hb_wrap(
    raw_full,
    scale=p_enc,
    slots=SLOTS,
    is_complex=True
)

print("\nfull-chain decrypt:")
print(hb_decrypt_slots(ct_full, key))


# ------------------------------------------------------------
# Rotation by one slot
# ------------------------------------------------------------
ct_rot1 = hb_left_rotate_fast(ct_full, 1, boot_keys)
dec_rot1 = hb_decrypt_slots(ct_rot1, key)

print("\nrotation by 1:")
print("got     =", dec_rot1)
print("expected=", vals[1:] + vals[:1])


# ------------------------------------------------------------
# Conjugation
# ------------------------------------------------------------
ct_conj = hb_conjugate(ct_full, boot_keys)
dec_conj = hb_decrypt_slots(ct_conj, key)

expected_conj = [CC(v).conjugate() for v in vals]

print("\nconjugation:")
print("got     =", dec_conj)
print("expected=", expected_conj)


# ------------------------------------------------------------
# Multiplication by i
# ------------------------------------------------------------
ct_i = hb_imult(ct_full)
dec_i = hb_decrypt_slots(ct_i, key)

expected_i = [CC(I) * CC(v) for v in vals]

print("\nmultiply by i:")
print("got     =", dec_i)
print("expected=", expected_i)

print("\n================================================")
print("TEST 2: normalized C2S -> S2C")
print("================================================")

logNh = ZZ(POLDEG//2).nbits() - 1

# This is the normalization HEAANBOOT performs immediately
# before CoeffToSlot.
ct_norm = hb_div_by_po2(ct_full, logNh)

print("after /2^logNh:")
print("got     =", hb_decrypt_slots(ct_norm, key))
print("expected=", [v / SLOTS for v in vals])

ct_c2s = heaan_coeff_to_slot(
    ct_norm,
    boot_context,
    boot_keys
)

ct_back = heaan_slot_to_coeff(
    ct_c2s,
    boot_context,
    boot_keys
)

got = hb_decrypt_slots(ct_back, key)

print("\nafter normalized C2S -> S2C:")
print("got     =", got)
print("expected=", vals)

print("\nerrors:")
for i in range(SLOTS):
    print(
        i,
        abs(complex(got[i]) - complex(vals[i]))
    )

print("\n================================================")
print("TEST 3: square -> relinearize -> RNS rescale")
print("================================================")

ct_sq = hb_square_rescale(ct_full, relin_key)

got_sq = hb_decrypt_slots(ct_sq, key)
expected_sq = [CC(v)**2 for v in vals]

print("cipher  =", ct_sq)
print("got     =", got_sq)
print("expected=", expected_sq)

print("\nerrors:")
for i in range(SLOTS):
    print(
        i,
        abs(complex(got_sq[i]) - complex(expected_sq[i]))
    )

print("\n================================================")
print("TEST 3B: repeated square -> relin -> rescale")
print("================================================")

test_vals = [
    0.80 + 0.05*I,
   -0.70 + 0.10*I,
    0.75 - 0.08*I,
   -0.65 - 0.04*I,
]

raw_test = encryption_ckks(
    test_vals,
    key,
    scale=p_enc,
    method="fft"
)

ct_test = hb_wrap(
    raw_test,
    scale=p_enc,
    slots=SLOTS,
    is_complex=True
)

expected = [CC(v) for v in test_vals]

for r in range(5):

    ct_test = hb_square_rescale(
        ct_test,
        relin_key
    )

    expected = [v*v for v in expected]

    got = hb_decrypt_slots(
        ct_test,
        key
    )

    errs = [
        abs(
            complex(got[i])
            -
            complex(expected[i])
        )
        for i in range(SLOTS)
    ]

    print(
        "round =", r+1,
        "level =", ct_test.level(),
        "maxerr =", max(errs)
    )

print("\n================================================")
print("TEST 3C: multiplication starting below full level")
print("================================================")

# Start from the correct full-chain ciphertext.
ct_low = hb_mod_down_to_level(
    ct_full,
    ct_full.level() - 3
)

print("full level =", ct_full.level())
print("test level =", ct_low.level())

# Modulus-down alone must not change the plaintext.
dec_low = hb_decrypt_slots(ct_low, key)

print("\nafter ModDown only:")
print("got     =", dec_low)
print("expected=", vals)

errs_low = [
    abs(complex(dec_low[i]) - complex(vals[i]))
    for i in range(SLOTS)
]
print("max ModDown error =", max(errs_low))


# Now multiply at THIS reduced level.
ct_low_sq = hb_square_rescale(
    ct_low,
    relin_key
)

got = hb_decrypt_slots(ct_low_sq, key)
expected = [CC(v)**2 for v in vals]

print("\nafter square/relin/rescale at lower level:")
print("cipher  =", ct_low_sq)
print("got     =", got)
print("expected=", expected)

errs = [
    abs(complex(got[i]) - complex(expected[i]))
    for i in range(SLOTS)
]

print("\nerrors:")
for i in range(SLOTS):
    print(i, errs[i])

print("max error =", max(errs))

# q0-only client encryption.  This is the ModHE-RNS analogue of encrypting at
# the low modulus before bootstrapping raises the ciphertext modulus.
raw_q0 = client_encrypt_q0_ckks(vals, key, scale=p_enc, method="fft")
ct = hb_wrap(raw_q0, scale=p_enc, slots=SLOTS, is_complex=True)

print("\n=== Before bootstrap ===")
print("cipher:", ct)
before = hb_decrypt_slots(ct, key)
print("decrypt:", before)

# Full sequence:
#   normalize/ModRaise -> SubSum -> CoeffToSlot -> EvalExp -> SlotToCoeff.
ct_boot = heaan_bootstrap(ct, boot_context, boot_keys, relin_key,
                          logT=LOGT, logI=LOGI,
                          output_scale=p_enc,
                          active_raise=L,
                          verbose=True)

print("\n=== After bootstrap ===")
print("cipher:", ct_boot)
after = hb_decrypt_slots(ct_boot, key)
print("decrypt:", after)

print("\n=== Slot errors ===")
maxerr = 0.0
for i in range(SLOTS):
    err = abs(complex(after[i]) - complex(vals[i]))
    maxerr = max(maxerr, err)
    print("slot", i, "expected=", vals[i], "got=", after[i], "abs_err=", err)
print("max abs error =", maxerr)

print("\nOutput active RNS shares =", ct_boot.active())
print("Output level             =", ct_boot.level())
print("Expected post-boot chain to expose the original evaluation primes.")
