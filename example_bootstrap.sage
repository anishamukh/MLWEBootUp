# -----------------------------------------------------------------------------
# example_bootstrap_almost.sage
#
# Almost-valid traditional MLWE-CKKS bootstrap scaffold test.
#
# This tests:
#
#   q0 ciphertext
#       -> ModRaise
#       -> true two-ciphertext CoeffToSlot
#       -> optional EvalMod scaled-sine prototype
#       -> true two-ciphertext SlotToCoeff
#
# For first debugging, use small POLDEG in constants.sage, e.g. 8, 16, or 32.
# -----------------------------------------------------------------------------

load("constants.sage")
load("auxiliary_functions.sage")
load("KeyGen.sage")
load("Encrypt-Decrypt.sage")
load("encode.sage")
load("encode-encrypt.sage")
load("HMult.sage")
load("Relin.sage")
load("bootstrap_keys.sage")
load("bootstrap_ops.sage")
load("linear_transform.sage")
load("bootstrap_more_ops.sage")

print("POLDEG=", POLDEG, "slots=", ckks_num_slots(), "MLWE=", MLWE)
print("Testing almost-valid traditional MLWE-CKKS bootstrap scaffold.")


# -----------------------------------------------------------------------------
# Helper: keep only the q0 share of a normal ciphertext.
# This emulates a bottom-level bootstrap input.
# -----------------------------------------------------------------------------

def restrict_ct_to_q0(ct):
    """Return q0-only version of a ModHE ciphertext."""
    ct = normalize_ct(ct)
    c0_q0 = [ct[0][0]]
    c1_q0 = []
    for i in range(len(ct[1])):
        c1_q0.append([ct[1][i][0]])
    return [c0_q0, c1_q0]


# -----------------------------------------------------------------------------
# 1. Generate true bootstrap C2S/S2C diagonals.
# -----------------------------------------------------------------------------

print("\n[generating true bootstrap CoeffToSlot diagonals]")
c2s_boot_diags = generate_coeff_to_slot_bootstrap_diagonals(
    diag_scale=p_enc,
    factor=1,
    method="fft"
)

print("C2S lo A=", len(c2s_boot_diags["lo"]["A"]),
      "lo B=", len(c2s_boot_diags["lo"]["B"]))
print("C2S hi A=", len(c2s_boot_diags["hi"]["A"]),
      "hi B=", len(c2s_boot_diags["hi"]["B"]))


print("\n[generating true bootstrap SlotToCoeff diagonals]")
s2c_boot_diags = generate_slot_to_coeff_bootstrap_diagonals(
    diag_scale=p_enc,
    factor=1,
    method="fft"
)

print("S2C lo A=", len(s2c_boot_diags["lo"]["A"]),
      "lo B=", len(s2c_boot_diags["lo"]["B"]))
print("S2C hi A=", len(s2c_boot_diags["hi"]["A"]),
      "hi B=", len(s2c_boot_diags["hi"]["B"]))


# -----------------------------------------------------------------------------
# 2. Key generation.
# -----------------------------------------------------------------------------

print("\n[key generation]")
key = public_private_keygen()

# Option A: if you added boot_keygen_from_diagonal_bundles() to bootstrap_keys.sage.
try:
    boot_keys = boot_keygen_from_diagonal_bundles(
        key,
        c2s_boot_diags,
        s2c_boot_diags,
        include_relin=True
    )
    print("boot keys generated using boot_keygen_from_diagonal_bundles()")

# Option B: fallback using linear_transform.sage helper.
except NameError:
    rot_exps = merge_rotation_exponents_deep(c2s_boot_diags, s2c_boot_diags)
    print("rotation exponents needed=", len(rot_exps))
    boot_keys = boot_keygen(
        key,
        rotation_exponents=rot_exps,
        include_relin=True
    )
    print("boot keys generated using boot_keygen() fallback")


# -----------------------------------------------------------------------------
# 3. Encrypt normal CKKS slots, then restrict to q0.
# -----------------------------------------------------------------------------

slots = [
    CCK(0.25),
    CCK(-0.5) + CCK(0.125)*CCK(I),
    CCK(0.0625) - CCK(0.2)*CCK(I)
]

print("\n[encrypt original CKKS slots]")
ct_full = encryption_ckks(slots, key, p_enc, method="fft")

print("decrypt full ciphertext:")
orig_full = decryption_ckks(
    ct_full,
    key,
    0,
    0,
    p_enc,
    len(slots),
    method="fft"
)

for i in range(len(slots)):
    print(i, "in=", slots[i], "dec=", orig_full[i], "err=", abs(orig_full[i]-slots[i]))


ct_q0 = restrict_ct_to_q0(ct_full)

print("\n[decrypt q0-only bootstrap input]")
orig_q0 = decryption_ckks(
    ct_q0,
    key,
    0,
    0,
    p_enc,
    len(slots),
    method="fft"
)

for i in range(len(slots)):
    print(i, "in=", slots[i], "q0dec=", orig_q0[i], "err=", abs(orig_q0[i]-slots[i]))

# -----------------------------------------------------------------------------
# 4. First test: bootstrap scaffold without EvalMod.
#
# This tests:
#
#     q0 ciphertext
#       -> ModRaise
#       -> CoeffToSlot
#       -> SlotToCoeff
#
# No EvalMod is applied here.
# Also no normalization is applied here.
# -----------------------------------------------------------------------------

print("\n[bootstrap_almost_ckks test: run_evalmod=False]")

boot_params_no_evalmod = {
    "c2s_boot_diags": c2s_boot_diags,
    "s2c_boot_diags": s2c_boot_diags,

    "do_modraise": True,
    "old_basis": [q_list[0]],
    "new_basis": [q_list[j] for j in range(L)],
    "centered_modraise": True,

    "divide_by_diag_scale": True,

    # Important:
    # If EvalMod is disabled, do not normalize before C2S.
    "normalize_before_c2s": False,

    "run_evalmod": False,
    "method": "fft"
}

ct_no_evalmod = bootstrap_almost_ckks(
    ct_q0,
    boot_keys,
    boot_params_no_evalmod
)

out_no_evalmod = decryption_ckks(
    ct_no_evalmod,
    key,
    0,
    0,
    p_enc,
    len(slots),
    method="fft"
)

print("\n[result without EvalMod]")
max_err = 0
for i in range(len(slots)):
    err = abs(out_no_evalmod[i] - slots[i])
    max_err = max(max_err, err)
    print(i, "in=", slots[i], "out=", out_no_evalmod[i], "err=", err)

print("max err without EvalMod=", max_err)

"""
# -----------------------------------------------------------------------------
# 4. First test: bootstrap scaffold without EvalMod.
#
# This should approximately preserve the input, because it only tests:
#
#     ModRaise -> C2S -> S2C
#
# -----------------------------------------------------------------------------

print("\n[bootstrap_almost_ckks test: run_evalmod=False]")

boot_params_no_evalmod = {
    "c2s_boot_diags": c2s_boot_diags,
    "s2c_boot_diags": s2c_boot_diags,

    # q0-only input -> full RNS basis
    "old_basis": [q_list[0]],
    "new_basis": [q_list[j] for j in range(L)],
    "centered_modraise": True,

    # Prototype scale correction for diagonal plaintexts.
    "divide_by_diag_scale": True,

    # Disable EvalMod first.
    "run_evalmod": False,

    "method": "fft"
}

ct_no_evalmod = bootstrap_almost_ckks(
    ct_q0,
    boot_keys,
    boot_params_no_evalmod
)

out_no_evalmod = decryption_ckks(
    ct_no_evalmod,
    key,
    0,
    0,
    p_enc,
    len(slots),
    method="fft"
)

print("\n[result without EvalMod]")
max_err = 0
for i in range(len(slots)):
    err = abs(out_no_evalmod[i] - slots[i])
    max_err = max(max_err, err)
    print(i, "in=", slots[i], "out=", out_no_evalmod[i], "err=", err)

print("max err without EvalMod=", max_err)
"""
# -----------------------------------------------------------------------------
# 5. Second test: bootstrap scaffold with OpenFHE-style normalized EvalMod.
#
# This tests:
#
#     q0 ciphertext
#       -> ModRaise
#       -> normalize by p_enc/(K*q0)
#       -> CoeffToSlot
#       -> normalized EvalMod prototype
#       -> SlotToCoeff
# -----------------------------------------------------------------------------

print("\n[bootstrap_almost_ckks test: run_evalmod=True]")

boot_params_with_evalmod = {
    "c2s_boot_diags": c2s_boot_diags,
    "s2c_boot_diags": s2c_boot_diags,

    "do_modraise": True,
    "old_basis": [q_list[0]],
    "new_basis": [q_list[j] for j in range(L)],
    "centered_modraise": True,

    "divide_by_diag_scale": True,

    # OpenFHE-style normalization before CoeffToSlot.
    "normalize_before_c2s": True,

    # This tells eval_exp_or_sine_ckks_prototype()
    # to call eval_mod_sine_openfhe_prototype().
    "normalized_evalmod": True,

    "q0": q_list[0],
    "K": 1,
    "post": RR(q_list[0]) / RR(p_enc),

    "eval_degree": 7,
    "const_scale": p_enc,
    "method": "fft",

    # Keep this false for the prototype.
    "do_rescale": False,

    "run_evalmod": True
}

ct_boot = bootstrap_almost_ckks(
    ct_q0,
    boot_keys,
    boot_params_with_evalmod
)

out_boot = decryption_ckks(
    ct_boot,
    key,
    0,
    0,
    p_enc,
    len(slots),
    method="fft"
)

print("\n[result with OpenFHE-style EvalMod prototype]")
max_err = 0
for i in range(len(slots)):
    err = abs(out_boot[i] - slots[i])
    max_err = max(max_err, err)
    print(i, "in=", slots[i], "out=", out_boot[i], "err=", err)

print("max err with OpenFHE-style EvalMod prototype=", max_err)

"""
# -----------------------------------------------------------------------------
# 5. Second test: bootstrap scaffold with EvalMod prototype.
#
# This is NOT expected to be as clean yet, because scale/level management and
# EvalMod approximation are still prototype-level.
# -----------------------------------------------------------------------------

print("\n[bootstrap_almost_ckks test: run_evalmod=True]")

boot_params_with_evalmod = {
    "c2s_boot_diags": c2s_boot_diags,
    "s2c_boot_diags": s2c_boot_diags,

    "old_basis": [q_list[0]],
    "new_basis": [q_list[j] for j in range(L)],
    "centered_modraise": True,

    "divide_by_diag_scale": True,

    "q0": q_list[0],
    "eval_degree": 5,
    "const_scale": p_enc,
    "method": "fft",

    # Still prototype. Keep false first.
    "do_rescale": False,

    # Enable scaled-sine EvalMod prototype.
    "run_evalmod": True
}

ct_boot = bootstrap_almost_ckks(
    ct_q0,
    boot_keys,
    boot_params_with_evalmod
)
"""
out_boot = decryption_ckks(
    ct_boot,
    key,
    0,
    0,
    p_enc,
    len(slots),
    method="fft"
)

print("\n[result with EvalMod prototype]")
max_err = 0
for i in range(len(slots)):
    err = abs(out_boot[i] - slots[i])
    max_err = max(max_err, err)
    print(i, "in=", slots[i], "out=", out_boot[i], "err=", err)

print("max err with EvalMod prototype=", max_err)


print("\nNOTE:")
print("The run_evalmod=False test should be the first correctness target.")
print("The run_evalmod=True test wires the EvalMod arithmetic path, but it is")
print("still not production CKKS bootstrapping because scale/rescale/level")
print("management and the EvalMod approximation are only prototype-level.")

eval_mod_double_angle_selftest(logT=1, base_degree=3)