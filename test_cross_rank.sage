# =============================================================================
# crossrank_primitive_error_bench.sage
#
# Primitive/noise experiment for Section 5 (cross-rank C2S optimization).
#
# We compare, on the SAME rank-1/q0 client ciphertext:
#
#   Baseline local primitive at the actual C2S entry:
#       ct_C --RankUp_{1->r}@q0--> ct_up
#            --ModRaise + full-pack normalization--> C2S input
#            --Rot_{r->r,j}-->
#
#   Optimized local primitive at the actual C2S entry:
#       ct_C --ModRaise + full-pack normalization--> C2S input
#            --XRotUp_{1->r,j}-->
#
# The experiment reports three *added-error* quantities:
#
#   eps_RU       = ||e_RU||_can,inf / Delta_client
#   eps_Rot      = ||e_Rot(r->r)||_can,inf / Delta_boot
#   eps_XRot     = ||e_XRot(1->r)||_can,inf / Delta_boot
#
# where Delta_client is the ordinary CKKS client scale and Delta_boot=16*q0 is
# the scale used by the HEAANBOOT port at the C2S entry.  This stage-aware
# normalization is important: eps_RU, eps_Rot and eps_XRot must NOT simply be
# added as if they were errors at the same scale/stage.
#
# Thus this primitive experiment answers the local question:
#   "Does a 1->r cross-rank baby rotation itself introduce more switching
#    error than the ordinary r->r baby rotation it replaces?"
#
# The stronger pipeline question -- whether replacing RankUp+bootstrap-C2S by
# XRotUp-C2S changes final post-bootstrap precision -- is deliberately left to
# the subsequent paired end-to-end bootstrap experiment.  These local terms
# explain that result, but are not by themselves an exact additive prediction
# of it because RankUp noise is propagated through ModRaise/SubSum/C2S.
#
# IMPORTANT IMPLEMENTATION NOTE
# -----------------------------
# The existing proof-of-concept Flexibility.sage uses a generic ModDown whose
# representative convention is convenient for correctness testing but is not
# the nearest/centered ModDown used in the paper's RankUp error derivation.
# For a clean numerical experiment, this bench instruments the *paper-aligned*
# nearest-ModDown generalized RNS key-switch error directly.
#
# For a generalized source-rank r_src -> target-rank r_tgt key switch with one
# RNS digit per q_j active at the C2S entry, the exact added decryption error is
#
#   e_KS = P^{-1} sum_{u,j} c_{u,j} e_{u,j}
#          + rho_0 + <rho_1, s_tgt>,
#
# where rho_0,rho_1 are the nearest-ModDown rounding errors.  To obtain these
# rounding errors exactly, only the P-residues of the uniform key polynomials
# are required.  The omitted q_i residues do not affect e_KS.  Therefore this
# bench samples the exact marginal distribution relevant to the decryption
# error without materializing the enormous full P*Q switching keys.  This is
# essential for making n=2^12 (and exploratory larger n) practical in Sage.
#
# Parameter intent
# ----------------
#   q0       : the project's 60-bit client modulus.
#   P_RU     : a distinct NTT-friendly 60-bit prime >= q0 (d=1 RankUp).
#   q_i boot : generated near 16*q0, as in the current HEAANBOOT port.
#   P_BOOT   : an NTT-friendly prime strictly larger than every active q_i,
#              so the single-prime gadget digits satisfy P_BOOT >= max q_i.
#
# These are FUNCTIONAL/NUMERICAL experiment parameters, not certified security
# parameters.  Concrete security must be estimated separately at P_BOOT*Qboot.
#
# Typical runs:
#   CR_LOGN=12 CR_RANK=2 sage crossrank_primitive_error_bench.sage
#   CR_LOGN=12 CR_RANK=4 sage crossrank_primitive_error_bench.sage
#
# A first exploratory n=2^15 run should use very few trials, e.g.
#   CR_LOGN=15 CR_RANK=2 CR_KEY_TRIALS=1 CR_CT_TRIALS=1 \
#     sage crossrank_primitive_error_bench.sage
#
# By default we test one representative non-zero C2S baby rotation j=1.
# Additional baby indices may be supplied, e.g. CR_ROTATIONS="1,31".
# =============================================================================

import os
import csv
import random as _py_random

# -----------------------------------------------------------------------------
# 0. Base project parameters, ring override, and HEAANBOOT-style Qboot chain
# -----------------------------------------------------------------------------

load("constants.sage")

EXP_LOGN = int(os.environ.get("CR_LOGN", "12"))
EXP_RANK = int(os.environ.get("CR_RANK", "2"))
POLDEG = ZZ(2)**ZZ(EXP_LOGN)
MLWE = ZZ(EXP_RANK)

R.<x> = QQ[]
firreducible = x**POLDEG + 1

# Build the same type of temporary modulus chain used by the working bootstrap.
# LOGT/LOGI only determine how many temporary primes Qboot contains here; this
# bench does not evaluate EvalMod itself.
load("bootstrap_heaanboot_chain.sage")

LOGT = int(os.environ.get("CR_LOGT", "5"))
LOGI = int(os.environ.get("CR_LOGI", "4"))

boot_cfg = configure_heaanboot_chain(
    logT=LOGT,
    logI=LOGI,
    extra_guard_primes=0,
    keep_existing_eval_primes=True,
    verbose=False
)

# Helpers must be loaded only after POLDEG/q_list have their experimental values.
load("auxiliary_functions.sage")
load("encode.sage")
# Only public ModRaise + full-pack normalization are used from this file.
# The local key-switch errors themselves are instrumented below.
load("bootstrap_heaanboot_modhe.sage")

from sage.stats.distributions.discrete_gaussian_integer import \
    DiscreteGaussianDistributionIntegerSampler

R_CLIENT = 1
R_EVAL = int(MLWE)
SLOTS = int(POLDEG // 2)
DELTA = ZZ(p_enc)
Q0 = ZZ(q_list[0])

# Optional cap is useful only as a scalability diagnostic.  For the paper's
# actual Qboot primitive experiment leave CR_ACTIVE_Q unset so every Qboot digit
# is included.
ACTIVE_Q = int(os.environ.get("CR_ACTIVE_Q", str(len(q_list))))
if ACTIVE_Q < 2 or ACTIVE_Q > len(q_list):
    raise ValueError("CR_ACTIVE_Q must lie in [2, len(q_list)]")
ACTIVE_QS = [ZZ(q_list[j]) for j in range(ACTIVE_Q)]


def _ntt_prime_strictly_above(lower, used=None):
    """First prime p > lower with p = 1 mod 2*POLDEG."""
    if used is None:
        used = set()
    lower = ZZ(lower)
    step = ZZ(2*POLDEG)
    k = lower // step
    while True:
        k += 1
        cand = k*step + 1
        if cand not in used and cand.is_prime():
            return cand


# Nominal q0-only RankUp auxiliary modulus: same word size as q0, distinct,
# NTT-friendly for the actual tested degree, and >= q0.
if "CR_P_RU" in os.environ:
    P_RU = ZZ(os.environ["CR_P_RU"])
else:
    P_RU = _ntt_prime_strictly_above(Q0, used=set(ACTIVE_QS))

# High-modulus rotation/cross-rank switching modulus.  We require the *value*,
# not merely the bit-size, to exceed every active RNS digit.
if "CR_P_BOOT" in os.environ:
    P_BOOT = ZZ(os.environ["CR_P_BOOT"])
else:
    P_BOOT = _ntt_prime_strictly_above(
        max(ACTIVE_QS), used=set(ACTIVE_QS + [P_RU])
    )

# Keep global special-modulus helpers consistent in case a loaded project helper
# touches them.  Q/q_hat/B do not change because only p is replaced here.
p = ZZ(P_BOOT)
moduli = list(q_list) + [p]
p_inv_q = [ZZ(p).inverse_mod(ZZ(qi)) for qi in q_list]

QRAISE = ZZ(1)
for qi in ACTIVE_QS:
    QRAISE *= qi

# The fully-packed SubSum normalization in the current RNS port consumes one
# top q-prime before CoeffToSlot.  Therefore C2S key switching sees one fewer
# active digit than ModRaise itself.
C2S_ACTIVE_Q = ACTIVE_Q - 1
C2S_QS = ACTIVE_QS[:C2S_ACTIVE_Q]
QC2S = ZZ(1)
for qi in C2S_QS:
    QC2S *= qi
BOOT_SCALE = ZZ(16) * Q0

# -----------------------------------------------------------------------------
# 1. Concrete distributions and statistical configuration
# -----------------------------------------------------------------------------

SECRET_H = int(os.environ.get("CR_SECRET_H", "64"))
SIGMA = RR(os.environ.get("CR_SIGMA", "3.2"))
ZO_RHO = RR(os.environ.get("CR_ZO_RHO", "0.5"))

KEY_TRIALS = int(os.environ.get("CR_KEY_TRIALS", "3"))
CT_TRIALS = int(os.environ.get("CR_CT_TRIALS", "5"))
SEED = int(os.environ.get("CR_SEED", "20260815"))

# C2S baby-step size from the current HEAANBOOT BSGS split.
logSlots = ZZ(SLOTS).nbits() - 1
logk = logSlots // 2
BABY_K = 1 << int(logk)

rot_env = os.environ.get("CR_ROTATIONS", "1")
ROTATIONS = [int(v.strip()) for v in rot_env.split(",") if v.strip() != ""]
if len(ROTATIONS) == 0:
    raise ValueError("CR_ROTATIONS must contain at least one baby rotation")
for rot in ROTATIONS:
    if rot <= 0 or rot >= BABY_K:
        raise ValueError(
            "rotation %d is not a non-zero baby index; require 1 <= j < %d" %
            (rot, BABY_K)
        )

OUTPUT_CSV = os.environ.get(
    "CR_OUTPUT_CSV",
    "crossrank_primitive_trials_n2^%d_r%d.csv" % (EXP_LOGN, EXP_RANK)
)
SUMMARY_CSV = os.environ.get(
    "CR_SUMMARY_CSV", "crossrank_primitive_paper_summary.csv"
)

if R_EVAL < 2:
    raise ValueError("Use r >= 2 for the cross-rank experiment")
if SECRET_H > POLDEG:
    raise ValueError("secret Hamming weight exceeds POLDEG")
if not Q0.is_prime():
    raise ValueError("q0 is not prime")
if not P_RU.is_prime() or not P_BOOT.is_prime():
    raise ValueError("P_RU/P_BOOT must be prime")
if P_RU < Q0:
    raise ValueError("nominal q0 RankUp requires P_RU >= q0")
if P_RU == Q0:
    raise ValueError("P_RU must be distinct from q0")
if P_BOOT <= max(ACTIVE_QS):
    raise ValueError("require P_BOOT > max active q_i")
if (Q0 - 1) % (2*POLDEG) != 0:
    raise ValueError("q0 is not 1 mod 2n")
if (P_RU - 1) % (2*POLDEG) != 0:
    raise ValueError("P_RU is not 1 mod 2n")
if (P_BOOT - 1) % (2*POLDEG) != 0:
    raise ValueError("P_BOOT is not 1 mod 2n")
for qi in ACTIVE_QS:
    if not qi.is_prime() or (qi - 1) % (2*POLDEG) != 0:
        raise ValueError("an active q_i is not an NTT-friendly prime")

set_random_seed(int(SEED))
_py_random.seed(int(SEED))
_DG = DiscreteGaussianDistributionIntegerSampler(sigma=SIGMA)

# -----------------------------------------------------------------------------
# 2. Concrete samplers (same distributions as the RankUp paper bench)
# -----------------------------------------------------------------------------


def sample_hwt_secret_poly(h=SECRET_H):
    positions = _py_random.sample(range(int(POLDEG)), int(h))
    coeffs = [ZZ(0)] * int(POLDEG)
    for i in positions:
        coeffs[i] = ZZ(1) if _py_random.getrandbits(1) else ZZ(-1)
    return R(coeffs)


def sample_error_poly():
    return R([ZZ(_DG()) for _ in range(int(POLDEG))])


def sample_zo_poly(rho=ZO_RHO):
    rho_f = float(rho)
    coeffs = []
    for _ in range(int(POLDEG)):
        u = _py_random.random()
        if u < rho_f/2.0:
            coeffs.append(ZZ(-1))
        elif u < rho_f:
            coeffs.append(ZZ(1))
        else:
            coeffs.append(ZZ(0))
    return R(coeffs)


def sample_uniform_poly(mod):
    mod = ZZ(mod)
    return R([ZZ.random_element(mod) for _ in range(int(POLDEG))])


def sample_uniform_poly_seeded(mod, seed):
    """Deterministic uniform P-residue for a fixed simulated switching key."""
    rng = _py_random.Random(int(seed))
    m = int(ZZ(mod))
    return R([ZZ(rng.randrange(m)) for _ in range(int(POLDEG))])

# -----------------------------------------------------------------------------
# 3. Rank-1 q0 encryption and specialized q0 RankUp (paper Algorithm 1)
# -----------------------------------------------------------------------------


def keygen_client_q0():
    s = sample_hwt_secret_poly()
    a = sample_uniform_poly(Q0)
    e_pk = sample_error_poly()
    b = poly_mod_prime((-a*s + e_pk) % firreducible, Q0)
    return [matrix(R, 1, 1, [s]), matrix(R, 1, 1, [b]), matrix(R, 1, 1, [a])]


def keygen_eval_secret(rank):
    s = matrix(R, 1, int(rank))
    for i in range(int(rank)):
        s[0, i] = sample_hwt_secret_poly()
    return s


def encrypt_client_q0(m_poly, client_key):
    b = client_key[1][0, 0]
    a = client_key[2][0, 0]
    v = sample_zo_poly()
    e0 = sample_error_poly()
    e1 = sample_error_poly()
    c0 = poly_mod_prime((v*b + m_poly + e0) % firreducible, Q0)
    c1 = poly_mod_prime((v*a + e1) % firreducible, Q0)
    return [[c0], [[c1]]]


def _center_coeff(a, q):
    a = ZZ(a) % ZZ(q)
    if a > ZZ(q)//2:
        a -= ZZ(q)
    return ZZ(a)


def moddown_nearest_one_q(q_share, p_share, P):
    """Nearest ModDown from P*q0 to q0; also return rho=-r/P."""
    P = ZZ(P)
    q_coeffs = _dense_coeffs(q_share, POLDEG)
    p_coeffs = _dense_coeffs(p_share, POLDEG)
    Pinv = inverse_mod(P, Q0)
    out = [ZZ(0)] * int(POLDEG)
    rho = [QQ(0)] * int(POLDEG)
    for i in range(int(POLDEG)):
        rem = _center_coeff(p_coeffs[i], P)
        out[i] = ZZ((Pinv * ((ZZ(q_coeffs[i]) - rem) % Q0)) % Q0)
        rho[i] = -QQ(rem) / QQ(P)
    return R(out), R(rho)


def rankup_keygen_q0(s_client, s_eval):
    rank = int(s_eval.ncols())
    M = ZZ(P_RU) * Q0
    a_raw = [sample_uniform_poly(M) for _ in range(rank)]
    e_ruk = sample_error_poly()
    inner = R(0)
    for i in range(rank):
        inner += a_raw[i] * s_eval[0, i]
    b_raw = poly_mod_prime(
        (-inner + e_ruk + ZZ(P_RU)*s_client) % firreducible, M
    )
    return {
        "b_q": poly_mod_prime(b_raw, Q0),
        "b_p": poly_mod_prime(b_raw, P_RU),
        "a_q": [poly_mod_prime(a_raw[i], Q0) for i in range(rank)],
        "a_p": [poly_mod_prime(a_raw[i], P_RU) for i in range(rank)],
        "e_ruk": e_ruk,
    }


def rank_up_q0(ct_client, ruk, s_eval):
    rank = int(s_eval.ncols())
    c0 = ct_client[0][0]
    c1_q = ct_client[1][0][0]
    c1_tilde = center_lift_poly_mod_q(c1_q, Q0)

    x0_q = poly_mod_prime((c1_tilde*ruk["b_q"]) % firreducible, Q0)
    x0_p = poly_mod_prime((c1_tilde*ruk["b_p"]) % firreducible, P_RU)
    d0, _ = moddown_nearest_one_q(x0_q, x0_p, P_RU)

    d1 = []
    for i in range(rank):
        xi_q = poly_mod_prime((c1_tilde*ruk["a_q"][i]) % firreducible, Q0)
        xi_p = poly_mod_prime((c1_tilde*ruk["a_p"][i]) % firreducible, P_RU)
        di, _ = moddown_nearest_one_q(xi_q, xi_p, P_RU)
        d1.append(di)

    out0 = poly_mod_prime(c0 + d0, Q0)
    out1 = [[poly_mod_prime(d1[i], Q0)] for i in range(rank)]
    return [[out0], out1]


def decrypt_q0_poly(ct, secret_matrix):
    temp = ct[0][0]
    rank = len(ct[1])
    for i in range(rank):
        temp = (temp + ct[1][i][0]*secret_matrix[0, i]) % firreducible
    return center_lift_poly_mod_q(poly_mod_prime(temp, Q0), Q0)


def centered_diff_q0(a, b):
    return center_lift_poly_mod_q(
        poly_mod_prime((a-b) % firreducible, Q0), Q0
    )

# -----------------------------------------------------------------------------
# 4. Ring automorphism and exact generalized RNS key-switch noise instrument
# -----------------------------------------------------------------------------


def automorphism_poly(poly, rot_slots):
    """Apply phi_j: X -> X^(5^j) in Z[X]/(X^n+1), also for QQ coeffs."""
    exponent = pow(ZZ(5), ZZ(rot_slots), ZZ(2*POLDEG))
    N = int(POLDEG)
    M = 2*N
    e = int(ZZ(exponent) % M)
    src = _dense_coeffs(poly, N)
    dst = [QQ(0)] * N
    for i in range(N):
        a = QQ(src[i])
        pos = (i*e) % M
        if pos >= N:
            dst[pos-N] -= a
        else:
            dst[pos] += a
    return R(dst) % firreducible


def rational_poly_scale(poly, denom, sign=1):
    coeffs = _dense_coeffs(poly, POLDEG)
    d = QQ(denom)
    return R([QQ(sign)*QQ(c)/d for c in coeffs])


def canonical_inf_norm(poly):
    vals = ckks_decode(poly, scale=1, slots=SLOTS, method="fft")
    if len(vals) == 0:
        return RR(0)
    return max([RR(abs(v)) for v in vals])


def eps_from_poly(poly):
    return RR(canonical_inf_norm(poly)) / RR(DELTA)


def precision_from_eps(eps):
    eps = RR(eps)
    if eps == 0:
        return RR(+Infinity)
    return -RR(log(eps)/log(RR(2)))


class KSNoiseState(object):
    """Fixed switching-key randomness relevant to the exact error distribution."""
    def __init__(self, r_src, r_tgt, d, master_seed):
        self.r_src = int(r_src)
        self.r_tgt = int(r_tgt)
        self.d = int(d)
        self.master_seed = int(master_seed)
        # One Gaussian switching-key error for every source component/digit.
        self.errors = []
        for u in range(self.r_src):
            row = []
            for j in range(self.d):
                row.append(sample_error_poly())
            self.errors.append(row)

    def a_seed(self, u, j, t):
        # Distinct deterministic streams.  This is only experiment randomness,
        # not a cryptographic PRNG construction.
        idx = ((int(u)*self.d + int(j))*self.r_tgt + int(t))
        return self.master_seed + 1000003*idx + 9176


def generalized_ks_error_at_c2s(c1_qshares, target_secret, rot_slots, state):
    """Exact added error of nearest-ModDown generalized RNS key switching.

    c1_qshares is a list over source components, with each component already in
    the exact RNS representation present immediately before CoeffToSlot, i.e.
    after ModRaise and the fully-packed SubSum normalization/rescale.

    Only P_BOOT-side uniform key residues are materialized because the exact
    decryption error depends on the q_i key shares only through the gadget
    reconstruction, while the nearest-ModDown rounding terms are determined by
    the P residue.  Independent uniform P residues are the exact marginal of a
    uniform key polynomial modulo P_BOOT*Q_C2S.
    """
    r_src = len(c1_qshares)
    r_tgt = int(target_secret.ncols())
    d = len(C2S_QS)
    if state.r_src != r_src or state.r_tgt != r_tgt or state.d != d:
        raise ValueError("KSNoiseState shape mismatch")
    for u in range(r_src):
        if len(c1_qshares[u]) != d:
            raise ValueError("unexpected number of active C2S q-shares")

    key_acc = R(0)
    acc1_p = [R(0) for _ in range(r_tgt)]

    for u in range(r_src):
        for j in range(d):
            qj = C2S_QS[j]
            # Rotation automorphism is applied component-wise before KS.
            digit = poly_mod_prime(
                automorphism_poly(c1_qshares[u][j], rot_slots), qj
            )

            key_acc = (key_acc + digit*state.errors[u][j]) % firreducible

            digit_p = poly_mod_prime(digit, P_BOOT)
            for t in range(r_tgt):
                a_p = sample_uniform_poly_seeded(
                    P_BOOT, state.a_seed(u, j, t)
                )
                acc1_p[t] = poly_mod_prime(
                    (acc1_p[t] + digit_p*a_p) % firreducible,
                    P_BOOT
                )

    # At modulus P the gadget term P*B_j*s_src vanishes.
    acc0_expr = key_acc
    for t in range(r_tgt):
        acc0_expr -= acc1_p[t]*target_secret[0, t]
    acc0_p = poly_mod_prime(acc0_expr % firreducible, P_BOOT)

    # Nearest ModDown rounding: rho=-centered(residue_P)/P.
    rem0 = center_lift_poly_mod_q(acc0_p, P_BOOT)
    rho0 = rational_poly_scale(rem0, P_BOOT, sign=-1)

    rho1 = []
    for t in range(r_tgt):
        remt = center_lift_poly_mod_q(acc1_p[t], P_BOOT)
        rho1.append(rational_poly_scale(remt, P_BOOT, sign=-1))

    key_term = rational_poly_scale(key_acc, P_BOOT, sign=1)
    round_term = R(rho0)
    for t in range(r_tgt):
        round_term = (round_term + rho1[t]*target_secret[0, t]) % firreducible

    e_ks = (key_term + round_term) % firreducible
    return e_ks, {
        "key_term": key_term,
        "round_term": round_term,
    }


def c2s_entry_from_q0_raw(raw_q0):
    """Execute the public prefix common to the fully-packed C2S entry.

    q0 -> ModRaise to ACTIVE_Q shares -> full-pack SubSum normalization.
    For full packing the SubSum rotation loop is empty, but hb_div_by_po2()
    performs the public 1/(N/2) normalization and consumes one top RNS prime in
    the current port.  The returned ciphertext therefore has C2S_ACTIVE_Q
    shares and scale 16*q0, exactly as the real C2S call sees it.
    """
    raised = heaan_mod_raise_q0(raw_q0, active=ACTIVE_Q)
    out = ModHEBootCipher(
        raised, scale=BOOT_SCALE, slots=SLOTS, is_complex=True
    )
    logNh = ZZ(POLDEG//2).nbits() - 1
    out = hb_div_by_po2(out, int(logNh))
    if out.active() != C2S_ACTIVE_Q:
        raise ArithmeticError("unexpected active level after full-pack normalization")
    return out

# -----------------------------------------------------------------------------
# 5. Fixed plaintext
# -----------------------------------------------------------------------------

base_vals = [
     CC(0.035) + CC(0.010)*I,
    -CC(0.020) + CC(0.015)*I,
     CC(0.010) - CC(0.025)*I,
    -CC(0.030) - CC(0.005)*I,
]
VALUES = [base_vals[i % len(base_vals)] for i in range(SLOTS)]
MU = ckks_encode(VALUES, DELTA, method="fft")

# -----------------------------------------------------------------------------
# 6. Diagnostics
# -----------------------------------------------------------------------------

print("\n============================================================")
print("CROSS-RANK PRIMITIVE ERROR BENCH")
print("============================================================")
print("n                              = 2^%d = %d" % (EXP_LOGN, POLDEG))
print("target rank r                  =", R_EVAL)
print("slots                          =", SLOTS)
print("C2S baby size B                =", BABY_K)
print("tested baby rotations          =", ROTATIONS)
print("Delta bits                     =", DELTA.nbits()-1)
print("q0 bits / value                =", Q0.nbits(), Q0)
print("P_RU bits / value              =", P_RU.nbits(), P_RU)
print("ModRaise active q digits         =", ACTIVE_Q)
print("ModRaise q_i bit sizes          =", [q.nbits() for q in ACTIVE_QS])
print("C2S active q digits             =", C2S_ACTIVE_Q)
print("C2S q_i bit sizes               =", [q.nbits() for q in C2S_QS])
print("log2(Q_raise)                   =", float(log(RR(QRAISE))/log(RR(2))))
print("log2(Q_C2S)                     =", float(log(RR(QC2S))/log(RR(2))))
print("C2S scale bits (approx)         =", float(log(RR(BOOT_SCALE))/log(RR(2))))
print("P_BOOT bits / value            =", P_BOOT.nbits(), P_BOOT)
print("P_BOOT > max(q_i)              =", bool(P_BOOT > max(ACTIVE_QS)))
print("log2(P_BOOT*Q_C2S)             =",
      float(log(RR(P_BOOT*QC2S))/log(RR(2))))
print("distribution                   = HWT(%d), DG(sigma=%s), ZO(%s)" %
      (SECRET_H, SIGMA, ZO_RHO))
print("statistics                     = %d key setups x %d ciphertexts" %
      (KEY_TRIALS, CT_TRIALS))
print("NOTE: this is a numerical/correctness parameter set, not a security claim.")

# -----------------------------------------------------------------------------
# 7. Main experiment
# -----------------------------------------------------------------------------

rows = []

for kt in range(KEY_TRIALS):
    client_key = keygen_client_q0()
    sC = client_key[0][0, 0]
    sE = keygen_eval_secret(R_EVAL)

    # q0-only standalone RankUp key, fixed for all ciphertext trials in setup.
    ruk = rankup_keygen_q0(sC, sE)

    # Fixed same-rank and cross-rank rotation-key randomness per baby index.
    # They are distinct keys, as in the real construction.
    rot_states = {}
    xrot_states = {}
    for rot in ROTATIONS:
        rot_states[rot] = KSNoiseState(
            r_src=R_EVAL, r_tgt=R_EVAL, d=C2S_ACTIVE_Q,
            master_seed=SEED + 100000000*kt + 100000*rot + 11
        )
        xrot_states[rot] = KSNoiseState(
            r_src=1, r_tgt=R_EVAL, d=C2S_ACTIVE_Q,
            master_seed=SEED + 100000000*kt + 100000*rot + 29
        )

    for tt in range(CT_TRIALS):
        ctC = encrypt_client_q0(MU, client_key)
        decC = decrypt_q0_poly(ctC, client_key[0])

        ctUp = rank_up_q0(ctC, ruk, sE)
        decUp = decrypt_q0_poly(ctUp, sE)
        eRU = centered_diff_q0(decUp, decC)
        E_RU = canonical_inf_norm(eRU)
        eps_RU = RR(E_RU) / RR(DELTA)

        # Build the exact public prefix seen immediately before C2S in each path.
        # Baseline has rank r before the prefix; optimized remains rank 1.
        c2s_base = c2s_entry_from_q0_raw(ctUp)
        c2s_opt = c2s_entry_from_q0_raw(ctC)
        if c2s_base.rank() != R_EVAL or c2s_opt.rank() != 1:
            raise ArithmeticError("unexpected rank at C2S entry")
        if abs(c2s_base.scale - HB_RR(BOOT_SCALE)) > HB_RR(0):
            raise ArithmeticError("unexpected baseline C2S scale")
        if abs(c2s_opt.scale - HB_RR(BOOT_SCALE)) > HB_RR(0):
            raise ArithmeticError("unexpected optimized C2S scale")

        for rot in ROTATIONS:
            # Local r->r key-switch error of an ordinary baby rotation.
            eRot, dRot = generalized_ks_error_at_c2s(
                c2s_base.raw[1], sE, rot, rot_states[rot]
            )
            E_Rot = canonical_inf_norm(eRot)
            eps_Rot = RR(E_Rot) / RR(BOOT_SCALE)

            # Local 1->r key-switch error of the fused cross-rank baby rotation.
            eX, dX = generalized_ks_error_at_c2s(
                c2s_opt.raw[1], sE, rot, xrot_states[rot]
            )
            E_X = canonical_inf_norm(eX)
            eps_X = RR(E_X) / RR(BOOT_SCALE)

            rows.append({
                "key_trial": kt,
                "ct_trial": tt,
                "logn": EXP_LOGN,
                "n": int(POLDEG),
                "r": R_EVAL,
                "rotation": int(rot),
                "baby_k": int(BABY_K),
                "modraise_active_q": int(ACTIVE_Q),
                "c2s_active_q": int(C2S_ACTIVE_Q),
                "log2_Qraise": float(log(RR(QRAISE))/log(RR(2))),
                "log2_Qc2s": float(log(RR(QC2S))/log(RR(2))),
                "q0_bits": int(Q0.nbits()),
                "P_RU_bits": int(P_RU.nbits()),
                "P_BOOT_bits": int(P_BOOT.nbits()),
                "P_BOOT_gt_max_q": int(P_BOOT > max(ACTIVE_QS)),
                "Delta_bits": int(DELTA.nbits()-1),
                "h": int(SECRET_H),
                "sigma": float(SIGMA),

                "E_RU_can": float(E_RU),
                "eps_RU": float(eps_RU),
                "prec_RU_bits": float(precision_from_eps(eps_RU)),

                "E_Rot_can": float(E_Rot),
                "eps_Rot": float(eps_Rot),
                "prec_Rot_bits": float(precision_from_eps(eps_Rot)),
                "E_Rot_keyterm_can": float(canonical_inf_norm(dRot["key_term"])),
                "E_Rot_roundterm_can": float(canonical_inf_norm(dRot["round_term"])),

                "E_XRot_can": float(E_X),
                "eps_XRot": float(eps_X),
                "prec_XRot_bits": float(precision_from_eps(eps_X)),
                "E_XRot_keyterm_can": float(canonical_inf_norm(dX["key_term"])),
                "E_XRot_roundterm_can": float(canonical_inf_norm(dX["round_term"])),

                # Local comparison only: both are measured at the same C2S scale.
                "xrot_over_rot_ratio": float(eps_X/eps_Rot) if eps_Rot != 0 else float('inf'),
                # positive => cross-rank primitive preserves MORE local precision bits
                "xrot_minus_rot_prec_bits": float(
                    precision_from_eps(eps_X) - precision_from_eps(eps_Rot)
                ),
            })

    # compact setup printout
    setup_rows = [z for z in rows if z["key_trial"] == kt]
    def _med_setup(field):
        vals = sorted([RR(z[field]) for z in setup_rows])
        return vals[len(vals)//2]
    print(
        "setup %d/%d: eps_RU=%.3e  eps_Rot=%.3e  eps_XRot=%.3e  X/Rot=%.3f" %
        (
            kt+1, KEY_TRIALS,
            float(_med_setup("eps_RU")),
            float(_med_setup("eps_Rot")),
            float(_med_setup("eps_XRot")),
            float(_med_setup("xrot_over_rot_ratio")),
        )
    )

# -----------------------------------------------------------------------------
# 8. CSV and paper-oriented medians
# -----------------------------------------------------------------------------

if len(rows) == 0:
    raise RuntimeError("no observations produced")

with open(OUTPUT_CSV, "w") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    for row in rows:
        writer.writerow(row)


def median_field(subrows, field):
    vals = sorted([RR(z[field]) for z in subrows])
    m = len(vals)
    if m % 2 == 1:
        return vals[m//2]
    return (vals[m//2-1] + vals[m//2]) / RR(2)

summary_rows = []
for rot in ROTATIONS:
    rr = [z for z in rows if z["rotation"] == rot]
    eps_ru = median_field(rr, "eps_RU")
    eps_rot = median_field(rr, "eps_Rot")
    eps_x = median_field(rr, "eps_XRot")

    srow = {
        "logn": EXP_LOGN,
        "n": int(POLDEG),
        "r": R_EVAL,
        "rotation": int(rot),
        "baby_k": int(BABY_K),
        "modraise_active_q": int(ACTIVE_Q),
        "c2s_active_q": int(C2S_ACTIVE_Q),
        "log2_Qraise": float(log(RR(QRAISE))/log(RR(2))),
        "log2_Qc2s": float(log(RR(QC2S))/log(RR(2))),
        "q0_bits": int(Q0.nbits()),
        "P_RU_bits": int(P_RU.nbits()),
        "P_BOOT_bits": int(P_BOOT.nbits()),
        "trials": len(rr),
        "eps_RU_med": float(eps_ru),
        "prec_RU_med_bits": float(precision_from_eps(eps_ru)),
        "eps_Rot_med": float(eps_rot),
        "prec_Rot_med_bits": float(precision_from_eps(eps_rot)),
        "eps_XRot_med": float(eps_x),
        "prec_XRot_med_bits": float(precision_from_eps(eps_x)),
        "xrot_over_rot_med": float(eps_x/eps_rot) if eps_rot != 0 else float('inf'),
        "xrot_minus_rot_prec_med_bits": float(
            precision_from_eps(eps_x) - precision_from_eps(eps_rot)
        ),
    }
    summary_rows.append(srow)

summary_exists = os.path.exists(SUMMARY_CSV) and os.path.getsize(SUMMARY_CSV) > 0
with open(SUMMARY_CSV, "a") as f:
    writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
    if not summary_exists:
        writer.writeheader()
    for row in summary_rows:
        writer.writerow(row)

print("\n============================================================")
print("PAPER-ORIENTED MEDIAN SUMMARY")
print("============================================================")
for srow in summary_rows:
    print(
        "j=%d: eps_RU=%.4e (%.2f client bits) | "
        "eps_Rot=%.4e (%.2f C2S bits)  eps_XRot=%.4e (%.2f C2S bits)  "
        "X/Rot=%.3f  DeltaPrec_local=%+.2f bits" %
        (
            srow["rotation"],
            srow["eps_RU_med"], srow["prec_RU_med_bits"],
            srow["eps_Rot_med"], srow["prec_Rot_med_bits"],
            srow["eps_XRot_med"], srow["prec_XRot_med_bits"],
            srow["xrot_over_rot_med"],
            srow["xrot_minus_rot_prec_med_bits"],
        )
    )

print("saved:", OUTPUT_CSV, "and", SUMMARY_CSV)
