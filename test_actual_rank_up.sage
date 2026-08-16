# =============================================================================
# rankup_error_precision_final.sage
#
# Focused experiment for the RankUp precision/error analysis in the paper.
#
# It measures, for a rank-1 q0 ciphertext ct_C and its rank-r RankUp output:
#
#   e_C   = Dec_{s_C}(ct_C)  - mu
#   e_up  = Dec_{s_E}(ct_up) - Dec_{s_C}(ct_C)
#   e_in  = Dec_{s_E}(ct_up) - mu = e_C + e_up       (mod q0)
#
# and converts canonical-norm errors into CKKS absolute precision bits
#
#   p(e) = log2(Delta / ||e||_can,infinity).
#
# It also compares against a fresh conventional rank-r ciphertext generated
# over the complete RNS evaluation modulus Q.
#
# IMPORTANT:
#   This bench uses the concrete HEAAN-style distributions:
#       secret polynomial : HWT(h), signed {-1,0,+1}, exact weight h
#       error polynomial  : coefficient-wise discrete Gaussian DG(sigma^2)
#       encryption v      : ZO(0.5)
#
#   The q0 RankUp below is deliberately specialized to the construction in the
#   paper: one key modulo P*q0, centered lift of c1, and nearest ModDown.
#   This makes the measured primitive match Algorithm 1 / Proposition 1 rather
#   than the more general full-RNS proof-of-concept key-switching path.
#
# Run from the project source directory, e.g.
#       RU_LOGN=10 RU_RANK=2 RU_P_BITS=60 sage rankup_error_precision_final.sage
#
# The script writes per-trial data to a parameter-specific CSV and appends
# one compact row to rankup_paper_summary.csv.
# =============================================================================

import os

# Load the project modulus chain, then override only the ring degree/rank for
# this isolated RankUp experiment.  This avoids editing constants.sage between
# runs.  The q-primes themselves are kept fixed across runs so that changing n
# measures the ring-degree effect rather than a simultaneous modulus change.
load("constants.sage")

EXP_LOGN = int(os.environ.get("RU_LOGN", "12"))
EXP_RANK = int(os.environ.get("RU_RANK", "2"))
POLDEG = ZZ(2)**ZZ(EXP_LOGN)
MLWE = ZZ(EXP_RANK)

# Rebuild the polynomial ring after overriding POLDEG.
R.<x> = QQ[]
firreducible = x**POLDEG + 1

# Only these helpers are needed; key generation/encryption are implemented
# locally with the experiment distributions below.
load("auxiliary_functions.sage")
load("encode.sage")

from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
import random as _py_random
import csv

# =============================================================================
# 0. Experiment configuration
# =============================================================================

R_CLIENT = 1
R_EVAL = int(MLWE)
SLOTS = int(POLDEG // 2)
DELTA = ZZ(p_enc)
Q0 = ZZ(q_list[0])

# Dedicated RankUp auxiliary primes for the precision-sensitivity experiment.
# Both are prime and satisfy P = 1 (mod 2*2^14), hence they are NTT-friendly
# for every tested n in {2^10, 2^12, 2^14}.
#
# P54 is an intentionally smaller sensitivity point.  It is below q0 and is
# NOT the nominal d=1 DPad-style choice; it is included only to demonstrate
# the dependence of the key-switching error on P.
#
# P60 is the nominal choice: it matches q0's 60-bit word size and satisfies
# P >= q0 for the single-digit RankUp specialization.
#
# These primes are LOCAL to this benchmark and do not replace the project's
# existing special modulus p used by the bootstrap/relinearization prototype.
P54 = ZZ(9007199379521537)
P60 = ZZ(576460752341598209)

P_BITS_CHOICE = int(os.environ.get("RU_P_BITS", "60"))
if "RU_P" in os.environ:
    # Optional exact override for reproducibility experiments.
    P_KS = ZZ(os.environ["RU_P"])
elif P_BITS_CHOICE == 54:
    P_KS = P54
elif P_BITS_CHOICE == 60:
    P_KS = P60
else:
    raise ValueError("RU_P_BITS must be 54 or 60 (or set exact RU_P)")

# Concrete distribution parameters.  h=64 and sigma=3.2 match the commonly
# used HEAAN bootstrapping experiment convention.  Change these only if the
# parameter table in the paper uses a different concrete instantiation.
SECRET_H = 64
SIGMA = RR(3.2)
ZO_RHO = RR(0.5)

# Statistical structure:
#   - regenerate keys / RankUp key KEY_TRIALS times;
#   - for each fixed key setup, sample CT_TRIALS independent encryptions.
# This separates variability due to setup material from per-ciphertext noise.
KEY_TRIALS = int(os.environ.get("RU_KEY_TRIALS", "5"))
CT_TRIALS = int(os.environ.get("RU_CT_TRIALS", "20"))
SEED = int(os.environ.get("RU_SEED", "20260814"))

OUTPUT_CSV = os.environ.get(
    "RU_OUTPUT_CSV",
    "rankup_trials_n2^%d_r%d_P%d.csv" %
        (EXP_LOGN, EXP_RANK, P_KS.nbits())
)
SUMMARY_CSV = os.environ.get("RU_SUMMARY_CSV", "rankup_paper_summary.csv")

if R_EVAL < 2:
    raise ValueError("Set MLWE >= 2 for a non-trivial RankUp experiment")
if SECRET_H > POLDEG:
    raise ValueError(
        "SECRET_H=%d exceeds POLDEG=%d. Increase POLDEG or choose the h used "
        "by the intended concrete parameter set." % (SECRET_H, POLDEG)
    )
if SLOTS < 1:
    raise ValueError("POLDEG must be at least 2")
if not Q0.is_prime():
    raise ValueError("q0 must be prime")
if not P_KS.is_prime():
    raise ValueError("P_KS must be prime")
# P54 deliberately violates P >= q0; this is allowed only as an experimental
# sensitivity point.  The arithmetic remains well-defined, but the 54-bit row
# should not be presented as the nominal DPad-aligned parameter choice.
P_UNDERSIZED = bool(P_KS < Q0)
if (Q0 - 1) % (2*POLDEG) != 0:
    raise ValueError("q0 is not 1 mod 2n for this POLDEG")
if (P_KS - 1) % (2*POLDEG) != 0:
    raise ValueError("P_KS is not 1 mod 2n for this POLDEG")

set_random_seed(int(SEED))
_py_random.seed(int(SEED))
_DG = DiscreteGaussianDistributionIntegerSampler(sigma=SIGMA)


# =============================================================================
# 1. Concrete samplers used by this experiment
# =============================================================================

def sample_hwt_secret_poly(h=SECRET_H):
    """Uniform signed binary polynomial of exact Hamming weight h."""
    h = int(h)
    positions = _py_random.sample(range(int(POLDEG)), h)
    coeffs = [ZZ(0)] * int(POLDEG)
    for i in positions:
        coeffs[i] = ZZ(1) if _py_random.getrandbits(1) else ZZ(-1)
    return R(coeffs)


def sample_error_poly():
    """Coefficient-wise discrete Gaussian DG(SIGMA^2)."""
    return R([ZZ(_DG()) for _ in range(int(POLDEG))])


def sample_zo_poly(rho=ZO_RHO):
    """ZO(rho): +/-1 with probability rho/2 each, 0 otherwise."""
    rho_f = float(rho)
    coeffs = []
    for _ in range(int(POLDEG)):
        u = _py_random.random()
        if u < rho_f / 2.0:
            coeffs.append(ZZ(-1))
        elif u < rho_f:
            coeffs.append(ZZ(1))
        else:
            coeffs.append(ZZ(0))
    return R(coeffs)


def sample_uniform_poly(mod):
    mod = ZZ(mod)
    return R([ZZ.random_element(mod) for _ in range(int(POLDEG))])


# =============================================================================
# 2. Paper-aligned key generation / encryption helpers
# =============================================================================

def keygen_client_q0():
    """Rank-1 public/secret key generated directly modulo q0."""
    s = sample_hwt_secret_poly()
    a = sample_uniform_poly(Q0)
    e_pk = sample_error_poly()
    b = poly_mod_prime((-a*s + e_pk) % firreducible, Q0)

    # Full-key-like layout so the secret is easy to access consistently.
    smat = matrix(R, 1, 1, [s])
    bmat = matrix(R, 1, 1, [b])
    Amat = matrix(R, 1, 1, [a])
    return [smat, bmat, Amat], e_pk


def keygen_eval_full(rank):
    """Conventional rank-r MLWE key over the complete RNS modulus Q."""
    rank = int(rank)
    s = matrix(R, 1, rank)
    for j in range(rank):
        s[0, j] = sample_hwt_secret_poly()

    A = matrix(R, rank, rank)
    for i in range(rank):
        for j in range(rank):
            A[i, j] = sample_uniform_poly(Q)

    e_pk = matrix(R, rank, 1)
    b = matrix(R, rank, 1)
    for i in range(rank):
        e_pk[i, 0] = sample_error_poly()
        inner = R(0)
        for j in range(rank):
            inner += A[i, j] * s[0, j]
        b[i, 0] = poly_mod_prime((-inner + e_pk[i, 0]) % firreducible, Q)

    return [s, b, A], e_pk


def encrypt_client_q0(m_poly, client_key):
    """Paper's rank-1, q0-only client encryption."""
    b = client_key[1][0, 0]
    a = client_key[2][0, 0]

    v = sample_zo_poly()
    e0 = sample_error_poly()
    e1 = sample_error_poly()

    c0 = poly_mod_prime((v*b + m_poly + e0) % firreducible, Q0)
    c1 = poly_mod_prime((v*a + e1) % firreducible, Q0)
    return [[c0], [[c1]]]


def encrypt_fresh_rank_r_full(m_poly, eval_key):
    """Conventional fresh rank-r encryption carrying all L RNS q-shares."""
    rank = int(eval_key[0].ncols())
    b = eval_key[1]
    A = eval_key[2]

    v = [sample_zo_poly() for _ in range(rank)]
    e0 = sample_error_poly()
    e1 = [sample_error_poly() for _ in range(rank)]

    c0 = R(0)
    for i in range(rank):
        c0 += v[i] * b[i, 0]
    c0 = poly_mod_prime((c0 + m_poly + e0) % firreducible, Q)

    c1 = []
    for j in range(rank):
        acc = R(0)
        for i in range(rank):
            acc += v[i] * A[i, j]
        c1.append(poly_mod_prime((acc + e1[j]) % firreducible, Q))

    c0_qshares = shares_in_mod(c0, L, q_list)
    c1_qshares = [shares_in_mod(c1[j], L, q_list) for j in range(rank)]
    return [c0_qshares, c1_qshares]


# =============================================================================
# 3. Specialized q0 RankUp matching Algorithm 1 in the paper
# =============================================================================

def rankup_keygen_q0(s_client, s_eval):
    """Generate ruk_{C->E} = (b_up, a_up) only modulo P*q0.

    b_up = -<a_up,s_E> + e_ruk + P*s_C  (mod P*q0).

    We keep q0/P shares plus the sampled e_ruk for diagnostics.
    """
    rank = int(s_eval.ncols())
    M = P_KS * Q0

    a_raw = [sample_uniform_poly(M) for _ in range(rank)]
    e_ruk = sample_error_poly()

    inner = R(0)
    for i in range(rank):
        inner += a_raw[i] * s_eval[0, i]

    b_raw = poly_mod_prime((-inner + e_ruk + P_KS*s_client) % firreducible, M)

    return {
        "b_q": poly_mod_prime(b_raw, Q0),
        "b_p": poly_mod_prime(b_raw, P_KS),
        "a_q": [poly_mod_prime(a_raw[i], Q0) for i in range(rank)],
        "a_p": [poly_mod_prime(a_raw[i], P_KS) for i in range(rank)],
        "e_ruk": e_ruk,
    }


def _center_coeff(a, q):
    a = ZZ(a) % ZZ(q)
    if a > ZZ(q)//2:
        a -= ZZ(q)
    return ZZ(a)


def moddown_nearest_q0(q_share, p_share):
    """Nearest ModDown from P*q0 to q0, coefficient-wise.

    If x = P*y + r with centered remainder r in (-P/2,P/2], return y mod q0.
    Also return rho = -r/P so that y = x/P + rho.
    """
    q_coeffs = _dense_coeffs(q_share, POLDEG)
    p_coeffs = _dense_coeffs(p_share, POLDEG)
    Pinv = inverse_mod(P_KS, Q0)

    out = [ZZ(0)] * int(POLDEG)
    rho = [QQ(0)] * int(POLDEG)
    for i in range(int(POLDEG)):
        r = _center_coeff(p_coeffs[i], P_KS)
        out[i] = ZZ((Pinv * ((ZZ(q_coeffs[i]) - r) % Q0)) % Q0)
        rho[i] = -QQ(r) / QQ(P_KS)

    return R(out), R(rho)


def rank_up_q0_instrumented(ct_client, ruk, s_eval):
    """RankUp plus exact decomposition diagnostics for e_up."""
    rank = int(s_eval.ncols())
    c0 = ct_client[0][0]
    c1_q = ct_client[1][0][0]

    # Algorithm 1: public coefficient-wise centered lift of c1.
    c1_tilde = center_lift_poly_mod_q(c1_q, Q0)

    x0_q = poly_mod_prime((c1_tilde * ruk["b_q"]) % firreducible, Q0)
    x0_p = poly_mod_prime((c1_tilde * ruk["b_p"]) % firreducible, P_KS)
    d0, rho0 = moddown_nearest_q0(x0_q, x0_p)

    d1 = []
    rho1 = []
    for i in range(rank):
        xi_q = poly_mod_prime((c1_tilde * ruk["a_q"][i]) % firreducible, Q0)
        xi_p = poly_mod_prime((c1_tilde * ruk["a_p"][i]) % firreducible, P_KS)
        di, rhoi = moddown_nearest_q0(xi_q, xi_p)
        d1.append(di)
        rho1.append(rhoi)

    out0 = poly_mod_prime(c0 + d0, Q0)
    out1 = [[poly_mod_prime(d1[i], Q0)] for i in range(rank)]

    # Proposition 1 decomposition:
    # e_up = P^{-1} c1_tilde e_ruk + rho0 + <rho1,s_E>  (mod q0).
    key_term = R((c1_tilde * ruk["e_ruk"]) % firreducible) / QQ(P_KS)
    round_term = R(rho0)
    for i in range(rank):
        round_term += (rho1[i] * s_eval[0, i]) % firreducible
    formula = (key_term + round_term) % firreducible

    return [ [out0], out1 ], {
        "c1_tilde": c1_tilde,
        "rho0": rho0,
        "rho1": rho1,
        "key_term": key_term,
        "round_term": round_term,
        "formula": formula,
    }


# =============================================================================
# 4. Decryption / norm / precision helpers
# =============================================================================

def decrypt_q0_poly(ct, secret_matrix):
    """Centered q0 decryption polynomial, with no decoding."""
    temp = ct[0][0]
    rank = len(ct[1])
    for i in range(rank):
        temp = (temp + ct[1][i][0] * secret_matrix[0, i]) % firreducible
    return center_lift_poly_mod_q(poly_mod_prime(temp, Q0), Q0)


def centered_diff_q0(a, b):
    return center_lift_poly_mod_q(
        poly_mod_prime((a - b) % firreducible, Q0), Q0
    )


def rational_poly_center_mod_q(poly, q=Q0):
    """Reduce a QQ-coefficient polynomial modulo q and center-lift it."""
    coeffs = _dense_coeffs(poly, POLDEG)
    out = []
    for c in coeffs:
        c = QQ(c)
        num = ZZ(c.numerator())
        den = ZZ(c.denominator())
        z = (num * inverse_mod(den % q, q)) % q
        out.append(_center_coeff(z, q))
    return R(out)


def canonical_inf_norm(poly):
    """Canonical embedding infinity norm for a real-coefficient ring element."""
    vals = ckks_decode(poly, scale=1, slots=SLOTS, method="fft")
    if len(vals) == 0:
        return RR(0)
    return max([RR(abs(v)) for v in vals])


def slot_vector(poly, scale=DELTA):
    return ckks_decode(poly, scale=scale, slots=SLOTS, method="fft")


def max_slot_error(got, expected):
    return max([RR(abs(got[i] - expected[i])) for i in range(len(expected))])


def max_slot_difference(a, b):
    return max([RR(abs(a[i] - b[i])) for i in range(len(a))])


def precision_bits_from_noise(E, scale=DELTA):
    E = RR(E)
    if E == 0:
        return RR(+Infinity)
    return RR(log(RR(scale) / E) / log(RR(2)))


def precision_bits_from_slot_error(err):
    err = RR(err)
    if err == 0:
        return RR(+Infinity)
    return RR(-log(err) / log(RR(2)))


def quantile(values, q):
    vals = sorted([RR(v) for v in values])
    if len(vals) == 1:
        return vals[0]
    pos = RR(q) * (len(vals) - 1)
    lo = int(floor(pos))
    hi = int(ceil(pos))
    if lo == hi:
        return vals[lo]
    alpha = pos - lo
    return (RR(1)-alpha)*vals[lo] + alpha*vals[hi]


# =============================================================================
# 5. Fixed plaintext and encoding-only error
# =============================================================================

base_vals = [
     CC(0.035) + CC(0.010)*I,
    -CC(0.020) + CC(0.015)*I,
     CC(0.010) - CC(0.025)*I,
    -CC(0.030) - CC(0.005)*I,
]
VALUES = [base_vals[i % len(base_vals)] for i in range(SLOTS)]
MU = ckks_encode(VALUES, DELTA, method="fft")
ENCODED_SLOTS = slot_vector(MU, DELTA)
ENCODING_SLOT_ERROR = max_slot_error(ENCODED_SLOTS, VALUES)

print("\nRankUp precision experiment")
print("n=2^%d, r=%d, q0=%d bits, P=%d bits, Delta=2^%d, %d setups x %d ciphertexts" %
      (EXP_LOGN, R_EVAL, Q0.nbits(), P_KS.nbits(), DELTA.nbits()-1,
       KEY_TRIALS, CT_TRIALS))


# =============================================================================
# 6. Main experiment
# =============================================================================

rows = []

for kt in range(int(KEY_TRIALS)):
    client_key, client_pk_error = keygen_client_q0()
    eval_key, eval_pk_error = keygen_eval_full(R_EVAL)

    sC = client_key[0][0, 0]
    sE = eval_key[0] # just means the client rank-r secret; not to be confused with the server public evaluation keys
    ruk = rankup_keygen_q0(sC, sE)

    B_ruk = canonical_inf_norm(ruk["e_ruk"])
    B_s = max([canonical_inf_norm(sE[0, i]) for i in range(R_EVAL)])

    setup_start = len(rows)

    for tt in range(int(CT_TRIALS)):
        # ---------------------------------------------------------------------
        # A. Rank-1 q0 client encryption
        # ---------------------------------------------------------------------
        ctC = encrypt_client_q0(MU, client_key)
        decC = decrypt_q0_poly(ctC, client_key[0])
        eC = centered_diff_q0(decC, MU)
        E_C = canonical_inf_norm(eC)

        # ---------------------------------------------------------------------
        # B. Specialized q0 RankUp
        # ---------------------------------------------------------------------
        ctUp, diag = rank_up_q0_instrumented(ctC, ruk, sE)
        decUp = decrypt_q0_poly(ctUp, sE)

        # This is the key quantity: RankUp-added error only.
        eUp = centered_diff_q0(decUp, decC)
        E_up = canonical_inf_norm(eUp)

        # Total error entering bootstrap after RankUp.
        eIn = centered_diff_q0(decUp, MU)
        E_in = canonical_inf_norm(eIn)

        # Check e_in = e_C + e_up modulo q0.
        triangle_relation = centered_diff_q0(
            eIn,
            poly_mod_prime((eC + eUp) % firreducible, Q0)
        )
        triangle_relation_coeff = max([abs(ZZ(c)) for c in _dense_coeffs(triangle_relation, POLDEG)])
        if triangle_relation_coeff != 0:
            raise ArithmeticError("e_in != e_C + e_up modulo q0")

        # Check the exact Proposition-1 decomposition modulo q0.
        formula_mod = rational_poly_center_mod_q(diag["formula"], Q0)
        formula_diff = centered_diff_q0(eUp, formula_mod)
        formula_diff_coeff = max([abs(ZZ(c)) for c in _dense_coeffs(formula_diff, POLDEG)])
        if formula_diff_coeff != 0:
            raise ArithmeticError(
                "instrumented e_up decomposition mismatch (max coeff %s)" % formula_diff_coeff
            )

        # Components entering the generic bound (Eq. 6 in the draft), using
        # observed norms instead of inserting a loose closed-form estimate.
        B_c1 = canonical_inf_norm(diag["c1_tilde"])
        B_md = max(
            [canonical_inf_norm(diag["rho0"])] +
            [canonical_inf_norm(r) for r in diag["rho1"]]
        )
        generic_bound = (B_c1 * B_ruk / RR(P_KS)) + B_md * (RR(1) + RR(R_EVAL)*B_s)

        key_term_norm = canonical_inf_norm(diag["key_term"])
        round_term_norm = canonical_inf_norm(diag["round_term"])

        if E_up > generic_bound * (RR(1) + RR(1e-20)):
            raise ArithmeticError("observed RankUp error exceeds generic norm bound")

        # ---------------------------------------------------------------------
        # C. Fresh conventional rank-r ciphertext over the complete RNS Q
        # ---------------------------------------------------------------------
        ctFresh = encrypt_fresh_rank_r_full(MU, eval_key)
        decFresh = decrypt_q0_poly(ctFresh, sE)
        eFresh = centered_diff_q0(decFresh, MU)
        E_fresh = canonical_inf_norm(eFresh)

        # ---------------------------------------------------------------------
        # D. Slot-domain accuracy (includes CKKS encoding rounding)
        # ---------------------------------------------------------------------
        slotsC = slot_vector(decC, DELTA)
        slotsUp = slot_vector(decUp, DELTA)
        slotsFresh = slot_vector(decFresh, DELTA)

        slot_err_C = max_slot_error(slotsC, VALUES)
        slot_err_Up = max_slot_error(slotsUp, VALUES)
        slot_err_Fresh = max_slot_error(slotsFresh, VALUES)
        slot_perturb_Up = max_slot_difference(slotsUp, slotsC)

        # Direct correspondence between canonical RankUp error and slot perturbation.
        expected_slot_perturb = E_up / RR(DELTA)

        rows.append({
            "key_trial": kt,
            "ct_trial": tt,
            "n": int(POLDEG),
            "r": int(R_EVAL),
            "h": int(SECRET_H),
            "sigma": float(SIGMA),
            "q0_bits": int(Q0.nbits()),
            "q0_value": str(Q0),
            "P_bits": int(P_KS.nbits()),
            "P_value": str(P_KS),
            "P_below_q0": int(P_UNDERSIZED),
            "Delta_bits": int(DELTA.nbits()-1),
            "encoding_slot_error": float(ENCODING_SLOT_ERROR),

            "E_client_can": float(E_C),
            "eps_client": float(RR(E_C) / RR(DELTA)),
            "prec_client_bits": float(precision_bits_from_noise(E_C)),

            "E_rankup_added_can": float(E_up),
            "eps_rankup_added": float(RR(E_up) / RR(DELTA)),
            "prec_rankup_added_bits": float(precision_bits_from_noise(E_up)),

            "E_after_rankup_can": float(E_in),
            "eps_after_rankup": float(RR(E_in) / RR(DELTA)),
            "prec_after_rankup_bits": float(precision_bits_from_noise(E_in)),
            "precision_loss_bits": float(precision_bits_from_noise(E_C) - precision_bits_from_noise(E_in)),
            "prec_after_triangle_bits": float(precision_bits_from_noise(E_C + E_up)),

            "E_fresh_rankr_can": float(E_fresh),
            "eps_fresh_rankr": float(RR(E_fresh) / RR(DELTA)),
            "prec_fresh_rankr_bits": float(precision_bits_from_noise(E_fresh)),

            "slot_err_client": float(slot_err_C),
            "slot_prec_client_bits": float(precision_bits_from_slot_error(slot_err_C)),
            "slot_err_after_rankup": float(slot_err_Up),
            "slot_prec_after_rankup_bits": float(precision_bits_from_slot_error(slot_err_Up)),
            "slot_err_fresh_rankr": float(slot_err_Fresh),
            "slot_prec_fresh_rankr_bits": float(precision_bits_from_slot_error(slot_err_Fresh)),
            "slot_rankup_perturb": float(slot_perturb_Up),
            "slot_rankup_perturb_bits": float(precision_bits_from_slot_error(slot_perturb_Up)),
            "slot_rankup_perturb_from_can": float(expected_slot_perturb),

            "B_c1_observed": float(B_c1),
            "B_ruk_observed": float(B_ruk),
            "B_md_observed": float(B_md),
            "B_s_observed": float(B_s),
            "rankup_generic_bound": float(generic_bound),
            "rankup_bound_precision_bits": float(precision_bits_from_noise(generic_bound)),
            "rankup_key_term_can": float(key_term_norm),
            "rankup_round_term_can": float(round_term_norm),
        })


    # One concise line per key setup: medians over the CT_TRIALS ciphertexts
    # generated under this fixed set of client/evaluation/RankUp keys.
    setup_rows = rows[setup_start:]

    eps_ru_setup = RR(quantile(
        [RR(x["slot_rankup_perturb"]) for x in setup_rows], RR(0.5)
    ))
    p_C_setup = RR(quantile(
        [RR(x["slot_prec_client_bits"]) for x in setup_rows], RR(0.5)
    ))
    p_RU_setup = precision_bits_from_slot_error(eps_ru_setup)
    p_after_setup = RR(quantile(
        [RR(x["slot_prec_after_rankup_bits"]) for x in setup_rows], RR(0.5)
    ))
    p_fresh_setup = RR(quantile(
        [RR(x["slot_prec_fresh_rankr_bits"]) for x in setup_rows], RR(0.5)
    ))

    print(
        "setup %d/%d: eps_RU=%.4e  p_C=%.2f  p_RU=%.2f  "
        "p_after=%.2f  p_fresh=%.2f"
        %
        (
            kt + 1, KEY_TRIALS,
            float(eps_ru_setup),
            float(p_C_setup),
            float(p_RU_setup),
            float(p_after_setup),
            float(p_fresh_setup)
        )
    )


# =============================================================================
# 7. CSV output
# =============================================================================

fieldnames = list(rows[0].keys())
with open(OUTPUT_CSV, "w") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        writer.writerow(row)


# =============================================================================
# 8. Compact paper-oriented summary
# =============================================================================

# For the paper we deliberately report ONE central value per parameter set:
# the median over all observations.  This follows the common CKKS experimental
# convention of repeating randomized trials and reporting median precision,
# without adding percentile columns to the main table.

"""
def med(field):
    return RR(quantile([RR(r[field]) for r in rows], RR(0.5)))


EPS_RU_MED = med("slot_rankup_perturb")
P_RU_MED = precision_bits_from_slot_error(EPS_RU_MED)

paper_row = {
    "logn": int(EXP_LOGN),
    "n": int(POLDEG),
    "r": int(R_EVAL),
    "log2_q0": int(Q0.nbits()),
    "q0_value": str(Q0),
    "log2_P": int(P_KS.nbits()),
    "P_value": str(P_KS),
    "P_below_q0": int(P_UNDERSIZED),
    "log2_Delta": int(DELTA.nbits()-1),
    "h": int(SECRET_H),
    "sigma": float(SIGMA),
    "trials": int(len(rows)),
    "eps_rankup": float(EPS_RU_MED),
    "prec_rankup_bits": float(P_RU_MED),
    "prec_input_bits": float(med("slot_prec_client_bits")),
    "prec_output_bits": float(med("slot_prec_after_rankup_bits")),
    "prec_fresh_rankr_bits": float(med("slot_prec_fresh_rankr_bits")),
}
"""

def maxf(field):
    return max([RR(r[field]) for r in rows])

def minf(field):
    return min([RR(r[field]) for r in rows])

EPS_RU_MAX = maxf("slot_rankup_perturb")

paper_row = {
    "logn": int(EXP_LOGN),
    "n": int(POLDEG),
    "r": int(R_EVAL),
    "log2_q0": int(Q0.nbits()),
    "log2_P": int(P_KS.nbits()),
    "log2_Delta": int(DELTA.nbits()-1),
    "trials": int(len(rows)),

    "eps_rankup_max": float(EPS_RU_MAX),
    "prec_rankup_min_bits":
        float(precision_bits_from_slot_error(EPS_RU_MAX)),

    "prec_input_min_bits": float(minf("slot_prec_client_bits")),
    "prec_output_min_bits": float(minf("slot_prec_after_rankup_bits")),
    "prec_fresh_min_bits": float(minf("slot_prec_fresh_rankr_bits")),
}

"""
print("median/%d observations: eps_RU=%.4e  p_C=%.2f  p_RU=%.2f  "
      "p_after=%.2f  p_fresh=%.2f" %
      (
          len(rows),
          float(paper_row["eps_rankup"]),
          float(paper_row["prec_input_bits"]),
          float(paper_row["prec_rankup_bits"]),
          float(paper_row["prec_output_bits"]),
          float(paper_row["prec_fresh_rankr_bits"])
      ))

# Append one compact row per (n,r,P) run.  The exact q0 and P values are kept
# in the CSV for reproducibility even if the paper table only shows bit sizes.
summary_fields = list(paper_row.keys())
summary_exists = os.path.exists(SUMMARY_CSV) and os.path.getsize(SUMMARY_CSV) > 0
with open(SUMMARY_CSV, "a") as f:
    writer = csv.DictWriter(f, fieldnames=summary_fields)
    if not summary_exists:
        writer.writeheader()
    writer.writerow(paper_row)

print("saved:", OUTPUT_CSV, "and", SUMMARY_CSV)
"""
print(
    "maximum/%d observations: eps_RU=%.4e  "
    "p_RU(min)=%.2f  p_C(min)=%.2f  "
    "p_after(min)=%.2f  p_fresh(min)=%.2f"
    %
    (
        len(rows),
        float(paper_row["eps_rankup_max"]),
        float(paper_row["prec_rankup_min_bits"]),
        float(paper_row["prec_input_min_bits"]),
        float(paper_row["prec_output_min_bits"]),
        float(paper_row["prec_fresh_min_bits"]),
    )
)


# Append one compact row per (n,r,P) run.
summary_fields = list(paper_row.keys())

summary_exists = (
    os.path.exists(SUMMARY_CSV)
    and os.path.getsize(SUMMARY_CSV) > 0
)

with open(SUMMARY_CSV, "a") as f:
    writer = csv.DictWriter(f, fieldnames=summary_fields)
    if not summary_exists:
        writer.writeheader()
    writer.writerow(paper_row)

print("saved:", OUTPUT_CSV, "and", SUMMARY_CSV)