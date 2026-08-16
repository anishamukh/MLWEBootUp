# -----------------------------------------------------------------------------
# BootstrapOps.sage
# Traditional CKKS-style bootstrapping scaffold adapted to the current ModHE code.
#
# This file integrates:
#   1. automorphism on polynomials/ciphertexts
#   2. generic ModHE key switching
#   3. rotation/conjugation
#   4. plaintext/scalar multiplication
#   5. linear_transform, coeff_to_slot, slot_to_coeff
#   6. mod_raise
#   7. eval_exp/eval_mod placeholder via scaled sine polynomial
#   8. img_extract
#   9. bootstrap wrapper
#
# It is a correctness/integration scaffold.  The real CKKS bootstrapping
# constants (CoeffToSlot/SlotToCoeff diagonals and EvalMod approximation) must be
# supplied separately for production-precision experiments.
# -----------------------------------------------------------------------------


# ------------------------- representation helpers ----------------------------

def ct0_as_list(ct):
    try:
        return [ct[0][0, j] for j in range(ct[0].ncols())]
    except Exception:
        return list(ct[0])


def ct1_as_lists(ct):
    try:
        return [[ct[1][i, j] for j in range(ct[1].ncols())] for i in range(ct[1].nrows())]
    except Exception:
        return [list(row) for row in ct[1]]


def normalize_ct(ct):
    return [ct0_as_list(ct), ct1_as_lists(ct)]


def ct_rank(ct):
    return len(ct1_as_lists(ct))


def ct_active(ct):
    return len(ct0_as_list(ct))


def zero_ct(active=None, rank=None):
    if active is None:
        active = L
    if rank is None:
        rank = MLWE
    return [[R(0) for _ in range(active)], [[R(0) for _ in range(active)] for _ in range(rank)]]


def copy_ct(ct):
    ct = normalize_ct(ct)
    return [[ct[0][j] for j in range(len(ct[0]))],
            [[ct[1][i][j] for j in range(len(ct[0]))] for i in range(len(ct[1]))]]


def _coeff(poly, i):
    coeffs = _dense_coeffs(poly, POLDEG)
    return coeffs[i]


def _output_mod(out_idx, active):
    if out_idx < active:
        return q_list[out_idx]
    return p


def _key_share_index(out_idx, active):
    # key shares are stored over q_0,...,q_{L-1},p, so p is always index L
    if out_idx < active:
        return out_idx
    return L


# ----------------------------- plaintext helpers -----------------------------

def scalar_to_plain_poly(alpha, encode=True):
    """Convert a Python/Sage scalar to a plaintext polynomial.

    encode=True matches the toy current code, which encodes real constants as
    round(alpha*p_enc) in the constant coefficient.  For already-scaled EvalMod
    coefficients, use encode=False or pass a polynomial directly.
    """
    if encode:
        return R(Integer(round(float(alpha)*p_enc)))
    return R(alpha)


def as_plain_poly(plain, encode=True):
    try:
        plain.degree()
        return R(plain) % firreducible
    except Exception:
        return scalar_to_plain_poly(plain, encode=encode)


def _scalar_mod_q(alpha, q):
    """Reduce an integer/rational scalar modulo q."""
    a = QQ(alpha)
    num = Integer(a.numerator())
    den = Integer(a.denominator())
    return (num % q) * (den.inverse_mod(Integer(q))) % q


def poly_mod_q_allow_rationals(poly, q):
    """Reduce a QQ polynomial coefficient-wise modulo q, allowing rationals."""
    coeffs = _dense_coeffs(poly % firreducible, POLDEG)
    return R([_scalar_mod_q(coeffs[i], q) for i in range(POLDEG)])


def poly_to_rns_shares(poly, basis):
    return [poly_mod_q_allow_rationals(poly, Integer(basis[j])) for j in range(len(basis))]


# ------------------------------ automorphisms --------------------------------

def automorphism_poly(poly, k, modulus=None):
    """Apply kappa_k: a(x) -> a(x^k) mod x^N+1.

    This implementation maps coefficients directly.  For x^e in Z[x]/(x^N+1),
    reduce e modulo 2N; if e >= N, use -x^{e-N}.
    """
    k = Integer(k) % (2*POLDEG)
    if gcd(k, 2*POLDEG) != 1:
        raise ValueError("automorphism exponent k must be odd/coprime to 2*POLDEG")

    coeffs = _dense_coeffs(poly, POLDEG)
    out = [0 for _ in range(POLDEG)]
    for i in range(POLDEG):
        c = coeffs[i]
        if c == 0:
            continue
        e = (Integer(i)*k) % (2*POLDEG)
        if e >= POLDEG:
            out[int(e-POLDEG)] -= c
        else:
            out[int(e)] += c
    res = R(out) % firreducible
    if modulus is not None:
        res = poly_mod_prime(res, modulus)
    return res


def automorphism_secret(s_or_key, k):
    ss = _bs_secret_to_list(s_or_key)
    return _bs_secret_matrix([automorphism_poly(si, k) for si in ss])


def automorphism_ct(ct, k):
    ct = normalize_ct(ct)
    active = len(ct[0])
    rank = len(ct[1])

    c0 = [automorphism_poly(ct[0][j], k, q_list[j]) for j in range(active)]
    c1 = []
    for i in range(rank):
        c1.append([automorphism_poly(ct[1][i][j], k, q_list[j]) for j in range(active)])
    return [c0, c1]


# ------------------------- basic homomorphic ops -----------------------------

def add_ct(ct1, ct2):
    ct1 = normalize_ct(ct1)
    ct2 = normalize_ct(ct2)
    active = min(len(ct1[0]), len(ct2[0]))
    rank = min(len(ct1[1]), len(ct2[1]))

    c0 = [poly_mod_prime((ct1[0][j] + ct2[0][j]) % firreducible, q_list[j]) for j in range(active)]
    c1 = []
    for i in range(rank):
        c1.append([poly_mod_prime((ct1[1][i][j] + ct2[1][i][j]) % firreducible, q_list[j]) for j in range(active)])
    return [c0, c1]


def sub_ct(ct1, ct2):
    ct1 = normalize_ct(ct1)
    ct2 = normalize_ct(ct2)
    active = min(len(ct1[0]), len(ct2[0]))
    rank = min(len(ct1[1]), len(ct2[1]))

    c0 = [poly_mod_prime((ct1[0][j] - ct2[0][j]) % firreducible, q_list[j]) for j in range(active)]
    c1 = []
    for i in range(rank):
        c1.append([poly_mod_prime((ct1[1][i][j] - ct2[1][i][j]) % firreducible, q_list[j]) for j in range(active)])
    return [c0, c1]


def neg_ct(ct):
    return scalar_multiply(ct, -1, encode=False)


def scalar_multiply(ct, alpha, encode=False):
    """Multiply all ciphertext components by an integer/rational scalar.

    Rational scalars are interpreted modulo each q_j using denominator inverses.
    """
    ct = normalize_ct(ct)
    active = len(ct[0])
    rank = len(ct[1])

    c0 = []
    for j in range(active):
        aa = _scalar_mod_q(alpha, q_list[j])
        c0.append(poly_mod_prime((aa*ct[0][j]) % firreducible, q_list[j]))

    c1 = []
    for i in range(rank):
        row = []
        for j in range(active):
            aa = _scalar_mod_q(alpha, q_list[j])
            row.append(poly_mod_prime((aa*ct[1][i][j]) % firreducible, q_list[j]))
        c1.append(row)
    return [c0, c1]


def plaintext_multiply(ct, plain, encode=True):
    """Multiply ciphertext by a known plaintext polynomial/scalar.

    For linear transforms, pass plaintext polynomials directly.  For a toy scalar
    message constant, pass encode=True.  For already-scaled polynomial constants,
    pass encode=False.
    """
    ct = normalize_ct(ct)
    active = len(ct[0])
    rank = len(ct[1])
    pt = as_plain_poly(plain, encode=encode)
    pt_shares = poly_to_rns_shares(pt, [q_list[j] for j in range(active)])

    c0 = [poly_mod_prime((ct[0][j]*pt_shares[j]) % firreducible, q_list[j]) for j in range(active)]
    c1 = []
    for i in range(rank):
        c1.append([poly_mod_prime((ct[1][i][j]*pt_shares[j]) % firreducible, q_list[j]) for j in range(active)])
    return [c0, c1]


def add_plain(ct, plain, encode=True):
    ct = copy_ct(ct)
    active = len(ct[0])
    pt = as_plain_poly(plain, encode=encode)
    pt_shares = poly_to_rns_shares(pt, [q_list[j] for j in range(active)])
    for j in range(active):
        ct[0][j] = poly_mod_prime((ct[0][j] + pt_shares[j]) % firreducible, q_list[j])
    return ct


def rescale_ct(ct, level_in=None, level_out=None):
    ct = normalize_ct(ct)
    if level_in is None:
        level_in = len(ct[0]) - 1
    if level_out is None:
        level_out = level_in - 1
    if level_out < 0:
        raise ValueError("cannot rescale below q0")
    c0 = rescaling_ciphertext(ct[0], level_in, level_out)
    c1 = [rescaling_ciphertext(ct[1][i], level_in, level_out) for i in range(len(ct[1]))]
    return [c0, c1]


# ----------------------------- key switching ---------------------------------

def keyswitch(ct, ksk):
    """Generic ModHE key-switch.

    Input ct is decryptable under ksk['source']; output is decryptable under the
    target secret used in ks_keygen().  This is the shared primitive for
    rotations, conjugation, and later bootstrapping transforms.
    """
    if isinstance(ksk, dict) and "ksk" in ksk:
        ksk = ksk["ksk"]

    ct = normalize_ct(ct)
    c0 = ct[0]
    c1 = ct[1]
    active = len(c0)
    r_from = ksk["rank_from"]
    r_to = ksk["rank_to"]

    if len(c1) != r_from:
        raise ValueError("keyswitch: ciphertext rank does not match key source rank")

    acc0 = [R(0) for _ in range(active+1)]
    acc1 = [[R(0) for _ in range(active+1)] for _ in range(r_to)]

    for out_idx in range(active+1):
        mod = _output_mod(out_idx, active)
        key_idx = _key_share_index(out_idx, active)
        a0 = R(0)
        a1 = [R(0) for _ in range(r_to)]

        for u in range(r_from):
            for j in range(active):
                cuj = c1[u][j]
                a0 = a0 + poly_mod_prime((cuj*ksk["b"][u][j][key_idx]) % firreducible, mod)
                for t in range(r_to):
                    a1[t] = a1[t] + poly_mod_prime((cuj*ksk["a"][u][t][j][key_idx]) % firreducible, mod)

        acc0[out_idx] = poly_mod_prime(a0, mod)
        for t in range(r_to):
            acc1[t][out_idx] = poly_mod_prime(a1[t], mod)

    sw0 = ModDown_poly(acc0, active)
    sw1 = [ModDown_poly(acc1[t], active) for t in range(r_to)]

    out0 = [poly_mod_prime(c0[j] + sw0[j], q_list[j]) for j in range(active)]
    out1 = []
    for t in range(r_to):
        out1.append([poly_mod_prime(sw1[t][j], q_list[j]) for j in range(active)])

    return [out0, out1]


def rotate(ct, rotation_key, k=None):
    if k is None:
        k = rotation_key["k"]
    return keyswitch(automorphism_ct(ct, k), rotation_key["ksk"])


def conjugate(ct, conjugation_key):
    return rotate(ct, conjugation_key, -1)


# ------------------------------ ModRaise -------------------------------------

def _center_lift_int(a, q):
    a = Integer(a) % Integer(q)
    if a > Integer(q)//2:
        a -= Integer(q)
    return a


def _rns_reconstruct_coeff(residues, basis, centered=False):
    residues = [Integer(residues[i]) % Integer(basis[i]) for i in range(len(basis))]
    val = crt(residues, [Integer(q) for q in basis])
    if centered:
        M = prod([Integer(q) for q in basis])
        if val > M//2:
            val -= M
    return val


def _rns_reconstruct_poly(poly_shares, basis, centered=False):
    out = [0 for _ in range(POLDEG)]
    coeff_lists = [_dense_coeffs(poly_shares[j], POLDEG) for j in range(len(basis))]
    for i in range(POLDEG):
        out[i] = _rns_reconstruct_coeff([coeff_lists[j][i] for j in range(len(basis))], basis, centered=centered)
    return R(out) % firreducible


def _base_extend_poly(poly_shares, old_basis, new_basis, centered=False):
    poly = _rns_reconstruct_poly(poly_shares, old_basis, centered=centered)
    return [poly_mod_prime(poly, q) for q in new_basis]


def mod_raise(ct, old_basis=None, new_basis=None, centered=True):
    """Materialize a ciphertext on a larger RNS basis.

    For a q0-only bootstrap input, old_basis=[q0] and new_basis=q_list[:L].  This
    is the ModRaise/base-extension scaffold used before CoeffToSlot.
    """
    ct = normalize_ct(ct)
    active_old = len(ct[0])
    if old_basis is None:
        old_basis = [q_list[j] for j in range(active_old)]
    if new_basis is None:
        new_basis = [q_list[j] for j in range(L)]

    c0 = _base_extend_poly(ct[0], old_basis, new_basis, centered=centered)
    c1 = [_base_extend_poly(ct[1][i], old_basis, new_basis, centered=centered) for i in range(len(ct[1]))]
    return [c0, c1]


# -------------------------- linear transformations ---------------------------

def linear_transform(ct, diagonals, rotation_keys=None, encode_plain=False):
    """Evaluate sum_k diag_k * kappa_k(ct).

    diagonals can be {k: plaintext_poly_or_scalar} or [(k, pt), ...].
    k is an automorphism exponent.  k=1 is the identity and needs no key.
    """
    if rotation_keys is None:
        rotation_keys = {}
    try:
        items = list(diagonals.items())
    except Exception:
        items = list(diagonals)

    result = None
    for k, pt in items:
        kk = Integer(k) % (2*POLDEG)
        if kk == 1:
            term = copy_ct(ct)
        else:
            if kk not in rotation_keys:
                raise ValueError("missing rotation key for automorphism exponent %s" % kk)
            term = rotate(ct, rotation_keys[kk], kk)
        term = plaintext_multiply(term, pt, encode=encode_plain)
        result = term if result is None else add_ct(result, term)

    if result is None:
        return zero_ct(ct_active(ct), ct_rank(ct))
    return result


def coeff_to_slot(ct, coeff_to_slot_diagonals, rotation_keys=None):
    return linear_transform(ct, coeff_to_slot_diagonals, rotation_keys, encode_plain=False)


def slot_to_coeff(ct, slot_to_coeff_diagonals, rotation_keys=None):
    return linear_transform(ct, slot_to_coeff_diagonals, rotation_keys, encode_plain=False)


# ---------------------------- EvalMod / EvalExp ------------------------------

def ct_mul_relin(ct1, ct2, relin_key, do_rescale=False):
    level = min(ct_active(ct1), ct_active(ct2)) - 1
    prod = homomorphic_mult(normalize_ct(ct1), normalize_ct(ct2), level)
    out = relinearisation(prod, relin_key)
    if do_rescale and ct_active(out) > 1:
        out = rescale_ct(out)
    return normalize_ct(out)


def eval_poly(ct, coeffs, relin_key, do_rescale=False, encode_coeffs=False):
    """Naive polynomial evaluation sum_i coeffs[i]*ct^i.

    This is intentionally simple; after correctness debugging, replace with
    Paterson-Stockmeyer or the CHK+18 EvalExp recursive squaring.
    """
    ct = normalize_ct(ct)
    result = zero_ct(ct_active(ct), ct_rank(ct))
    if len(coeffs) == 0:
        return result

    result = add_plain(result, coeffs[0], encode=encode_coeffs)
    if len(coeffs) == 1:
        return result

    power = copy_ct(ct)
    for i in range(1, len(coeffs)):
        if coeffs[i] != 0:
            result = add_ct(result, plaintext_multiply(power, coeffs[i], encode=encode_coeffs))
        if i != len(coeffs)-1:
            power = ct_mul_relin(power, ct, relin_key, do_rescale=do_rescale)
    return result


def eval_exp_or_sine(ct, relin_key, coeffs=None, do_rescale=False):
    """EvalMod placeholder using a scaled-sine polynomial.

    Default: sin(x) Taylor x - x^3/6 + x^5/120.  This wires the arithmetic path;
    serious CKKS bootstrap experiments should replace coeffs by the chosen
    scaled-sine / EvalExp approximation and scaling convention.
    """
    if coeffs is None:
        coeffs = [0, 1, 0, -QQ(1)/QQ(6), 0, QQ(1)/QQ(120)]
    return eval_poly(ct, coeffs, relin_key, do_rescale=do_rescale, encode_coeffs=False)


def img_extract(ct, conjugation_key, scale_half=False):
    """Extract imaginary part algebraically as ct - Conj(ct).

    CKKS would use (ct-Conj(ct))/(2i).  The current real-toy encoding does not
    represent i, so this function returns the numerator by default.  Set
    scale_half=True to multiply by 1/2 in tests where that is desired.
    """
    out = sub_ct(ct, conjugate(ct, conjugation_key))
    if scale_half:
        out = scalar_multiply(out, QQ(1)/QQ(2), encode=False)
    return out

# -----------------------------------------------------------------------------
# CKKS bootstrap EvalMod / EvalExp prototype
# -----------------------------------------------------------------------------
#
# This is still not production CKKS EvalMod.  It is an "almost bootstrap" bridge:
#
#     - works on ciphertexts whose slots contain coefficient values t_j
#     - evaluates a scaled sine approximation to [t_j]_{q0}
#     - uses ciphertext-ciphertext multiplications and relinearization
#
# Mathematical target:
#
#     EvalMod(x) ≈ q0/(2*pi) * sin(2*pi*x/q0)
#
# Around x = q0*I + m with |m| << q0, this approximates m.
#
# The coefficients contain powers of pi, so we encode slot constants as CKKS
# plaintext polynomials rather than trying to represent them exactly modulo q.
# -----------------------------------------------------------------------------


def _ckks_slot_constant_plain(alpha, scale=p_enc, method="fft"):
    """Encode the same scalar alpha in all CKKS slots as a plaintext polynomial."""
    n = ckks_num_slots()
    return ckks_encode([CCK(alpha) for _ in range(n)], scale, method=method)


def plaintext_multiply_slot_constant(ct,
                                     alpha,
                                     const_scale=p_enc,
                                     method="fft",
                                     divide_by_const_scale=True):
    """Multiply ciphertext slot-wise by an approximate CKKS constant.

    The constant is encoded as a plaintext polynomial using ckks_encode().
    Because the plaintext is scaled by const_scale, we optionally divide the
    ciphertext by const_scale afterwards.  This mirrors the prototype correction
    used in linear_transform.sage for diagonal plaintexts.
    """
    pt = _ckks_slot_constant_plain(alpha, const_scale, method)
    out = plaintext_multiply(ct, pt, encode=False)

    if divide_by_const_scale:
        out = scalar_multiply(out, QQ(1)/QQ(const_scale), encode=False)

    return out


def eval_mod_sine_coeffs(q0=None, degree=5):
    """Return approximate coefficients for q0/(2*pi)*sin(2*pi*x/q0).

    degree must be odd.  For degree=5:

        x
        - (2*pi)^2/(3! q0^2) x^3
        + (2*pi)^4/(5! q0^4) x^5

    Coefficients are returned as high-precision complex/real numbers and should
    be encoded as CKKS plaintext constants, not as exact modular rationals.
    """
    if q0 is None:
        q0 = q_list[0]

    if degree % 2 == 0:
        raise ValueError("degree must be odd")

    RR_high = RealField(160)
    pi_high = RR_high.pi()
    q0_high = RR_high(q0)

    coeffs = {}
    max_k = (degree - 1)//2

    for k in range(max_k + 1):
        power = 2*k + 1
        sign = -1 if (k % 2 == 1) else 1

        # q0/(2*pi) * (2*pi/q0)^power / power!
        c = sign * ( (2*pi_high)**(power-1) ) / ( factorial(power) * (q0_high**(power-1)) )
        coeffs[power] = CCK(c)

    return coeffs


def eval_mod_sine_prototype(ct,
                            relin_key,
                            q0=None,
                            degree=5,
                            const_scale=p_enc,
                            method="fft",
                            do_rescale=False):
    """Prototype slot-wise EvalMod using scaled sine Taylor approximation.

    Input:
        ct encrypts slots approximately containing t_j = q0*I_j + m_j.

    Output:
        ciphertext encrypting slots approximately containing m_j.

    This is intended for small-N correctness experiments.  It does not yet
    implement minimax approximation, double-angle EvalExp, precise scale-chain
    scheduling, or production CKKS rescaling.
    """
    if relin_key is None:
        raise ValueError("eval_mod_sine_prototype requires a relinearization key")

    if q0 is None:
        q0 = q_list[0]

    coeffs = eval_mod_sine_coeffs(q0=q0, degree=degree)

    # result = c1*x
    result = plaintext_multiply_slot_constant(
        ct,
        coeffs[1],
        const_scale=const_scale,
        method=method,
        divide_by_const_scale=True
    )

    if degree == 1:
        return result

    # x2 = x*x
    x2 = ct_mul_relin(ct, ct, relin_key, do_rescale=do_rescale)

    # power = x
    power = copy_ct(ct)

    for d in range(3, degree + 1, 2):
        # power = x^(d-2) * x^2 = x^d
        power = ct_mul_relin(power, x2, relin_key, do_rescale=do_rescale)

        term = plaintext_multiply_slot_constant(
            power,
            coeffs[d],
            const_scale=const_scale,
            method=method,
            divide_by_const_scale=True
        )

        result = add_ct(result, term)

    return result

# -----------------------------------------------------------------------------
# OpenFHE-style normalized EvalMod prototype
# -----------------------------------------------------------------------------
#
# OpenFHE-style idea:
#
#   After ModRaise, the ciphertext decrypts to
#
#       t = m + q0*I
#
#   Instead of evaluating
#
#       q0/(2*pi) * sin(2*pi*t/q0)
#
#   directly on huge t, we first normalize:
#
#       y = t / (K*q0)
#         = I/K + m/(K*q0)
#
#   Then evaluate
#
#       post/(2*pi) * sin(2*pi*K*y)
#
#   Since
#
#       sin(2*pi*K*y)
#         = sin(2*pi*I + 2*pi*m/q0)
#         ≈ 2*pi*m/q0,
#
#   the result is approximately
#
#       post * m/q0.
#
#   If post ≈ q0/p_enc, then the result is approximately m/p_enc,
#   which matches the decoded CKKS coefficient scale convention in this Sage
#   prototype.
# -----------------------------------------------------------------------------


def eval_mod_sine_openfhe_coeffs(post, K=1, degree=7):
    """Coefficients of post/(2*pi) * sin(2*pi*K*x).

    For degree 7:

        post*K*x
        - post*(2*pi)^2*K^3/3! * x^3
        + post*(2*pi)^4*K^5/5! * x^5
        - post*(2*pi)^6*K^7/7! * x^7

    These coefficients are numerically meaningful after normalization, unlike
    the direct q0-sine coefficients, whose nonlinear terms contain powers of
    1/q0.
    """
    if degree % 2 == 0:
        raise ValueError("degree must be odd")

    RR_high = RealField(160)
    pi_high = RR_high.pi()
    post_high = RR_high(post)
    K_high = RR_high(K)

    coeffs = {}
    max_k = (degree - 1)//2

    for r in range(max_k + 1):
        power = 2*r + 1
        sign = -1 if (r % 2 == 1) else 1

        # post/(2*pi) * (2*pi*K)^power / power!
        c = sign * post_high * ((2*pi_high*K_high)**(power - 1)) * K_high / factorial(power)

        coeffs[power] = CCK(c)

    return coeffs


def eval_mod_sine_openfhe_prototype(ct,
                                    relin_key,
                                    post=None,
                                    K=1,
                                    degree=7,
                                    const_scale=p_enc,
                                    method="fft",
                                    do_rescale=False):
    """OpenFHE-style normalized sine EvalMod prototype.

    Input:
        ct encrypts normalized coefficient slots

            y_j = t_j / (K*q0)

        where t_j = m_j + q0*I_j.

    Output:
        approximate encryption of

            post/(2*pi) * sin(2*pi*K*y_j).

    With post ≈ q0/p_enc, this approximates m_j/p_enc in decoded scale.

    This is still a prototype:
        - Taylor sine, not Chebyshev/minimax.
        - No true CKKS scale-chain management.
        - No double-angle EvalExp yet.
    """
    if relin_key is None:
        raise ValueError("eval_mod_sine_openfhe_prototype requires a relinearization key")

    if post is None:
        post = RR(q_list[0]) / RR(p_enc)

    coeffs = eval_mod_sine_openfhe_coeffs(
        post=post,
        K=K,
        degree=degree
    )

    # result = c1*x
    result = plaintext_multiply_slot_constant(
        ct,
        coeffs[1],
        const_scale=const_scale,
        method=method,
        divide_by_const_scale=True
    )

    if degree == 1:
        return result

    # x2 = x*x
    x2 = ct_mul_relin(
        ct,
        ct,
        relin_key,
        do_rescale=do_rescale
    )

    # power = x
    power = copy_ct(ct)

    for d in range(3, degree + 1, 2):
        # power = x^(d-2) * x^2 = x^d
        power = ct_mul_relin(
            power,
            x2,
            relin_key,
            do_rescale=do_rescale
        )

        term = plaintext_multiply_slot_constant(
            power,
            coeffs[d],
            const_scale=const_scale,
            method=method,
            divide_by_const_scale=True
        )

        result = add_ct(result, term)

    return result

"""
def eval_exp_or_sine_ckks_prototype(ct,
                                    relin_key,
                                    q0=None,
                                    degree=5,
                                    const_scale=p_enc,
                                    method="fft",
                                    do_rescale=False):
    #Compatibility wrapper replacing the old eval_exp_or_sine() placeholder.

    #For now this calls the scaled-sine EvalMod prototype.
    #Later, this can be replaced by the CHK+18-style complex EvalExp with
    #recursive squaring and ImgExt.
    
    return eval_mod_sine_prototype(
        ct,
        relin_key,
        q0=q0,
        degree=degree,
        const_scale=const_scale,
        method=method,
        do_rescale=do_rescale
    )
"""
def eval_exp_or_sine_ckks_prototype(ct,
                                    relin_key,
                                    q0=None,
                                    degree=7,
                                    const_scale=p_enc,
                                    method="fft",
                                    do_rescale=False,
                                    K=1,
                                    post=None,
                                    normalized=True):
    """Compatibility wrapper for the bootstrap EvalMod/EvalExp stage.

    Default path:
        OpenFHE-style normalized EvalMod:

            post/(2*pi) * sin(2*pi*K*x)

        where x is expected to already be normalized:

            x = (m + q0*I)/(K*q0).

    Fallback path:
        old direct EvalMod:

            q0/(2*pi) * sin(2*pi*x/q0)

        This is mostly useful only for toy experiments with small q0.
    """
    if q0 is None:
        q0 = q_list[0]

    if normalized:
        if post is None:
            post = RR(q0) / RR(p_enc)

        return eval_mod_sine_openfhe_prototype(
            ct,
            relin_key,
            post=post,
            K=K,
            degree=degree,
            const_scale=const_scale,
            method=method,
            do_rescale=do_rescale
        )

    return eval_mod_sine_prototype(
        ct,
        relin_key,
        q0=q0,
        degree=degree,
        const_scale=const_scale,
        method=method,
        do_rescale=do_rescale
    )

# -----------------------------------------------------------------------------
# Two-branch almost-CKKS bootstrap wrapper
# -----------------------------------------------------------------------------
#
# This wrapper assumes linear_transform.sage has been loaded and provides:
#
#     coeff_to_slot_bootstrap_two_ct()
#     slot_to_coeff_bootstrap_two_ct()
#
# Pipeline:
#
#     ct
#       -> ModRaise
#       -> CoeffToSlot: ct_lo, ct_hi
#       -> EvalMod on ct_lo and ct_hi
#       -> SlotToCoeff merge
#
# This is intentionally called "almost" because scale/level management and the
# EvalMod approximation are still prototype-level.
# -----------------------------------------------------------------------------

"""
def bootstrap_almost_ckks(ct, boot_keys, boot_params):
    #Almost-valid traditional MLWE-CKKS bootstrap scaffold.

    #Required boot_params:
        #"c2s_boot_diags"
        #"s2c_boot_diags"

    #Optional boot_params:
        #old_basis
        #new_basis
       # centered_modraise
        #divide_by_diag_scale
       # q0
        #eval_degree
        #const_scale
        #method
        #do_rescale
       # run_evalmod

   # If run_evalmod=False, the function only tests:
       # ModRaise -> C2S -> S2C
    
    if boot_params is None:
        raise ValueError("bootstrap_almost_ckks requires boot_params")

    if "c2s_boot_diags" not in boot_params:
        raise ValueError("missing boot_params['c2s_boot_diags']")

    if "s2c_boot_diags" not in boot_params:
        raise ValueError("missing boot_params['s2c_boot_diags']")

    relin_key = boot_keys.get("relin_key", None)
    run_evalmod = boot_params.get("run_evalmod", True)

    if run_evalmod and relin_key is None:
        raise ValueError("EvalMod requires boot_keys['relin_key']; regenerate boot keys with include_relin=True")

    # 1. ModRaise / base extension.
    ct_top = mod_raise(
        ct,
        old_basis=boot_params.get("old_basis", None),
        new_basis=boot_params.get("new_basis", None),
        centered=boot_params.get("centered_modraise", True)
    )

    # 2. True bootstrap CoeffToSlot:
    #    one ciphertext -> two ciphertexts containing coefficient halves.
    ct_lo, ct_hi = coeff_to_slot_bootstrap_two_ct(
        ct_top,
        boot_params["c2s_boot_diags"],
        boot_keys,
        divide_by_diag_scale=boot_params.get("divide_by_diag_scale", True)
    )

    # 3. EvalMod / EvalExp prototype on both halves.
    if run_evalmod:
        ct_lo = eval_exp_or_sine_ckks_prototype(
            ct_lo,
            relin_key,
            q0=boot_params.get("q0", q_list[0]),
            degree=boot_params.get("eval_degree", 5),
            const_scale=boot_params.get("const_scale", p_enc),
            method=boot_params.get("method", "fft"),
            do_rescale=boot_params.get("do_rescale", False)
        )

        ct_hi = eval_exp_or_sine_ckks_prototype(
            ct_hi,
            relin_key,
            q0=boot_params.get("q0", q_list[0]),
            degree=boot_params.get("eval_degree", 5),
            const_scale=boot_params.get("const_scale", p_enc),
            method=boot_params.get("method", "fft"),
            do_rescale=boot_params.get("do_rescale", False)
        )

    # 4. True bootstrap SlotToCoeff:
    #    two coefficient-slot ciphertexts -> one native polynomial ciphertext.
    ct_out = slot_to_coeff_bootstrap_two_ct(
        ct_lo,
        ct_hi,
        boot_params["s2c_boot_diags"],
        boot_keys,
        divide_by_diag_scale=boot_params.get("divide_by_diag_scale", True)
    )

    return ct_out
"""
def bootstrap_almost_ckks(ct, boot_keys, boot_params):
    """Almost-valid traditional MLWE-CKKS bootstrap scaffold.

    Required boot_params:
        "c2s_boot_diags"
        "s2c_boot_diags"

    Optional boot_params:
        old_basis
        new_basis
        centered_modraise
        do_modraise
        normalize_before_c2s
        divide_by_diag_scale
        q0
        K
        post
        eval_degree
        const_scale
        method
        do_rescale
        run_evalmod
        normalized_evalmod

    Pipeline:

        ct
          -> optional ModRaise
          -> optional OpenFHE-style normalization
          -> true two-ciphertext CoeffToSlot
          -> optional normalized EvalMod prototype
          -> true two-ciphertext SlotToCoeff
    """
    if boot_params is None:
        raise ValueError("bootstrap_almost_ckks requires boot_params")

    if "c2s_boot_diags" not in boot_params:
        raise ValueError("missing boot_params['c2s_boot_diags']")

    if "s2c_boot_diags" not in boot_params:
        raise ValueError("missing boot_params['s2c_boot_diags']")

    relin_key = boot_keys.get("relin_key", None)
    run_evalmod = boot_params.get("run_evalmod", True)

    if run_evalmod and relin_key is None:
        raise ValueError("EvalMod requires boot_keys['relin_key']; regenerate boot keys with include_relin=True")

    # -------------------------------------------------------------------------
    # 1. ModRaise / base extension.
    # -------------------------------------------------------------------------
    if boot_params.get("do_modraise", True):
        ct_top = mod_raise(
            ct,
            old_basis=boot_params.get("old_basis", None),
            new_basis=boot_params.get("new_basis", None),
            centered=boot_params.get("centered_modraise", True)
        )
    else:
        ct_top = copy_ct(ct)

    # -------------------------------------------------------------------------
    # 1b. OpenFHE-style normalization before CoeffToSlot.
    #
    # After ModRaise:
    #
    #     Dec(ct_top) = m + q0*I.
    #
    # We want the coefficient slots after CoeffToSlot to contain roughly:
    #
    #     y = (m + q0*I)/(K*q0).
    #
    # Since decryption/decoding in this Sage prototype divides by p_enc, we
    # multiply the plaintext polynomial by:
    #
    #     beta = p_enc/(K*q0).
    #
    # This should only be enabled when EvalMod is enabled.  For the no-EvalMod
    # C2S/S2C preservation test, keep normalize_before_c2s=False.
    # -------------------------------------------------------------------------
    '''
    if boot_params.get("normalize_before_c2s", False):
        q0_norm = boot_params.get("q0", q_list[0])
        K_norm = boot_params.get("K", 1)

        beta = QQ(p_enc) / QQ(K_norm * q0_norm)

        ct_top = scalar_multiply(
            ct_top,
            beta,
            encode=False
        )
    '''

    if boot_params.get("normalize_before_c2s", False):
        q0_norm = boot_params.get("q0", q_list[0])
        K_norm = boot_params.get("K", 1)

    # beta is a real CKKS plaintext constant, not a modular rational.
    #
    # Do NOT use scalar_multiply(ct_top, beta), because beta has denominator q0,
    # and q0 is not invertible modulo the q0 RNS limb.
        beta = RR(p_enc) / RR(K_norm * q0_norm)

        ct_top = plaintext_multiply_slot_constant(
            ct_top,
            beta,
            const_scale=boot_params.get("const_scale", p_enc),
            method=boot_params.get("method", "fft"),
            divide_by_const_scale=True
        )

    # -------------------------------------------------------------------------
    # 2. True bootstrap CoeffToSlot:
    #
    #     one ciphertext -> two ciphertexts
    #
    #     ct_lo encrypts coefficient slots:
    #         (t0, ..., t_{N/2-1})
    #
    #     ct_hi encrypts coefficient slots:
    #         (t_{N/2}, ..., t_{N-1})
    # -------------------------------------------------------------------------
    ct_lo, ct_hi = coeff_to_slot_bootstrap_two_ct(
        ct_top,
        boot_params["c2s_boot_diags"],
        boot_keys,
        divide_by_diag_scale=boot_params.get("divide_by_diag_scale", True)
    )

    # -------------------------------------------------------------------------
    # 3. EvalMod / EvalExp prototype on both halves.
    # -------------------------------------------------------------------------
    if run_evalmod:
        q0_eval = boot_params.get("q0", q_list[0])
        post_eval = boot_params.get("post", RR(q0_eval) / RR(p_enc))

        ct_lo = eval_exp_or_sine_ckks_prototype(
            ct_lo,
            relin_key,
            q0=q0_eval,
            degree=boot_params.get("eval_degree", 7),
            const_scale=boot_params.get("const_scale", p_enc),
            method=boot_params.get("method", "fft"),
            do_rescale=boot_params.get("do_rescale", False),
            K=boot_params.get("K", 1),
            post=post_eval,
            normalized=boot_params.get("normalized_evalmod", True)
        )

        ct_hi = eval_exp_or_sine_ckks_prototype(
            ct_hi,
            relin_key,
            q0=q0_eval,
            degree=boot_params.get("eval_degree", 7),
            const_scale=boot_params.get("const_scale", p_enc),
            method=boot_params.get("method", "fft"),
            do_rescale=boot_params.get("do_rescale", False),
            K=boot_params.get("K", 1),
            post=post_eval,
            normalized=boot_params.get("normalized_evalmod", True)
        )

    # -------------------------------------------------------------------------
    # 4. True bootstrap SlotToCoeff:
    #
    #     two coefficient-slot ciphertexts -> one native polynomial ciphertext.
    # -------------------------------------------------------------------------
    ct_out = slot_to_coeff_bootstrap_two_ct(
        ct_lo,
        ct_hi,
        boot_params["s2c_boot_diags"],
        boot_keys,
        divide_by_diag_scale=boot_params.get("divide_by_diag_scale", True)
    )

    return ct_out

# ------------------------------- bootstrap -----------------------------------

def bootstrap(ct, boot_keys, boot_params=None):
    """Traditional CKKS-style ModHE bootstrap scaffold.

    boot_params keys:
      old_basis, new_basis, centered_modraise
      coeff_to_slot_diagonals   default {1: 1}
      slot_to_coeff_diagonals   default {1: 1}
      eval_coeffs               default sine Taylor placeholder
      do_rescale                default False
      extract_imag              default True
      scale_half_img            default False
    """
    if boot_params is None:
        boot_params = {}

    rot_keys = boot_keys.get("rot_keys", {})
    conj_key = boot_keys.get("conj_key", None)
    relin_key = boot_keys.get("relin_key", None)
    if relin_key is None:
        raise ValueError("bootstrap requires boot_keys['relin_key']")

    ct1 = mod_raise(ct,
                    old_basis=boot_params.get("old_basis", None),
                    new_basis=boot_params.get("new_basis", None),
                    centered=boot_params.get("centered_modraise", True))

    ct2 = coeff_to_slot(ct1,
                        boot_params.get("coeff_to_slot_diagonals", {1: R(1)}),
                        rot_keys)

    ct3 = eval_exp_or_sine(ct2,
                           relin_key,
                           coeffs=boot_params.get("eval_coeffs", None),
                           do_rescale=boot_params.get("do_rescale", False))

    if boot_params.get("extract_imag", True):
        if conj_key is None:
            raise ValueError("bootstrap requires boot_keys['conj_key'] when extract_imag=True")
        ct4 = img_extract(ct3, conj_key, scale_half=boot_params.get("scale_half_img", False))
    else:
        ct4 = ct3

    ct5 = slot_to_coeff(ct4,
                        boot_params.get("slot_to_coeff_diagonals", {1: R(1)}),
                        rot_keys)
    return ct5
