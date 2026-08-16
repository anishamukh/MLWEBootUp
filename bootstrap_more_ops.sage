# -----------------------------------------------------------------------------
# eval_mod_double_angle.sage
#
# HEAAN-style double-angle EvalMod for the ModHE real-encoding prototype.
#
# This is the ModHE analogue of HEAAN's evalExpAndEqual + exp2piAndEqual
# (Scheme.cpp).  Because the ModHE toy encoder is REAL-valued (it cannot
# represent multiplication by i / has no imult), we cannot square a complex
# exponential the way HEAAN does.  Instead we carry a real (sin, cos) PAIR and
# apply the real double-angle identities
#
#       sin(2*theta) = 2 * sin(theta) * cos(theta)
#       cos(2*theta) = cos(theta)^2 - sin(theta)^2
#
# logT times.  This is the faithful real-domain counterpart of HEAAN's
#
#       for i in range(logI + logT): squareAndEqual(cipher); reScale(...)
#
# doubling loop.  Argument reduction by 2^logT keeps the base Taylor series of
# sin/cos accurate, exactly as HEAAN reduces before exp2piAndEqual.
#
# Mathematical target (identical to eval_mod_sine_openfhe_prototype):
#
#       EvalMod(t) ~= post/(2*pi) * sin(2*pi*K*y),   with  y = t/(K*q0),
#                     t = m + q0*I,  post ~= q0/p_enc  =>  result ~= m/p_enc.
#
# Load AFTER the standard ModHE stack, e.g.:
#     load("constants.sage"); load("auxiliary_functions.sage")
#     load("KeyGen.sage");    load("Encrypt-Decrypt.sage")
#     load("encode.sage");    load("encode-encrypt.sage")
#     load("HMult.sage");     load("Relin.sage")
#     load("bootstrap_keys.sage"); load("bootstrap_ops.sage")
#     load("linear_transform.sage")
#     load("eval_mod_double_angle.sage")
#
# Dependencies used from the existing code base (all in bootstrap_ops.sage
# unless noted):
#     normalize_ct, copy_ct, zero_ct, ct_active, ct_rank
#     add_ct, sub_ct, scalar_multiply, add_plain
#     ct_mul_relin
#     plaintext_multiply_slot_constant, _ckks_slot_constant_plain
#     encryption_ckks / decryption_ckks (encode-encrypt.sage)
#     public_private_keygen, relin_keygen (KeyGen.sage / Relin.sage)
#     CCK, RRK, p_enc, q_list, POLDEG, L, ckks_num_slots
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# High-precision constants
# -----------------------------------------------------------------------------

def _RRh():
    return RealField(160)


def _two_pi_high():
    return _RRh()(2) * _RRh().pi()


# -----------------------------------------------------------------------------
# Base Taylor coefficients for sin/cos of a reduced argument
# -----------------------------------------------------------------------------

def sincos_small_taylor_coeffs(alpha, degree=7):
    """Taylor coefficients of sin(alpha*x) and cos(alpha*x) about x=0.

    Returns (sin_coeffs, cos_coeffs), each a dict {power: CCK constant}.

        sin(alpha*x) = sum_{k>=0} (-1)^k alpha^(2k+1) x^(2k+1)/(2k+1)!
        cos(alpha*x) = sum_{k>=0} (-1)^k alpha^(2k)   x^(2k)  /(2k)!

    'alpha' is the angular scale multiplying x inside sin/cos.  'degree' is the
    highest power kept (both series are truncated at this total degree).
    """
    RRh = _RRh()
    a = RRh(alpha)

    sin_c = {}
    cos_c = {}
    for pwr in range(degree + 1):
        k = pwr // 2
        sign = -1 if (k % 2 == 1) else 1
        term = sign * (a ** pwr) / factorial(pwr)
        if pwr % 2 == 1:
            sin_c[pwr] = CCK(term)
        else:
            cos_c[pwr] = CCK(term)
    return sin_c, cos_c


# -----------------------------------------------------------------------------
# Slot-constant ciphertext helper
# -----------------------------------------------------------------------------

def slot_constant_ct(like_ct, alpha, const_scale=p_enc, method="fft"):
    """Return a ciphertext-shaped object encoding constant 'alpha' in all slots.

    Built as add_plain onto a zero ciphertext, so it lives at the same active
    level (number of q-shares) as 'like_ct'.  The plaintext is scaled by
    const_scale; callers that need the un-scaled constant must divide later.
    """
    z = zero_ct(ct_active(like_ct), ct_rank(like_ct))
    pt = _ckks_slot_constant_plain(alpha, const_scale, method)
    return add_plain(z, pt, encode=False)


# -----------------------------------------------------------------------------
# Sparse even/odd polynomial evaluation given x and x^2
# -----------------------------------------------------------------------------

def _eval_sparse_poly(ct, x2, coeffs, relin_key,
                      const_scale=p_enc, method="fft", do_rescale=True):
    """Evaluate sum_p coeffs[p] * x^p, where 'ct' encodes x and 'x2' encodes x^2.

    coeffs is a dict {power: CCK constant}.  Constant (power 0) and both even
    and odd powers are supported.  Powers above 2 are built by repeatedly
    multiplying the running odd/even power ciphertext by x2.

    Scale/level note: each coeffs[p]*x^p term is formed with
    plaintext_multiply_slot_constant(divide_by_const_scale=True), which is a
    pure scalar multiply and does NOT consume a level.  Levels are consumed only
    by ct_mul_relin (when do_rescale=True) while building the powers.  All terms
    are therefore added at the level of their respective power; because higher
    powers sit at deeper levels, we add each term into 'result' by re-basing to
    the lowest active level via add_ct (which truncates to min active shares).
    """
    if len(coeffs) == 0:
        return zero_ct(ct_active(ct), ct_rank(ct))

    powers = sorted(coeffs.keys())
    result = None

    # p = 0 (constant)
    if 0 in coeffs:
        result = slot_constant_ct(ct, coeffs[0], const_scale, method)

    # p = 1
    if 1 in coeffs:
        term = plaintext_multiply_slot_constant(
            ct, coeffs[1], const_scale=const_scale, method=method,
            divide_by_const_scale=True)
        result = term if result is None else add_ct(result, term)

    # p = 2
    if 2 in coeffs:
        term = plaintext_multiply_slot_constant(
            x2, coeffs[2], const_scale=const_scale, method=method,
            divide_by_const_scale=True)
        result = term if result is None else add_ct(result, term)

    # p >= 3: walk up in unit steps, keeping separate odd/even running powers
    max_p = powers[-1]
    if max_p >= 3:
        power_odd = copy_ct(ct)    # x^1
        power_even = copy_ct(x2)   # x^2
        for pwr in range(3, max_p + 1):
            if pwr % 2 == 1:
                power_odd = ct_mul_relin(power_odd, x2, relin_key, do_rescale=do_rescale)
                cur = power_odd
            else:
                power_even = ct_mul_relin(power_even, x2, relin_key, do_rescale=do_rescale)
                cur = power_even
            if pwr in coeffs:
                term = plaintext_multiply_slot_constant(
                    cur, coeffs[pwr], const_scale=const_scale, method=method,
                    divide_by_const_scale=True)
                result = add_ct(result, term)

    if result is None:
        return zero_ct(ct_active(ct), ct_rank(ct))
    return result


# -----------------------------------------------------------------------------
# Level-budget estimate
# -----------------------------------------------------------------------------

def eval_mod_double_angle_level_cost(logT, base_degree=7, do_rescale=True):
    """Estimate the number of q-levels consumed by eval_mod_double_angle.

    Returns a dict with a breakdown.  Only ct_mul_relin(do_rescale=True) steps
    consume a level; scalar/plaintext-constant multiplies do not.

    Cost model:
        - x2 = x*x                                  : 1 level
        - building base sin/cos powers up to
          base_degree from x2                       : (max_p - 2) levels once,
          where max_p = base_degree (odd) or base_degree (shared x2 reused,
          so the incremental multiplies are floor((base_degree-1)/2)*... );
          we bound it simply below.
        - each doubling iteration does:
              S*C, C*C, S*S                          : 3 levels
          times logT.
    This is an UPPER-BOUND-ish estimate meant to warn about chain underflow.
    """
    if not do_rescale:
        return {"total": 0, "note": "do_rescale=False consumes no levels but "
                                    "will overflow ciphertext scale; not for real runs"}

    # x^2
    lvl_x2 = 1
    # incremental multiplies to reach base_degree using shared x2:
    #   odd chain: x^3,x^5,...  even chain: x^4,x^6,...
    #   each new power is one ct_mul_relin with x2
    incr = 0
    for pwr in range(3, base_degree + 1):
        incr += 1
    lvl_base = incr

    # doubling loop
    lvl_double = 3 * logT

    total = lvl_x2 + lvl_base + lvl_double
    return {
        "x2": lvl_x2,
        "base_powers": lvl_base,
        "doubling": lvl_double,
        "total": total,
    }


def warn_if_chain_too_short(start_active, logT, base_degree=7, do_rescale=True):
    """Print a warning if the RNS chain is likely too short for EvalMod.

    'start_active' is the number of active q-shares of the ciphertext entering
    EvalMod (e.g. ct_active(ct)).  We need at least (level cost + 1) shares so
    that at least q0 remains at the end.
    """
    cost = eval_mod_double_angle_level_cost(logT, base_degree, do_rescale)
    need = cost["total"] + 1
    print("[evalmod] level budget: start_active=%d, estimated cost=%d, need>=%d"
          % (start_active, cost["total"], need))
    if start_active < need:
        print("[evalmod] WARNING: RNS chain likely too short. "
              "Reduce logT/base_degree, or raise L in constants.sage. "
              "(breakdown: %s)" % cost)
        return False
    return True


# -----------------------------------------------------------------------------
# Double-angle EvalMod
# -----------------------------------------------------------------------------

def eval_mod_double_angle(ct,
                          relin_key,
                          post=None,
                          K=1,
                          logT=3,
                          base_degree=7,
                          const_scale=p_enc,
                          method="fft",
                          do_rescale=True,
                          check_budget=True,
                          verbose=False):
    """HEAAN-style double-angle EvalMod on a real sin/cos pair.

    Input:
        ct encrypts normalized slots  y_j = t_j/(K*q0),  t_j = m_j + q0*I_j.
        relin_key from relin_keygen(key).

    Parameters:
        post        scaling applied at the end; default q0/p_enc so the decoded
                    output matches the message scale m/p_enc.
        K           range factor used in the normalization y = t/(K*q0).
        logT        number of double-angle doublings (HEAAN's logT).  The base
                    argument is theta0 = (2*pi*K*y)/2^logT.
        base_degree highest power in the base sin/cos Taylor series (odd; both
                    series truncate at this total degree).
        do_rescale  rescale after each ciphertext-ciphertext multiply (needed
                    for a realistic level chain; keep True).

    Output:
        approximate encryption of  post/(2*pi) * sin(2*pi*K*y_j)  ~=  post*m_j/q0.

    Structure (mirrors evalExpAndEqual):
        1. Reduce: theta0 = alpha0 * y with alpha0 = 2*pi*K / 2^logT.
        2. Base: S = sin(theta0), C = cos(theta0) via short Taylor series.
        3. Double logT times:  S <- 2SC,  C <- C^2 - S^2.
        4. Scale by post/(2*pi).
    """
    if relin_key is None:
        raise ValueError("eval_mod_double_angle requires a relinearization key")
    if post is None:
        post = RRK(q_list[0]) / RRK(p_enc)

    ct = normalize_ct(ct)

    if check_budget:
        warn_if_chain_too_short(ct_active(ct), logT, base_degree, do_rescale)

    RRh = _RRh()
    alpha0 = _two_pi_high() * RRh(K) / RRh(2 ** logT)

    sin_c, cos_c = sincos_small_taylor_coeffs(alpha0, degree=base_degree)

    # shared x^2 for both series
    x2 = ct_mul_relin(ct, ct, relin_key, do_rescale=do_rescale)

    S = _eval_sparse_poly(ct, x2, sin_c, relin_key,
                          const_scale=const_scale, method=method, do_rescale=do_rescale)
    C = _eval_sparse_poly(ct, x2, cos_c, relin_key,
                          const_scale=const_scale, method=method, do_rescale=do_rescale)

    if verbose:
        print("[evalmod] base sin/cos built; S active=%d C active=%d; doubling %d times"
              % (ct_active(S), ct_active(C), logT))

    # Double-angle doubling loop (HEAAN squareAndEqual analogue)
    for i in range(logT):
        SC = ct_mul_relin(S, C, relin_key, do_rescale=do_rescale)
        newS = scalar_multiply(SC, 2, encode=False)

        CC = ct_mul_relin(C, C, relin_key, do_rescale=do_rescale)
        SS = ct_mul_relin(S, S, relin_key, do_rescale=do_rescale)
        newC = sub_ct(CC, SS)

        S, C = newS, newC
        if verbose:
            print("[evalmod] doubling %d/%d done; S active=%d" % (i + 1, logT, ct_active(S)))

    # final scale: post/(2*pi) * S
    scale_const = RRK(post) / RRK(_two_pi_high())
    result = plaintext_multiply_slot_constant(
        S, CCK(scale_const), const_scale=const_scale, method=method,
        divide_by_const_scale=True)

    return result


# -----------------------------------------------------------------------------
# Cleartext reference (for validating the encrypted result)
# -----------------------------------------------------------------------------

def eval_mod_double_angle_clear(y, post=None, K=1):
    """Cleartext target value post/(2*pi) * sin(2*pi*K*y) for a slot value y."""
    if post is None:
        post = RRK(q_list[0]) / RRK(p_enc)
    RRh = _RRh()
    val = RRh(post) / _two_pi_high() * (_two_pi_high() * RRh(K) * RRh(y)).sin()
    return CCK(val)


# -----------------------------------------------------------------------------
# Self-test: encrypted double-angle EvalMod vs cleartext sine
# -----------------------------------------------------------------------------

def eval_mod_double_angle_selftest(y_slots=None,
                                   K=1,
                                   logT=3,
                                   base_degree=7,
                                   post=None,
                                   do_rescale=True,
                                   verbose=True):
    """Encrypt normalized slot values y_j, run double-angle EvalMod, and compare
    against cleartext post/(2*pi)*sin(2*pi*K*y_j).

    Intended for SMALL POLDEG (e.g. 8 or 16).  Set POLDEG in constants.sage
    accordingly before loading.

    Because y_j = I_j/K + m_j/(K*q0), a realistic test uses small |y| near an
    integer multiple of 1/K so that sin(2*pi*K*y) ~= 2*pi*(fractional part).
    For a pure functional test of the sine evaluation, small |y_j| in
    (-1/(2K), 1/(2K)) is the cleanest.
    """
    nslots = ckks_num_slots()

    if y_slots is None:
        # small arguments where sin is well approximated; a few sign/magnitude cases
        base = [0.00, 0.10, -0.15, 0.20][:nslots]
        y_slots = base + [0.0] * (nslots - len(base))
    else:
        y_slots = list(y_slots) + [0.0] * (nslots - len(y_slots))

    if post is None:
        post = RRK(q_list[0]) / RRK(p_enc)

    print("=" * 70)
    print("eval_mod_double_angle self-test")
    print("POLDEG=%d  slots=%d  L=%d  K=%d  logT=%d  base_degree=%d"
          % (POLDEG, nslots, L, K, logT, base_degree))
    print("post = q0/p_enc =", RRK(post))
    print("=" * 70)

    warn_if_chain_too_short(L, logT, base_degree, do_rescale)

    # keys
    key = public_private_keygen()
    relin_key = relin_keygen(key)

    # encrypt the normalized y-slots at full chain
    ct = encryption_ckks([CCK(v) for v in y_slots], key, p_enc, method="fft")

    # run EvalMod
    ct_out = eval_mod_double_angle(
        ct, relin_key,
        post=post, K=K, logT=logT, base_degree=base_degree,
        do_rescale=do_rescale, check_budget=False, verbose=verbose)

    # decrypt/decode
    out_slots = decryption_ckks(ct_out, key, 0, 0, p_enc, nslots)

    print("\nslot :         y_j          EvalMod(enc)             cleartext           abs_err")
    max_err = RRK(0)
    for j in range(min(nslots, 8)):
        target = eval_mod_double_angle_clear(y_slots[j], post=post, K=K)
        got = CCK(out_slots[j])
        err = abs(got - target)
        max_err = max(max_err, RRK(err))
        print("%4d : %12.6f   %20s   %20s   %10.3e"
              % (j, float(y_slots[j]), _fmt_c(got), _fmt_c(target), float(err)))

    print("\nmax abs error over reported slots =", float(max_err))
    print("=" * 70)
    return ct_out, out_slots, float(max_err)


def _fmt_c(z):
    z = CCK(z)
    return "%+.6f%+.6fi" % (float(z.real()), float(z.imag()))


# -----------------------------------------------------------------------------
# Optional: comparison against the flat-Taylor prototype already in the tree
# -----------------------------------------------------------------------------

def compare_double_angle_vs_flat(y_slots=None, K=1, logT=3, base_degree=7,
                                 flat_degree=7, do_rescale=True):
    """Side-by-side accuracy of double-angle vs eval_mod_sine_openfhe_prototype.

    Both target post/(2*pi)*sin(2*pi*K*y).  Useful to see the accuracy gain of
    argument reduction + doubling, especially as |y| grows toward 1/(2K).
    """
    nslots = ckks_num_slots()
    if y_slots is None:
        base = [0.05, 0.20, 0.35, 0.45][:nslots]
        y_slots = base + [0.0] * (nslots - len(base))
    else:
        y_slots = list(y_slots) + [0.0] * (nslots - len(y_slots))

    post = RRK(q_list[0]) / RRK(p_enc)

    key = public_private_keygen()
    relin_key = relin_keygen(key)
    ct = encryption_ckks([CCK(v) for v in y_slots], key, p_enc, method="fft")

    ct_da = eval_mod_double_angle(ct, relin_key, post=post, K=K, logT=logT,
                                  base_degree=base_degree, do_rescale=do_rescale,
                                  check_budget=False, verbose=False)
    ct_flat = eval_mod_sine_openfhe_prototype(ct, relin_key, post=post, K=K,
                                              degree=flat_degree, do_rescale=do_rescale)

    out_da = decryption_ckks(ct_da, key, 0, 0, p_enc, nslots)
    out_flat = decryption_ckks(ct_flat, key, 0, 0, p_enc, nslots)

    print("=" * 78)
    print("double-angle vs flat-Taylor EvalMod  (POLDEG=%d, K=%d, logT=%d)" % (POLDEG, K, logT))
    print("slot :        y_j        err_double_angle      err_flat_taylor")
    for j in range(min(nslots, 8)):
        target = eval_mod_double_angle_clear(y_slots[j], post=post, K=K)
        e_da = abs(CCK(out_da[j]) - target)
        e_fl = abs(CCK(out_flat[j]) - target)
        print("%4d : %10.5f      %12.3e        %12.3e"
              % (j, float(y_slots[j]), float(e_da), float(e_fl)))
    print("=" * 78)
