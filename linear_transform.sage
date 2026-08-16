# -----------------------------------------------------------------------------
# LinearTransforms.sage
#
# CKKS DFT/iDFT diagonal generation and homomorphic linear-transform wrappers
# for the current ModHE / MLWE-CKKS Sage prototype.
#
# Load order:
#   load("constants.sage")
#   load("auxiliary_functions.sage")
#   load("KeyGen.sage")
#   load("Encrypt-Decrypt.sage")
#   load("encode.sage")
#   load("encode-encrypt.sage")
#   load("HMult.sage")
#   load("Relin.sage")
#   load("BootstrapKeys.sage")
#   load("BootstrapOps.sage")
#   load("LinearTransforms.sage")
#
# Scope:
#   - Message remains in the BASE RING R = Z[X]/(X^POLDEG+1).
#   - No packing across the MLWE module rank is used here.
#   - CoeffToSlot and SlotToCoeff are implemented as ordinary CKKS slot linear
#     transforms, evaluated by the diagonal method using rotations and plaintext
#     multiplications.
#   - Since these transforms use Re/Im packing, they are real-linear maps:
#
#         T(z) = A*z + B*conj(z).
#
#     Therefore we generate diagonals for both A and B.  The B part is evaluated
#     by first applying homomorphic conjugation and then the same diagonal method.
# -----------------------------------------------------------------------------

try:
    CKKS_LT_TOL
except NameError:
    CKKS_LT_TOL = RealField(CKKS_PREC)(2)**(-80)


# -----------------------------------------------------------------------------
# Unrounded CKKS embedding helpers
# -----------------------------------------------------------------------------

def ckks_coeffs_to_slots_unrounded(coeffs, scale=1, slots=None):
    """Return Dcd(coeffs) without coefficient rounding.

    coeffs is a length-POLDEG list of real/complex coefficients representing
    m(X).  The output is

        (m(zeta^(5^j)) / scale)_j.

    This is the same root ordering as ckks_decode_fft(), but it accepts complex
    coefficient lists directly.  It is used only to build the DFT/iDFT matrices.
    """
    if slots is None:
        slots = ckks_num_slots()
    if len(coeffs) != POLDEG:
        raise ValueError("expected exactly POLDEG coefficients")
    if slots > ckks_num_slots():
        raise ValueError("requested too many slots")

    zeta = _ckks_zeta()
    pow_zeta = CCK(1)
    twisted = []
    for i in range(POLDEG):
        twisted.append(CCK(coeffs[i]) * pow_zeta)
        pow_zeta *= zeta

    evals = ckks_fft_forward(twisted)
    slot_idx = ckks_slot_indices_standard_order()
    return [evals[slot_idx[j]] / CCK(scale) for j in range(slots)]


def ckks_slots_to_coeffs_unrounded(values, scale=1):
    """Return inverse CKKS embedding coefficients without rounding.

    Given N/2 CKKS slots, this fills the conjugate slots, applies the inverse
    negacyclic FFT, and returns a length-POLDEG complex coefficient vector.
    """
    evals = _ckks_standard_evaluations_from_slots(values, scale)
    twisted = ckks_fft_inverse(evals)
    zeta = _ckks_zeta()

    coeffs = []
    pow_zeta = CCK(1)
    for i in range(POLDEG):
        coeffs.append(twisted[i] / pow_zeta)
        pow_zeta *= zeta
    return coeffs


# -----------------------------------------------------------------------------
# Cleartext definitions of SlotToCoeff and CoeffToSlot
# -----------------------------------------------------------------------------

def slot_to_coeff_clear_slots(slots, factor=1):
    """Cleartext SlotToCoeff map used to generate diagonals.

    Input slots z in C^(N/2).  Output is the native CKKS slot representation of
    the coefficient-packed polynomial

        m(X) = sum_j Re(factor*z_j) X^j
             + sum_j Im(factor*z_j) X^(j+N/2).

    This corresponds to the S2C step in standard CKKS bootstrapping: slots are
    converted into a polynomial whose coefficients contain the real/imaginary
    parts of the slots.
    """
    n = ckks_num_slots()
    if len(slots) != n:
        raise ValueError("slot_to_coeff_clear_slots expects exactly N/2 slots")

    f = CCK(factor)
    coeffs = [CCK(0) for _ in range(POLDEG)]
    for j in range(n):
        v = f * CCK(slots[j])
        coeffs[j] = CCK(v.real())
        coeffs[j+n] = CCK(v.imag())

    return ckks_coeffs_to_slots_unrounded(coeffs, scale=1, slots=n)


def coeff_to_slot_clear_slots(slots, factor=1):
    """Cleartext CoeffToSlot map used to generate diagonals.

    Input is the native CKKS slot representation of a coefficient-packed
    polynomial.  The polynomial coefficients are recovered by the inverse
    embedding, and the output slot j is

        factor * (coeff[j] + i*coeff[j+N/2]).

    This is the inverse of slot_to_coeff_clear_slots() when the same/reciprocal
    scaling convention is used.
    """
    n = ckks_num_slots()
    if len(slots) != n:
        raise ValueError("coeff_to_slot_clear_slots expects exactly N/2 slots")

    f = CCK(factor)
    coeffs = ckks_slots_to_coeffs_unrounded(slots, scale=1)
    out = []
    for j in range(n):
        out.append(f * (coeffs[j] + CCK(I)*coeffs[j+n]))
    return out


# -----------------------------------------------------------------------------
# Real-linear matrix decomposition: T(z) = A*z + B*conj(z)
# -----------------------------------------------------------------------------

def real_linear_matrices_from_transform(transform_fn):
    """Recover A,B matrices for a real-linear transform T.

    For each basis vector e_j:
        T(e_j)   = A[:,j] + B[:,j]
        T(i e_j) = i*A[:,j] - i*B[:,j]

    Hence:
        A[:,j] = (T(e_j) - i*T(i e_j))/2
        B[:,j] = (T(e_j) + i*T(i e_j))/2.
    """
    n = ckks_num_slots()
    A = matrix(CCK, n, n)
    B = matrix(CCK, n, n)

    for j in range(n):
        e = [CCK(0) for _ in range(n)]
        ie = [CCK(0) for _ in range(n)]
        e[j] = CCK(1)
        ie[j] = CCK(I)

        Te = [CCK(v) for v in transform_fn(e)]
        Tie = [CCK(v) for v in transform_fn(ie)]

        for i in range(n):
            A[i, j] = (Te[i] - CCK(I)*Tie[i]) / CCK(2)
            B[i, j] = (Te[i] + CCK(I)*Tie[i]) / CCK(2)

    return A, B


def _shifted_diagonal(M, r):
    """Return shifted diagonal u_r with u_r[i] = M[i,(i+r) mod n]."""
    n = M.nrows()
    return [CCK(M[i, (i+r) % n]) for i in range(n)]


def _is_zero_complex_vector(v, tol=CKKS_LT_TOL):
    return all(abs(CCK(x)) <= tol for x in v)


def matrix_to_ckks_diagonal_plaintexts(M, diag_scale=p_enc, method="fft", skip_zero=True):
    """Convert a slot matrix into CKKS plaintext diagonals.

    Returns a dictionary keyed by automorphism exponent k = 5^r mod 2N:

        { k : ckks_encode( shifted_diagonal_r(M), diag_scale ) }.

    The corresponding homomorphic evaluation is

        sum_r EncodeDiag(u_r) * Rot(ct, r).

    The encoded diagonal contains the additional CKKS plaintext scale diag_scale.
    The wrapper linear_transform_with_conjugate() can divide by diag_scale after
    summation for prototype correctness.
    """
    n = ckks_num_slots()
    if M.nrows() != n or M.ncols() != n:
        raise ValueError("matrix must be (N/2) x (N/2)")

    diags = {}
    for r in range(n):
        u = _shifted_diagonal(M, r)
        if skip_zero and _is_zero_complex_vector(u):
            continue
        k = Integer(power_mod(5, Integer(r), 2*POLDEG))
        diags[k] = ckks_encode(u, diag_scale, method=method)
    return diags


def _all_rotation_steps_from_diags(diag_bundle):
    steps = set()
    n = ckks_num_slots()
    for label in ["A", "B"]:
        if label in diag_bundle:
            for k in diag_bundle[label].keys():
                # recover r such that 5^r = k mod 2N by brute force; small helper
                kk = Integer(k) % (2*POLDEG)
                cur = Integer(1)
                for r in range(n):
                    if cur == kk:
                        steps.add(r)
                        break
                    cur = (cur * 5) % (2*POLDEG)
    return sorted(list(steps))


# -----------------------------------------------------------------------------
# Public diagonal generators
# -----------------------------------------------------------------------------

def generate_slot_to_coeff_diagonals(diag_scale=p_enc, factor=1, method="fft", skip_zero=True):
    """Generate diagonals for SlotToCoeff/S2C.

    The returned bundle has A and B diagonals for the real-linear map

        z -> native_slots( coeff_pack(z) ).

    In CKKS bootstrapping, factor is often q0/gamma or a nearby scaling factor.
    For testing the inverse relation with CoeffToSlot, keep factor=1.
    """
    T = lambda z: slot_to_coeff_clear_slots(z, factor=factor)
    A, B = real_linear_matrices_from_transform(T)
    bundle = {
        "type": "slot_to_coeff",
        "diag_scale": diag_scale,
        "factor": factor,
        "A": matrix_to_ckks_diagonal_plaintexts(A, diag_scale, method, skip_zero),
        "B": matrix_to_ckks_diagonal_plaintexts(B, diag_scale, method, skip_zero),
    }
    bundle["rotation_steps"] = _all_rotation_steps_from_diags(bundle)
    bundle["rotation_exponents"] = sorted(list(set(list(bundle["A"].keys()) + list(bundle["B"].keys()))))
    return bundle


def generate_coeff_to_slot_diagonals(diag_scale=p_enc, factor=1, method="fft", skip_zero=True):
    """Generate diagonals for CoeffToSlot/C2S.

    The returned bundle has A and B diagonals for the inverse real-linear map

        native_slots(coeff_packed_poly) -> coeff[j] + i*coeff[j+N/2].

    For a bootstrap parameterization with SlotToCoeff factor=q0/gamma, the
    corresponding CoeffToSlot factor would usually be gamma/q0, modulo the
    precise EvalMod scaling convention.  For unit tests, keep factor=1.
    """
    T = lambda z: coeff_to_slot_clear_slots(z, factor=factor)
    A, B = real_linear_matrices_from_transform(T)
    bundle = {
        "type": "coeff_to_slot",
        "diag_scale": diag_scale,
        "factor": factor,
        "A": matrix_to_ckks_diagonal_plaintexts(A, diag_scale, method, skip_zero),
        "B": matrix_to_ckks_diagonal_plaintexts(B, diag_scale, method, skip_zero),
    }
    bundle["rotation_steps"] = _all_rotation_steps_from_diags(bundle)
    bundle["rotation_exponents"] = sorted(list(set(list(bundle["A"].keys()) + list(bundle["B"].keys()))))
    return bundle


# -----------------------------------------------------------------------------
# Homomorphic evaluation wrappers
# -----------------------------------------------------------------------------

def linear_transform_with_conjugate(ct, diag_bundle, rotation_keys=None, conjugation_key=None, divide_by_diag_scale=True):
    """Evaluate a real-linear CKKS transform T(z)=A*z+B*conj(z).

    diag_bundle must be generated by generate_*_diagonals().

    This uses the current prototype's plaintext_multiply(), rotate(), conjugate(),
    add_ct(), and scalar_multiply() from BootstrapOps.sage.  Because the encoded
    diagonal plaintexts carry an extra scale diag_scale, this function optionally
    multiplies the result by 1/diag_scale at the end.  This is a prototype-scale
    correction, not a replacement for production CKKS rescaling/level management.
    """
    if rotation_keys is None:
        rotation_keys = {}

    result = None

    if len(diag_bundle.get("A", {})) > 0:
        partA = linear_transform(ct, diag_bundle["A"], rotation_keys, encode_plain=False)
        result = partA if result is None else add_ct(result, partA)

    if len(diag_bundle.get("B", {})) > 0:
        if conjugation_key is None:
            raise ValueError("conjugation_key is required for the B*conj(z) part")
        ct_conj = conjugate(ct, conjugation_key)
        partB = linear_transform(ct_conj, diag_bundle["B"], rotation_keys, encode_plain=False)
        result = partB if result is None else add_ct(result, partB)

    if result is None:
        result = zero_ct(ct_active(ct), ct_rank(ct))

    if divide_by_diag_scale:
        result = scalar_multiply(result, QQ(1)/QQ(diag_bundle["diag_scale"]), encode=False)
    return result


def slot_to_coeff_dft(ct, s2c_diags, boot_keys, divide_by_diag_scale=True):
    """Homomorphic SlotToCoeff/S2C using generated DFT diagonals."""
    return linear_transform_with_conjugate(
        ct,
        s2c_diags,
        rotation_keys=boot_keys.get("rot_keys", {}),
        conjugation_key=boot_keys.get("conj_key", None),
        divide_by_diag_scale=divide_by_diag_scale,
    )


def coeff_to_slot_dft(ct, c2s_diags, boot_keys, divide_by_diag_scale=True):
    """Homomorphic CoeffToSlot/C2S using generated inverse-DFT diagonals."""
    return linear_transform_with_conjugate(
        ct,
        c2s_diags,
        rotation_keys=boot_keys.get("rot_keys", {}),
        conjugation_key=boot_keys.get("conj_key", None),
        divide_by_diag_scale=divide_by_diag_scale,
    )

# -----------------------------------------------------------------------------
# True two-ciphertext bootstrap CoeffToSlot / SlotToCoeff
# -----------------------------------------------------------------------------
#
# These functions implement the CKKS bootstrapping version of C2S/S2C.
#
# Existing paired transform in this file:
#     one ciphertext  <->  slots coeff[j] + i*coeff[j+N/2]
#
# Bootstrap transform needed here:
#     CoeffToSlot:
#         one ciphertext decrypting to t(X)
#             -> ct_lo encrypts (t0, ..., t_{N/2-1})
#             -> ct_hi encrypts (t_{N/2}, ..., t_{N-1})
#
#     SlotToCoeff:
#         ct_lo encrypting (m0, ..., m_{N/2-1})
#         ct_hi encrypting (m_{N/2}, ..., m_{N-1})
#             -> one ciphertext decrypting to m(X)
#
# Message remains in the base ring.  No MLWE-rank packing is introduced here.
# -----------------------------------------------------------------------------


def _diag_bundle_from_clear_transform(transform_fn,
                                      bundle_type,
                                      diag_scale=p_enc,
                                      factor=1,
                                      method="fft",
                                      skip_zero=True):
    """Generate A/B CKKS diagonal plaintexts for a real-linear slot transform.

    transform_fn maps C^(N/2) -> C^(N/2).  We recover matrices A,B such that

        T(z) = A*z + B*conj(z),

    then encode the shifted diagonals as plaintext polynomials.
    """
    A, B = real_linear_matrices_from_transform(transform_fn)
    bundle = {
        "type": bundle_type,
        "diag_scale": diag_scale,
        "factor": factor,
        "A": matrix_to_ckks_diagonal_plaintexts(A, diag_scale, method, skip_zero),
        "B": matrix_to_ckks_diagonal_plaintexts(B, diag_scale, method, skip_zero),
    }
    bundle["rotation_steps"] = _all_rotation_steps_from_diags(bundle)
    bundle["rotation_exponents"] = sorted(list(set(list(bundle["A"].keys()) +
                                                   list(bundle["B"].keys()))))
    return bundle


# -----------------------------------------------------------------------------
# Cleartext bootstrap CoeffToSlot maps
# -----------------------------------------------------------------------------

def coeff_to_slot_clear_half(slots, half=0, factor=1):
    """Cleartext bootstrap CoeffToSlot half-map.

    Input:
        slots = native CKKS slots of a polynomial t(X).

    Internally:
        recover coefficient vector (t0, ..., t_{N-1}) by inverse embedding.

    Output:
        if half == 0:
            factor * (t0, ..., t_{N/2-1})
        if half == 1:
            factor * (t_{N/2}, ..., t_{N-1})

    This is the true bootstrap C2S direction: coefficients become slots.
    """
    n = ckks_num_slots()
    if len(slots) != n:
        raise ValueError("coeff_to_slot_clear_half expects exactly N/2 slots")
    if half not in [0, 1]:
        raise ValueError("half must be 0 or 1")

    f = CCK(factor)
    coeffs = ckks_slots_to_coeffs_unrounded(slots, scale=1)
    off = half * n

    return [f * CCK(coeffs[off + j]) for j in range(n)]


def coeff_to_slot_clear_lo(slots, factor=1):
    return coeff_to_slot_clear_half(slots, half=0, factor=factor)


def coeff_to_slot_clear_hi(slots, factor=1):
    return coeff_to_slot_clear_half(slots, half=1, factor=factor)


# -----------------------------------------------------------------------------
# Cleartext bootstrap SlotToCoeff maps
# -----------------------------------------------------------------------------

def slot_to_coeff_clear_from_half(slots, half=0, factor=1):
    """Cleartext bootstrap SlotToCoeff half-map.

    Input:
        slots = coefficient values packed in CKKS slots.

    Output:
        native CKKS slots of the polynomial whose selected half of coefficients
        is filled with the input slot values.

    If half == 0:
        output corresponds to polynomial sum_j factor*slots[j] X^j.

    If half == 1:
        output corresponds to polynomial sum_j factor*slots[j] X^(j+N/2).

    The full SlotToCoeff is the sum of the two half maps.
    """
    n = ckks_num_slots()
    if len(slots) != n:
        raise ValueError("slot_to_coeff_clear_from_half expects exactly N/2 slots")
    if half not in [0, 1]:
        raise ValueError("half must be 0 or 1")

    f = CCK(factor)
    off = half * n
    coeffs = [CCK(0) for _ in range(POLDEG)]

    for j in range(n):
        coeffs[off + j] = f * CCK(slots[j])

    return ckks_coeffs_to_slots_unrounded(coeffs, scale=1, slots=n)


def slot_to_coeff_clear_from_lo(slots, factor=1):
    return slot_to_coeff_clear_from_half(slots, half=0, factor=factor)


def slot_to_coeff_clear_from_hi(slots, factor=1):
    return slot_to_coeff_clear_from_half(slots, half=1, factor=factor)


# -----------------------------------------------------------------------------
# Public diagonal generators for true bootstrap C2S/S2C
# -----------------------------------------------------------------------------

def generate_coeff_to_slot_bootstrap_diagonals(diag_scale=p_enc,
                                               factor=1,
                                               method="fft",
                                               skip_zero=True):
    """Generate true bootstrap CoeffToSlot diagonals.

    Returns two diagonal bundles:

        bundle["lo"] : t(X) -> slots (t0, ..., t_{N/2-1})
        bundle["hi"] : t(X) -> slots (t_{N/2}, ..., t_{N-1})

    Each half is represented as a real-linear transform

        T(z) = A*z + B*conj(z).

    In a full bootstrap, this should be applied immediately after ModRaise.
    """
    Tlo = lambda z: coeff_to_slot_clear_half(z, half=0, factor=factor)
    Thi = lambda z: coeff_to_slot_clear_half(z, half=1, factor=factor)

    out = {
        "type": "bootstrap_coeff_to_slot_two_ct",
        "diag_scale": diag_scale,
        "factor": factor,
        "lo": _diag_bundle_from_clear_transform(Tlo,
                                                "bootstrap_coeff_to_slot_lo",
                                                diag_scale,
                                                factor,
                                                method,
                                                skip_zero),
        "hi": _diag_bundle_from_clear_transform(Thi,
                                                "bootstrap_coeff_to_slot_hi",
                                                diag_scale,
                                                factor,
                                                method,
                                                skip_zero),
    }
    out["rotation_exponents"] = merge_rotation_exponents_deep(out)
    return out


def generate_slot_to_coeff_bootstrap_diagonals(diag_scale=p_enc,
                                               factor=1,
                                               method="fft",
                                               skip_zero=True):
    """Generate true bootstrap SlotToCoeff diagonals.

    Returns two diagonal bundles:

        bundle["lo"] : slots (m0, ..., m_{N/2-1}) -> native slots of low-half polynomial
        bundle["hi"] : slots (m_{N/2}, ..., m_{N-1}) -> native slots of high-half polynomial

    The full SlotToCoeff output is:

        Apply(lo bundle to ct_lo) + Apply(hi bundle to ct_hi).
    """
    Tlo = lambda z: slot_to_coeff_clear_from_half(z, half=0, factor=factor)
    Thi = lambda z: slot_to_coeff_clear_from_half(z, half=1, factor=factor)

    out = {
        "type": "bootstrap_slot_to_coeff_two_ct",
        "diag_scale": diag_scale,
        "factor": factor,
        "lo": _diag_bundle_from_clear_transform(Tlo,
                                                "bootstrap_slot_to_coeff_lo",
                                                diag_scale,
                                                factor,
                                                method,
                                                skip_zero),
        "hi": _diag_bundle_from_clear_transform(Thi,
                                                "bootstrap_slot_to_coeff_hi",
                                                diag_scale,
                                                factor,
                                                method,
                                                skip_zero),
    }
    out["rotation_exponents"] = merge_rotation_exponents_deep(out)
    return out


# -----------------------------------------------------------------------------
# Recursive rotation-exponent collection
# -----------------------------------------------------------------------------

def merge_rotation_exponents_deep(*objects):
    """Collect all non-identity automorphism exponents from nested diagonal bundles.

    Works for:
        ordinary bundle with A/B,
        true C2S/S2C bundle with lo/hi,
        lists/tuples of bundles.
    """
    exps = set()

    def visit(obj):
        if obj is None:
            return

        if isinstance(obj, (list, tuple)):
            for x in obj:
                visit(x)
            return

        if not isinstance(obj, dict):
            return

        # Ordinary A/B diagonal bundle.
        if "A" in obj or "B" in obj:
            for label in ["A", "B"]:
                if label in obj:
                    for k in obj[label].keys():
                        kk = Integer(k) % (2*POLDEG)
                        if kk != 1:
                            exps.add(kk)

        # Existing explicit rotation_exponents field.
        if "rotation_exponents" in obj:
            for k in obj["rotation_exponents"]:
                kk = Integer(k) % (2*POLDEG)
                if kk != 1:
                    exps.add(kk)

        # Nested two-branch bundles.
        for label in ["lo", "hi", "low", "high"]:
            if label in obj:
                visit(obj[label])

    for obj in objects:
        visit(obj)

    return sorted(list(exps))


# -----------------------------------------------------------------------------
# Homomorphic true bootstrap C2S/S2C wrappers
# -----------------------------------------------------------------------------

def coeff_to_slot_bootstrap_two_ct(ct,
                                   c2s_boot_diags,
                                   boot_keys,
                                   divide_by_diag_scale=True):
    """True bootstrap CoeffToSlot.

    Input:
        ct decrypts to polynomial t(X).

    Output:
        ct_lo encrypts slots (t0, ..., t_{N/2-1})
        ct_hi encrypts slots (t_{N/2}, ..., t_{N-1})
    """
    ct_lo = linear_transform_with_conjugate(
        ct,
        c2s_boot_diags["lo"],
        rotation_keys=boot_keys.get("rot_keys", {}),
        conjugation_key=boot_keys.get("conj_key", None),
        divide_by_diag_scale=divide_by_diag_scale,
    )

    ct_hi = linear_transform_with_conjugate(
        ct,
        c2s_boot_diags["hi"],
        rotation_keys=boot_keys.get("rot_keys", {}),
        conjugation_key=boot_keys.get("conj_key", None),
        divide_by_diag_scale=divide_by_diag_scale,
    )

    return ct_lo, ct_hi


def slot_to_coeff_bootstrap_two_ct(ct_lo,
                                   ct_hi,
                                   s2c_boot_diags,
                                   boot_keys,
                                   divide_by_diag_scale=True):
    """True bootstrap SlotToCoeff.

    Input:
        ct_lo encrypts slots (m0, ..., m_{N/2-1})
        ct_hi encrypts slots (m_{N/2}, ..., m_{N-1})

    Output:
        one ciphertext decrypting to

            m(X) = sum_j m_j X^j.
    """
    part_lo = linear_transform_with_conjugate(
        ct_lo,
        s2c_boot_diags["lo"],
        rotation_keys=boot_keys.get("rot_keys", {}),
        conjugation_key=boot_keys.get("conj_key", None),
        divide_by_diag_scale=divide_by_diag_scale,
    )

    part_hi = linear_transform_with_conjugate(
        ct_hi,
        s2c_boot_diags["hi"],
        rotation_keys=boot_keys.get("rot_keys", {}),
        conjugation_key=boot_keys.get("conj_key", None),
        divide_by_diag_scale=divide_by_diag_scale,
    )

    return add_ct(part_lo, part_hi)


# -----------------------------------------------------------------------------
# Cleartext sanity test for the true two-ciphertext maps
# -----------------------------------------------------------------------------

def ckks_bootstrap_c2s_s2c_clear_selftest(values=None,
                                          factor_c2s=1,
                                          factor_s2c=1):
    """Cleartext test for true bootstrap C2S followed by S2C.

    This tests the mathematical basis-change maps before encryption.

    Starting from slots of a native CKKS polynomial:
        z = tau(m)

    C2S gives:
        lo = (m0, ..., m_{N/2-1})
        hi = (m_{N/2}, ..., m_{N-1})

    S2C merges:
        lo,hi -> tau(m) again.
    """
    n = ckks_num_slots()

    if values is None:
        coeffs = [CCK(0) for _ in range(POLDEG)]
        coeffs[0] = CCK(0.25)
        if POLDEG > 1:
            coeffs[1] = CCK(-0.5)
        if POLDEG > 2:
            coeffs[2] = CCK(0.125)
        if POLDEG > n:
            coeffs[n] = CCK(0.0625)
        native_slots = ckks_coeffs_to_slots_unrounded(coeffs, scale=1, slots=n)
    else:
        native_slots = [CCK(v) for v in values] + [CCK(0)]*(n-len(values))

    lo = coeff_to_slot_clear_half(native_slots, half=0, factor=factor_c2s)
    hi = coeff_to_slot_clear_half(native_slots, half=1, factor=factor_c2s)

    back_lo = slot_to_coeff_clear_from_half(lo, half=0, factor=factor_s2c)
    back_hi = slot_to_coeff_clear_from_half(hi, half=1, factor=factor_s2c)
    back = [back_lo[i] + back_hi[i] for i in range(n)]

    max_err = max(abs(back[i] - native_slots[i]) for i in range(n))

    print("[clear true bootstrap CoeffToSlot -> SlotToCoeff self-test]")
    print("max_err=", max_err)
    for i in range(min(8, n)):
        print(i, "in=", native_slots[i], "out=", back[i], "err=", abs(back[i]-native_slots[i]))

    return lo, hi, back, max_err

# -----------------------------------------------------------------------------
# old merge rotation (for one ciphertext StC <-> CtS)
# -----------------------------------------------------------------------------

def merge_rotation_exponents(*diag_bundles):
    """Collect all automorphism exponents needed by the given diagonal bundles."""
    exps = set()
    for b in diag_bundles:
        for k in b.get("rotation_exponents", []):
            kk = Integer(k) % (2*POLDEG)
            if kk != 1:
                exps.add(kk)
    return sorted(list(exps))


# -----------------------------------------------------------------------------
# Old cleartext self-test (for one ciphertext StC <-> CtS)
# -----------------------------------------------------------------------------

def ckks_dft_diagonal_clear_selftest(values=None, factor_s2c=1, factor_c2s=1):
    """Check the cleartext SlotToCoeff -> CoeffToSlot maps before encryption."""
    n = ckks_num_slots()
    if values is None:
        values = [CCK(0) for _ in range(n)]
        values[0] = CCK(0.25)
        if n > 1:
            values[1] = CCK(-0.5) + CCK(0.125)*CCK(I)
        if n > 2:
            values[2] = CCK(0.0625) - CCK(0.2)*CCK(I)
    else:
        values = [CCK(v) for v in values] + [CCK(0)]*(n-len(values))

    y = slot_to_coeff_clear_slots(values, factor=factor_s2c)
    z = coeff_to_slot_clear_slots(y, factor=factor_c2s)

    print("[clear SlotToCoeff -> CoeffToSlot self-test]")
    upto = min(8, n)
    max_err = max(abs(z[i] - values[i]) for i in range(n))
    for i in range(upto):
        print(i, "in=", values[i], "out=", z[i], "err=", abs(z[i]-values[i]))
    print("max_err=", max_err)
    return y, z, max_err
