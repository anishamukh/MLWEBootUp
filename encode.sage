# -----------------------------------------------------------------------------
# encode.sage
# Base-ring CKKS encoding/decoding for the ModHE / MLWE-CKKS Sage.
# -----------------------------------------------------------------------------

try:
    CKKS_PREC
except NameError:
    CKKS_PREC = 160

CCK = ComplexField(CKKS_PREC)
RRK = RealField(CKKS_PREC)

try:
    CKKS_MATRIX_REF_MAX_N
except NameError:
    CKKS_MATRIX_REF_MAX_N = 128


# Basic CKKS root/index functions

def ckks_num_slots():
    return POLDEG // 2


def _ckks_zeta():
    return CCK(exp(CCK(I) * CCK.pi() / CCK(POLDEG)))


def _ckks_omega():
    return CCK(exp(CCK(2) * CCK(I) * CCK.pi() / CCK(POLDEG)))


def _pow5_mod_2N(j):
    return ZZ(pow(ZZ(5), ZZ(j), ZZ(2*POLDEG)))


def _root_exp_to_standard_index(exp_odd):
    e = ZZ(exp_odd) % ZZ(2*POLDEG)
    if e % 2 == 0:
        raise ValueError("root exponent must be odd")
    return ZZ((e - 1) // 2)


def ckks_slot_exponents():
    return [_pow5_mod_2N(j) for j in range(ckks_num_slots())]


def ckks_slot_indices_standard_order():
    return [_root_exp_to_standard_index(e) for e in ckks_slot_exponents()]


def ckks_conjugate_indices_standard_order():
    return [_root_exp_to_standard_index(-e) for e in ckks_slot_exponents()]


def ckks_primitive_root(j):
    return _ckks_zeta() ** _pow5_mod_2N(j)


def ckks_slot_roots():
    return [ckks_primitive_root(j) for j in range(ckks_num_slots())]


def _ckks_all_crt_roots():
    roots = ckks_slot_roots()
    return roots + [z.conjugate() for z in roots]


# helpers

def _ckks_coeffs_complex(poly, n=POLDEG):
    coeffs = _dense_coeffs(poly, n)
    return [CCK(c) for c in coeffs]


def ckks_eval_poly_complex(poly, z):
    coeffs = _ckks_coeffs_complex(poly, POLDEG)
    acc = CCK(0)
    for c in reversed(coeffs):
        acc = acc*z + c
    return acc

# Lightweight recursive FFT over ComplexField

def _ckks_fft_rec(a, inverse=False):
    n = len(a)
    if n == 1:
        return [CCK(a[0])]
    if n % 2 != 0:
        raise ValueError("FFT length must be a power of two")

    even = _ckks_fft_rec(a[0::2], inverse)
    odd = _ckks_fft_rec(a[1::2], inverse)

    sign = -1 if inverse else 1
    wn = CCK(exp(CCK(sign) * CCK(2) * CCK(I) * CCK.pi() / CCK(n)))
    w = CCK(1)
    out = [CCK(0)] * n
    half = n // 2
    for k in range(half):
        t = w * odd[k]
        out[k] = even[k] + t
        out[k + half] = even[k] - t
        w *= wn
    return out


def ckks_fft_forward(a):
    return _ckks_fft_rec([CCK(x) for x in a], inverse=False)


def ckks_fft_inverse(A):
    n = len(A)
    out = _ckks_fft_rec([CCK(x) for x in A], inverse=True)
    return [x / CCK(n) for x in out]


# FFT/iFFT CKKS encode/decode

def _ckks_standard_evaluations_from_slots(values, scale=p_enc):
    nslots = ckks_num_slots()
    if len(values) > nslots:
        raise ValueError("too many CKKS slots: got %s, max is %s" % (len(values), nslots))

    vals = [CCK(v) for v in values] + [CCK(0)]*(nslots - len(values))
    evals = [CCK(0)] * POLDEG

    slot_idx = ckks_slot_indices_standard_order()
    conj_idx = ckks_conjugate_indices_standard_order()

    for j in range(nslots):
        y = CCK(scale) * vals[j]
        evals[slot_idx[j]] = y
        evals[conj_idx[j]] = y.conjugate()
    return evals


def ckks_encode_fft(values, scale=p_enc):
    evals = _ckks_standard_evaluations_from_slots(values, scale)
    b = ckks_fft_inverse(evals)
    zeta = _ckks_zeta()

    coeffs_z = []
    max_imag = RRK(0)
    pow_zeta = CCK(1)
    for i in range(POLDEG):
        c = b[i] / pow_zeta
        max_imag = max(max_imag, abs(c.imag()))
        coeffs_z.append(ZZ(round(RRK(c.real()))))
        pow_zeta *= zeta

    return R(coeffs_z) % firreducible


def ckks_decode_fft(poly, scale=p_enc, slots=None):
    if slots is None:
        slots = ckks_num_slots()
    if slots > ckks_num_slots():
        raise ValueError("requested too many slots")

    coeffs = _ckks_coeffs_complex(poly, POLDEG)
    zeta = _ckks_zeta()
    b = []
    pow_zeta = CCK(1)
    for i in range(POLDEG):
        b.append(coeffs[i] * pow_zeta)
        pow_zeta *= zeta

    evals = ckks_fft_forward(b)
    slot_idx = ckks_slot_indices_standard_order()
    return [evals[slot_idx[j]] / CCK(scale) for j in range(slots)]


# hwlpers for debugging small POLDEG only

def ckks_encode_matrix(values, scale=p_enc):
    nslots = ckks_num_slots()
    if len(values) > nslots:
        raise ValueError("too many CKKS slots: got %s, max is %s" % (len(values), nslots))

    vals = [CCK(v) for v in values] + [CCK(0)]*(nslots - len(values))
    target = [CCK(scale)*v for v in vals] + [CCK(scale)*v.conjugate() for v in vals]
    roots = _ckks_all_crt_roots()

    coeffs_z = []
    for i in range(POLDEG):
        s = CCK(0)
        for j in range(POLDEG):
            # Inverse DFT/Vandermonde over roots of X^N+1:
            # coefficient_i = (1/N) * sum_j target_j * root_j^(-i)
            s += target[j] * (roots[j] ** (-i))
        c = s / CCK(POLDEG)
        coeffs_z.append(ZZ(round(RRK(c.real()))))

    return R(coeffs_z) % firreducible


def ckks_decode_matrix(poly, scale=p_enc, slots=None):
    if slots is None:
        slots = ckks_num_slots()
    roots = ckks_slot_roots()
    return [ckks_eval_poly_complex(poly, roots[j]) / CCK(scale) for j in range(slots)]

# Public wrappers

def ckks_encode(values, scale=p_enc, method="fft"):
    if method == "fft":
        return ckks_encode_fft(values, scale)
    elif method == "matrix":
        return ckks_encode_matrix(values, scale)
    else:
        raise ValueError("unknown CKKS encode method: %s" % method)


def ckks_decode(poly, scale=p_enc, slots=None, method="fft"):
    if method == "fft":
        return ckks_decode_fft(poly, scale, slots)
    elif method == "matrix":
        return ckks_decode_matrix(poly, scale, slots)
    else:
        raise ValueError("unknown CKKS decode method: %s" % method)


def ckks_encode_real(values, scale=p_enc, method="fft"):
    return ckks_encode([CCK(v) for v in values], scale, method)


# Center lifting and q0 decoding helpers

def _center_lift_int_mod_q(a, q):
    a = ZZ(a) % ZZ(q)
    if a > ZZ(q)//2:
        a -= ZZ(q)
    return ZZ(a)


def center_lift_poly_mod_q(poly_q, q=q_list[0]):
    coeffs = _dense_coeffs(poly_q, POLDEG)
    return R([_center_lift_int_mod_q(c, q) for c in coeffs]) % firreducible


def ckks_decode_from_q0(poly_q0, scale=p_enc, slots=None, method="fft"):
    return ckks_decode(center_lift_poly_mod_q(poly_q0, q_list[0]), scale, slots, method)


# Some self-tests

def ckks_codec_selftest(values=None, scale=p_enc, method="fft"):
    """Encode/decode self-test independent of encryption."""
    if values is None:
        values = [CCK(0.25), CCK(-0.5), CCK(0.125) + CCK(0.2)*CCK(I)]
    pt = ckks_encode(values, scale, method)
    dec = ckks_decode(pt, scale, len(values), method)
    print("CKKS codec self-test, method=", method)
    for i in range(len(values)):
        print(i, "in=", CCK(values[i]), "out=", dec[i], "err=", abs(dec[i]-CCK(values[i])))
    return pt, dec


def ckks_codec_compare_fft_matrix(values=None, scale=p_enc):
    if POLDEG > CKKS_MATRIX_REF_MAX_N:
        print("Skipping matrix reference: POLDEG=", POLDEG, ">", CKKS_MATRIX_REF_MAX_N)
        return None
    if values is None:
        values = [CCK(0.25), CCK(-0.5), CCK(0.125) + CCK(0.0625)*CCK(I)]

    pt_fft = ckks_encode(values, scale, method="fft")
    pt_mat = ckks_encode(values, scale, method="matrix")
    coeff_fft = _dense_coeffs(pt_fft, POLDEG)
    coeff_mat = _dense_coeffs(pt_mat, POLDEG)
    max_coeff_diff = max(abs(ZZ(coeff_fft[i]) - ZZ(coeff_mat[i])) for i in range(POLDEG))

    dec_fft = ckks_decode(pt_fft, scale, len(values), method="fft")
    dec_mat = ckks_decode(pt_mat, scale, len(values), method="matrix")
    max_dec_diff = max(abs(dec_fft[i] - dec_mat[i]) for i in range(len(values)))

    print("CKKS FFT vs matrix comparison")
    print("max coefficient difference =", max_coeff_diff)
    print("max decoded-slot difference =", max_dec_diff)
    return pt_fft, pt_mat, dec_fft, dec_mat
