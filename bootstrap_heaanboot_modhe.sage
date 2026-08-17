# -----------------------------------------------------------------------------
# bootstrap_heaanboot_modhe.sage
#
# Sage translation/adaptation of the bootstrapping routines in
#   kimandrik/HEAANBOOT, HEAAN/src/Scheme.cpp and Context.cpp
# for the ModHE / MLWE-CKKS Sage ciphertext representation used in this project.
#
# -----------------------------------------------------------------------------

import copy as _hb_copy_module

try:
    HB_PREC
except NameError:
    HB_PREC = 180

HB_CC = ComplexField(HB_PREC)
HB_RR = RealField(HB_PREC)


# Metadata wrapper (following heaanboot Ciphertext::logp/logq/slots/isComplex)

class ModHEBootCipher(object):
    def __init__(self, raw, scale=p_enc, slots=None, is_complex=True,
                 virtual_logq=None):
        self.raw = raw
        self.scale = HB_RR(scale)
        if slots is None:
            slots = POLDEG // 2
        self.slots = int(slots)
        self.is_complex = bool(is_complex)
        self.virtual_logq = virtual_logq

    def copy(self):
        return ModHEBootCipher(_hb_copy_module.deepcopy(self.raw), self.scale,
                              self.slots, self.is_complex, self.virtual_logq)

    def level(self):
        return len(self.raw[0]) - 1

    def active(self):
        return len(self.raw[0])

    def rank(self):
        return len(self.raw[1])

    def __repr__(self):
        return ("ModHEBootCipher(level=%s, rank=%s, slots=%s, scale~2^%.3f)" %
                (self.level(), self.rank(), self.slots,
                 float(log(self.scale)/log(HB_RR(2))) if self.scale > 0 else -1.0))


def hb_wrap(raw, scale=p_enc, slots=None, is_complex=True):
    return ModHEBootCipher(raw, scale, slots, is_complex)


def hb_unwrap(cipher):
    if isinstance(cipher, ModHEBootCipher):
        return cipher.raw
    return cipher


# Basic ciphertext/RNS helpers

def _hb_raw_copy(ct):
    return _hb_copy_module.deepcopy(ct)


def _hb_raw_level(ct):
    return len(ct[0]) - 1


def _hb_raw_rank(ct):
    return len(ct[1])


def _hb_raw_mod_down_to_level(ct, target_level):
    """RNS analogue of HEAAN modDownTo: drop high q-shares, no scaling."""
    target_level = int(target_level)
    if target_level < 0:
        raise ValueError("target_level must be >= 0")
    if target_level > _hb_raw_level(ct):
        raise ValueError("cannot modulus-down to a higher level")
    active = target_level + 1
    return [list(ct[0][:active]), [list(comp[:active]) for comp in ct[1]]]


def hb_mod_down_to_level(cipher, target_level):
    out = cipher.copy()
    out.raw = _hb_raw_mod_down_to_level(out.raw, target_level)
    return out


def _hb_align_raw_levels(ct1, ct2):
    level = min(_hb_raw_level(ct1), _hb_raw_level(ct2))
    return (_hb_raw_mod_down_to_level(ct1, level),
            _hb_raw_mod_down_to_level(ct2, level), level)


def _hb_raw_add(ct1, ct2):
    a, b, level = _hb_align_raw_levels(ct1, ct2)
    if _hb_raw_rank(a) != _hb_raw_rank(b):
        raise ValueError("ciphertext ranks do not match")
    rank = _hb_raw_rank(a)
    active = level + 1
    c0 = [poly_mod_prime(a[0][j] + b[0][j], q_list[j]) for j in range(active)]
    c1 = []
    for k in range(rank):
        c1.append([poly_mod_prime(a[1][k][j] + b[1][k][j], q_list[j])
                   for j in range(active)])
    return [c0, c1]


def _hb_raw_sub(ct1, ct2):
    a, b, level = _hb_align_raw_levels(ct1, ct2)
    if _hb_raw_rank(a) != _hb_raw_rank(b):
        raise ValueError("ciphertext ranks do not match")
    rank = _hb_raw_rank(a)
    active = level + 1
    c0 = [poly_mod_prime(a[0][j] - b[0][j], q_list[j]) for j in range(active)]
    c1 = []
    for k in range(rank):
        c1.append([poly_mod_prime(a[1][k][j] - b[1][k][j], q_list[j])
                   for j in range(active)])
    return [c0, c1]

def hb_sub_into_second(c1, c2):
    """
    HEAAN Ring2Utils::subAndEqual2:
        c2 <- c1 - c2

    IMPORTANT:
    the result retains c2's metadata.
    """
    if c1.level() != c2.level():
        raise ValueError(
            "subAndEqual2 requires equal levels: %d vs %d"
            % (c1.level(), c2.level())
        )

    out = c2.copy()                  # <-- retain SECOND operand metadata
    out.raw = _hb_raw_sub(c1.raw, c2.raw)

    # Do NOT overwrite out.scale with c1.scale.
    return out

def hb_add(c1, c2):
    out = c1.copy()
    out.raw = _hb_raw_add(c1.raw, c2.raw)
    out.scale = c1.scale
    return out


def hb_sub(c1, c2):
    out = c1.copy()
    out.raw = _hb_raw_sub(c1.raw, c2.raw)
    out.scale = c1.scale
    return out


def hb_negate(c):
    out = c.copy()
    active = out.active()
    out.raw[0] = [poly_mod_prime(-out.raw[0][j], q_list[j]) for j in range(active)]
    for k in range(out.rank()):
        out.raw[1][k] = [poly_mod_prime(-out.raw[1][k][j], q_list[j])
                         for j in range(active)]
    return out


# Centered CRT reconstruction and heaanboot-style bit shifts

def _hb_crt_integer(residues, moduli_list):
    """CRT without relying on Sage's overloaded crt() signatures."""
    mods = [ZZ(m) for m in moduli_list]
    M = ZZ(1)
    for m in mods:
        M *= m
    acc = ZZ(0)
    for a, m in zip(residues, mods):
        Mi = M // m
        acc += (ZZ(a) % m) * Mi * Mi.inverse_mod(m)
    return acc % M


def _hb_center_crt(residues, moduli_list):
    M = ZZ(1)
    for m in moduli_list:
        M *= ZZ(m)
    z = _hb_crt_integer(residues, moduli_list)
    if z > M // 2:
        z -= M
    return z


def _hb_round_div_signed(z, divisor):
    """Nearest-integer division, symmetric for negative integers."""
    z = ZZ(z)
    d = ZZ(divisor)
    if d <= 0:
        raise ValueError("divisor must be positive")
    if z >= 0:
        return (z + d//2) // d
    return -((-z + d//2) // d)


def _hb_reconstruct_component(poly_shares):
    active = len(poly_shares)
    mods = [ZZ(q_list[j]) for j in range(active)]
    coeff_by_share = [_dense_coeffs(poly_shares[j], POLDEG) for j in range(active)]
    coeffs = []
    for i in range(POLDEG):
        residues = [ZZ(coeff_by_share[j][i]) for j in range(active)]
        coeffs.append(_hb_center_crt(residues, mods))
    return coeffs


def _hb_reshare_integer_coeffs(coeffs, active):
    poly = R([ZZ(c) for c in coeffs]) % firreducible
    return [poly_mod_prime(poly, q_list[j]) for j in range(active)]


def _hb_raw_div_pow2_round(ct, bits):
    """Slow Sage analogue of Ring2Utils::rightShift on every ct component."""
    bits = int(bits)
    if bits <= 0:
        return _hb_raw_copy(ct)
    divisor = ZZ(1) << bits
    active = len(ct[0])

    c0z = _hb_reconstruct_component(ct[0])
    c0z = [_hb_round_div_signed(z, divisor) for z in c0z]
    c0 = _hb_reshare_integer_coeffs(c0z, active)

    c1 = []
    for comp in ct[1]:
        z = _hb_reconstruct_component(comp)
        z = [_hb_round_div_signed(v, divisor) for v in z]
        c1.append(_hb_reshare_integer_coeffs(z, active))
    return [c0, c1]


"""def hb_div_by_po2(cipher, bits):
    #HEAAN divByPo2: divide ct coefficients, keep plaintext scale metadata.
    out = cipher.copy()
    out.raw = _hb_raw_div_pow2_round(out.raw, bits)
    if out.virtual_logq is not None:
        out.virtual_logq -= int(bits)
    return out
"""

def hb_div_by_po2(cipher, bits):

    bits = int(bits)

    if bits <= 0:
        return cipher.copy()

    if cipher.level() <= 0:
        raise ValueError("no RNS level left for hb_div_by_po2")

    qdrop = ZZ(q_list[cipher.level()])
    two_pow = ZZ(1) << bits

    # Encode 1/2^bits at scale qdrop.
    # After dropping qdrop:
    #
    #   new scale = old_scale * qdrop / qdrop = old_scale
    #   new value = old_value / 2^bits.
    tmp = hb_mult_by_const(
        cipher,
        HB_RR(1) / HB_RR(two_pow),
        plain_scale=qdrop
    )

    return hb_rescale_one(tmp)

"""
def hb_virtual_rescale_by_po2(cipher, bits):
    #HEAAN reScaleBy(bits) for a small non-RNS-aligned bit count.

    #Both coefficients and scale are divided by 2^bits, so the represented plaintext value is unchanged.  Physical RNS primes are deliberately kept.
    bits = int(bits)
    out = cipher.copy()
    if bits <= 0:
        return out
    out.raw = _hb_raw_div_pow2_round(out.raw, bits)
    out.scale /= HB_RR(ZZ(1) << bits)
    if out.virtual_logq is not None:
        out.virtual_logq -= bits
    return out
"""
def hb_virtual_rescale_by_po2(cipher, bits):

    bits = int(bits)

    if bits <= 0:
        return cipher.copy()

    if cipher.level() <= 0:
        raise ValueError(
            "no RNS level left for hb_virtual_rescale_by_po2"
        )

    qdrop = ZZ(q_list[cipher.level()])
    two_pow = ZZ(1) << bits

    pscale = ZZ(round(
        HB_RR(qdrop) / HB_RR(two_pow)
    ))

    if pscale <= 0:
        raise ValueError("invalid plaintext scale in small rescale")

    # Plaintext being multiplied is 1.
    tmp = hb_mult_by_const(
        cipher,
        HB_RR(1),
        plain_scale=pscale
    )

    # After the ordinary qdrop rescale:
    #
    # scale_new = scale_old * pscale / qdrop
    #           ~= scale_old / 2^bits,
    #
    # while plaintext value stays unchanged.
    return hb_rescale_one(tmp)


# True ModHE RNS rescaling (one q-prime)

def hb_rescale_one(cipher):
    """Consume the current top RNS q-prime using ModHE rescaling."""
    level = cipher.level()
    if level <= 0:
        raise ValueError("no RNS level left for rescaling")

    qdrop = ZZ(q_list[level])
    out = cipher.copy()
    out.raw[0] = rescaling_ciphertext(cipher.raw[0], level, level-1)
    out.raw[1] = [rescaling_ciphertext(comp, level, level-1)
                  for comp in cipher.raw[1]]
    out.scale /= HB_RR(qdrop)
    return out


def hb_mod_down_one(cipher):
    """Drop one RNS prime without changing coefficients/scale."""
    if cipher.level() <= 0:
        raise ValueError("no RNS level left")
    return hb_mod_down_to_level(cipher, cipher.level()-1)



# Plaintext polynomial multiplication / constants / i and -i

def _hb_raw_mult_by_poly(ct, poly):
    active = len(ct[0])
    rank = len(ct[1])
    pshares = [poly_mod_prime(poly, q_list[j]) for j in range(active)]
    c0 = [poly_mod_prime((ct[0][j] * pshares[j]) % firreducible, q_list[j])
          for j in range(active)]
    c1 = []
    for k in range(rank):
        c1.append([poly_mod_prime((ct[1][k][j] * pshares[j]) % firreducible,
                                  q_list[j])
                   for j in range(active)])
    return [c0, c1]


def hb_mult_by_poly(cipher, poly, plain_scale):
    out = cipher.copy()
    out.raw = _hb_raw_mult_by_poly(cipher.raw, poly)
    out.scale *= HB_RR(plain_scale)
    return out


def _hb_constant_poly(value, scale):
    """HEAAN encodeSingle-style constant: real at 0, imag at N/2."""
    z = HB_CC(value)
    coeffs = [ZZ(0)] * POLDEG
    coeffs[0] = ZZ(round(HB_RR(z.real()) * HB_RR(scale)))
    if z.imag() != 0:
        coeffs[POLDEG//2] = ZZ(round(HB_RR(z.imag()) * HB_RR(scale)))
    return R(coeffs) % firreducible


"""def hb_add_const(cipher, value):
    out = cipher.copy()
    z = HB_CC(value)
    zr = ZZ(round(HB_RR(z.real()) * out.scale))
    zi = ZZ(round(HB_RR(z.imag()) * out.scale))
    active = out.active()
    for j in range(active):
        addpoly = R([zr] + [0]*(POLDEG-1))
        if zi != 0:
            coeffs = [ZZ(0)]*POLDEG
            coeffs[0] = zr
            coeffs[POLDEG//2] = zi
            addpoly = R(coeffs)
        out.raw[0][j] = poly_mod_prime(out.raw[0][j] + addpoly, q_list[j])
    return out"""

def hb_add_const(cipher, value, const_scale=None):
    """
    HEAAN addConst analogue.

    If const_scale is supplied, encode the added constant at that
    explicit scale, but retain the ciphertext's own scale metadata.
    """
    out = cipher.copy()

    if const_scale is None:
        const_scale = out.scale

    z = HB_CC(value)

    zr = ZZ(round(
        HB_RR(z.real()) * HB_RR(const_scale)
    ))
    zi = ZZ(round(
        HB_RR(z.imag()) * HB_RR(const_scale)
    ))

    active = out.active()

    coeffs = [ZZ(0)] * POLDEG
    coeffs[0] = zr
    if zi != 0:
        coeffs[POLDEG//2] = zi

    addpoly = R(coeffs) % firreducible

    for j in range(active):
        out.raw[0][j] = poly_mod_prime(
            out.raw[0][j] + addpoly,
            q_list[j]
        )

    # IMPORTANT:
    # out.scale is unchanged.
    return out

def hb_mult_by_const(cipher, value, plain_scale=None):
    """Multiply by a scalar/plain constant without rescaling."""
    if plain_scale is None:
        # RNS adaptation: choose the top q as the plaintext scale.  A following
        # hb_rescale_one then preserves the ciphertext scale exactly.
        plain_scale = ZZ(q_list[cipher.level()])
    poly = _hb_constant_poly(value, plain_scale)
    return hb_mult_by_poly(cipher, poly, plain_scale)


"""def hb_mult_by_const_rescale(cipher, value):
    return hb_rescale_one(hb_mult_by_const(cipher, value))"""

def hb_mult_by_const_rescale(cipher, value, plain_scale=None):
    return hb_rescale_one(
        hb_mult_by_const(
            cipher,
            value,
            plain_scale=plain_scale
        )
    )


def hb_mult_by_monomial(cipher, degree):
    out = cipher.copy()
    mon = (x ** ZZ(degree)) % firreducible
    out.raw = _hb_raw_mult_by_poly(out.raw, mon)
    # Monomial has scale 1.
    return out


def hb_imult(cipher):
    """Multiply slots by i, matching Scheme::imultAndEqual."""
    return hb_mult_by_monomial(cipher, POLDEG//2)


def hb_idiv(cipher):
    """Multiply slots by -i, matching Scheme::idivAndEqual."""
    return hb_mult_by_monomial(cipher, 3*(POLDEG//2))


# Homomorphic multiplication + relinearization (added here again because modhe mult+relin code was not friendly to parameterization)


def hb_mult(c1, c2, relin_key):
    a, b, level = _hb_align_raw_levels(c1.raw, c2.raw)
    nonlin = homomorphic_mult(a, b, level)
    raw = relinearisation(nonlin, relin_key)
    out = c1.copy()
    out.raw = raw
    out.scale = c1.scale * c2.scale
    return out


def hb_mult_rescale(c1, c2, relin_key):
    return hb_rescale_one(hb_mult(c1, c2, relin_key))


def hb_square(cipher, relin_key):
    return hb_mult(cipher, cipher, relin_key)


def hb_square_rescale(cipher, relin_key):
    return hb_rescale_one(hb_square(cipher, relin_key))


# Ring automorphisms, rotations and conjugation for MLWE ciphertexts

def hb_automorphism_poly(poly, exponent):
    """Apply X -> X^exponent in Z[X]/(X^N+1), exponent odd mod 2N."""
    N = POLDEG
    M = 2*N
    e = int(ZZ(exponent) % M)
    if e % 2 == 0:
        raise ValueError("automorphism exponent must be odd modulo 2N")

    src = _dense_coeffs(poly, N)
    dst = [ZZ(0)] * N
    for i in range(N):
        a = ZZ(src[i])
        t = (i * e) % M
        if t >= N:
            dst[t-N] -= a
        else:
            dst[t] += a
    return R(dst) % firreducible


def _hb_raw_automorphism(ct, exponent):
    active = len(ct[0])
    rank = len(ct[1])
    c0 = []
    for j in range(active):
        c0.append(poly_mod_prime(hb_automorphism_poly(ct[0][j], exponent), q_list[j]))
    c1 = []
    for k in range(rank):
        row = []
        for j in range(active):
            row.append(poly_mod_prime(hb_automorphism_poly(ct[1][k][j], exponent), q_list[j]))
        c1.append(row)
    return [c0, c1]


def hb_rotation_exponent(rot_slots):
    return pow(ZZ(5), ZZ(rot_slots), ZZ(2*POLDEG))


def hb_conjugation_exponent():
    return ZZ(2*POLDEG - 1)


def _hb_secret_automorphism(secret_or_key, exponent):
    slist = _secret_to_list(secret_or_key)
    return [hb_automorphism_poly(s, exponent) for s in slist]


class HEAANBootKeys(object):
    def __init__(self, rot_keys, conj_key, rotation_indices):
        self.rot_keys = rot_keys
        self.conj_key = conj_key
        self.rotation_indices = list(rotation_indices)


def heaan_boot_required_rotations(slots, exact_cpp_keyset=False):
    """Rotation indices needed by SubSum, BSGS C2S/S2C and sparse EvalExp."""
    slots = int(slots)
    if slots < 1 or slots > POLDEG//2 or (slots & (slots-1)) != 0:
        raise ValueError("slots must be a power of two in [1,N/2]")
    logSlots = ZZ(slots).nbits() - 1
    logNh = ZZ(POLDEG//2).nbits() - 1
    logk = logSlots // 2
    k = 1 << int(logk)
    m = 1 << int(logSlots - logk)

    rots = set()
    # SubSum rotations: HEAAN leftRotateByPo2(cipher, i).
    for i in range(int(logSlots), int(logNh)):
        rots.add(1 << i)
    # BSGS baby rotations.
    for j in range(1, k):
        rots.add(j)
    # BSGS giant rotations.
    for i in range(1, m):
        rots.add(i*k)
    # Sparse EvalExp uses leftRotateFast(..., slots).
    if slots < POLDEG//2:
        rots.add(slots)

    if exact_cpp_keyset:
        # Scheme::addBootKey first calls addLeftRotKeys(), i.e. all powers of 2.
        for i in range(int(logNh)):
            rots.add(1 << i)

    return sorted([r for r in rots if r > 0 and r < POLDEG//2])


def heaan_boot_keygen(key, slots, exact_cpp_keyset=False, verbose=True):
    """MLWE generalization of HEAAN addBootKey().

    A rotation/conjugation first applies an automorphism to every ciphertext
    component, which changes the secret from s to phi(s).  We then reuse the
    project's general rankup_keygen() as a same-rank key switch phi(s) -> s.
    """
    s = key[0]
    rot_indices = heaan_boot_required_rotations(slots, exact_cpp_keyset)
    rot_keys = {}

    if verbose:
        print("Generating", len(rot_indices), "HEAANBOOT rotation keys")
    for rot in rot_indices:
        expn = hb_rotation_exponent(rot)
        srot = _hb_secret_automorphism(s, expn)
        rot_keys[int(rot)] = rankup_keygen(srot, s)
        #print("\nRotation key=", rot_keys)

    sconj = _hb_secret_automorphism(s, hb_conjugation_exponent())
    conj_key = rankup_keygen(sconj, s)
    #print("\nConjugation key=", conj_key)
    return HEAANBootKeys(rot_keys, conj_key, rot_indices)


def hb_left_rotate_fast(cipher, rot_slots, boot_keys):
    rot_slots = int(rot_slots)
    if rot_slots == 0:
        return cipher.copy()
    if rot_slots not in boot_keys.rot_keys:
        raise KeyError("missing HEAANBOOT rotation key for rotation %s" % rot_slots)
    expn = hb_rotation_exponent(rot_slots)
    raw_auto = _hb_raw_automorphism(cipher.raw, expn)
    raw_sw = rank_up(raw_auto, boot_keys.rot_keys[rot_slots])
    out = cipher.copy()
    out.raw = raw_sw
    return out


def hb_conjugate(cipher, boot_keys):
    raw_auto = _hb_raw_automorphism(cipher.raw, hb_conjugation_exponent())
    raw_sw = rank_up(raw_auto, boot_keys.conj_key)
    out = cipher.copy()
    out.raw = raw_sw
    return out



# Exact HEAAN Context::fftSpecialInv and addBootContext translation

def _hb_bit_reverse(vals):
    vals = list(vals)
    size = len(vals)
    j = 0
    for i in range(1, size):
        bitv = size >> 1
        while j >= bitv and bitv > 0:
            j -= bitv
            bitv >>= 1
        j += bitv
        if i < j:
            vals[i], vals[j] = vals[j], vals[i]
    return vals


def _hb_ksi_pow(idx):
    M = ZZ(2*POLDEG)
    idx = ZZ(idx) % M
    return HB_CC(exp(HB_CC(2)*HB_CC.pi()*HB_CC(I)*HB_CC(idx)/HB_CC(M)))


def _hb_rot_group(size=None):
    Nh = POLDEG//2
    if size is None:
        size = Nh
    vals = []
    v = ZZ(1)
    M = ZZ(2*POLDEG)
    for _ in range(int(size)):
        vals.append(v)
        v = (v*5) % M
    return vals


def _hb_fft_special_inv(vals):
    """Literal translation of Context::fftSpecialInv()."""
    vals = [HB_CC(v) for v in vals]
    size = len(vals)
    if size < 1 or (size & (size-1)) != 0:
        raise ValueError("special FFT size must be a power of two")
    rotGroup = _hb_rot_group(size)
    M = ZZ(2*POLDEG)

    length = size
    while length >= 1:
        lenh = length >> 1
        lenq = length << 2
        for base in range(0, size, length):
            for j in range(lenh):
                idx = (lenq - (int(rotGroup[j]) % lenq)) * int(M) // lenq
                u = vals[base+j] + vals[base+j+lenh]
                v = vals[base+j] - vals[base+j+lenh]
                v *= _hb_ksi_pow(idx)
                vals[base+j] = u
                vals[base+j+lenh] = v
        length >>= 1

    vals = _hb_bit_reverse(vals)
    return [v / HB_CC(size) for v in vals]


def _hb_right_rotate_values(vals, rot):
    vals = list(vals)
    n = len(vals)
    if n == 0:
        return vals
    r = int(rot) % n
    if r == 0:
        return vals
    return vals[n-r:] + vals[:n-r]


def _hb_special_poly_unit_coeffs(transformed_vals, gap):
    """C++ Context encoding layout, but return unscaled real coefficients.

    The real parts go to coefficient i*gap and imaginary parts to N/2+i*gap.
    """
    coeffs = [HB_RR(0)] * POLDEG
    Nh = POLDEG//2
    for i, z in enumerate(transformed_vals):
        idx = i*int(gap)
        jdx = Nh + i*int(gap)
        if idx >= Nh or jdx >= POLDEG:
            raise ValueError("invalid HEAAN special polynomial layout")
        coeffs[idx] = HB_RR(HB_CC(z).real())
        coeffs[jdx] = HB_RR(HB_CC(z).imag())
    return coeffs


def _hb_scale_unit_poly(unit_coeffs, scale):
    coeffs = [ZZ(round(HB_RR(c)*HB_RR(scale))) for c in unit_coeffs]
    return R(coeffs) % firreducible


class HEAANBootContext(object):
    def __init__(self, slots, pvec_unit, pvec_inv_unit,
                 p1_unit=None, p2_unit=None, nominal_logp=None):
        self.slots = int(slots)
        self.logSlots = ZZ(slots).nbits() - 1
        self.pvec_unit = pvec_unit
        self.pvecInv_unit = pvec_inv_unit
        self.p1_unit = p1_unit
        self.p2_unit = p2_unit
        self.nominal_logp = nominal_logp

    def pvec(self, idx, scale):
        return _hb_scale_unit_poly(self.pvec_unit[int(idx)], scale)

    def pvecInv(self, idx, scale):
        return _hb_scale_unit_poly(self.pvecInv_unit[int(idx)], scale)

    def p1(self, scale):
        if self.p1_unit is None:
            raise ValueError("p1 only exists for sparse packing")
        return _hb_scale_unit_poly(self.p1_unit, scale)

    def p2(self, scale):
        if self.p2_unit is None:
            raise ValueError("p2 only exists for sparse packing")
        return _hb_scale_unit_poly(self.p2_unit, scale)


def heaan_add_boot_context(slots, nominal_logp=None, verbose=True):

    slots = int(slots)
    Nh = POLDEG//2
    M = 2*POLDEG
    logNh = ZZ(Nh).nbits() - 1
    if slots < 1 or slots > Nh or (slots & (slots-1)) != 0:
        raise ValueError("slots must be a power of two in [1,N/2]")
    logSlots = ZZ(slots).nbits() - 1
    dslots = 2*slots
    logk = logSlots >> 1
    k = 1 << int(logk)
    gap = Nh >> int(logSlots)
    rotGroup = _hb_rot_group(Nh)

    pvec = [None] * slots
    pvecInv = [None] * slots
    p1 = None
    p2 = None
    c = HB_RR(1) / (HB_RR(4)*HB_RR.pi())

    if logSlots < logNh:
        dgap = gap >> 1
        for ki in range(0, slots, k):
            for pos in range(ki, ki+k):
                pvals = [HB_CC(0)] * dslots
                for i in range(0, slots-pos):
                    deg = ((M - int(rotGroup[i+pos])) * i * gap) % M
                    z = _hb_ksi_pow(deg)
                    pvals[i] = z
                    pvals[i+slots] = (HB_CC(-z.imag()) + HB_CC(I)*HB_CC(z.real()))
                for i in range(slots-pos, slots):
                    deg = ((M - int(rotGroup[i+pos-slots])) * i * gap) % M
                    z = _hb_ksi_pow(deg)
                    pvals[i] = z
                    pvals[i+slots] = (HB_CC(-z.imag()) + HB_CC(I)*HB_CC(z.real()))
                pvals = _hb_right_rotate_values(pvals, ki)
                pvals = _hb_fft_special_inv(pvals)
                pvec[pos] = _hb_special_poly_unit_coeffs(pvals, dgap)

        pvals = [HB_CC(0)] * dslots
        for i in range(slots):
            pvals[i] = HB_CC(0)
            pvals[i+slots] = (HB_CC(I)*HB_CC(-c))
        pvals = _hb_fft_special_inv(pvals)
        p1 = _hb_special_poly_unit_coeffs(pvals, dgap)

        pvals = [HB_CC(0)] * dslots
        for i in range(slots):
            pvals[i] = HB_CC(c)
            pvals[i+slots] = HB_CC(0)
        pvals = _hb_fft_special_inv(pvals)
        p2 = _hb_special_poly_unit_coeffs(pvals, dgap)

    else:
        for ki in range(0, slots, k):
            for pos in range(ki, ki+k):
                pvals = [HB_CC(0)] * slots
                for i in range(0, slots-pos):
                    deg = ((M - int(rotGroup[i+pos])) * i * gap) % M
                    pvals[i] = _hb_ksi_pow(deg)
                for i in range(slots-pos, slots):
                    deg = ((M - int(rotGroup[i+pos-slots])) * i * gap) % M
                    pvals[i] = _hb_ksi_pow(deg)
                pvals = _hb_right_rotate_values(pvals, ki)
                pvals = _hb_fft_special_inv(pvals)
                pvec[pos] = _hb_special_poly_unit_coeffs(pvals, gap)

    # pvecInv is common to sparse/full cases in Context.cpp.
    for ki in range(0, slots, k):
        for pos in range(ki, ki+k):
            pvals = [HB_CC(0)] * slots
            for i in range(0, slots-pos):
                deg = (int(rotGroup[i]) * (i+pos) * gap) % M
                pvals[i] = _hb_ksi_pow(deg)
            for i in range(slots-pos, slots):
                deg = (int(rotGroup[i]) * (i+pos-slots) * gap) % M
                pvals[i] = _hb_ksi_pow(deg)
            pvals = _hb_right_rotate_values(pvals, ki)
            pvals = _hb_fft_special_inv(pvals)
            pvecInv[pos] = _hb_special_poly_unit_coeffs(pvals, gap)

    if verbose:
        print("HEAAN boot context generated: slots=", slots,
              "k=", k, "sparse=", bool(logSlots < logNh))

    return HEAANBootContext(slots, pvec, pvecInv, p1, p2, nominal_logp)



# CoeffToSlot / SlotToCoeff (following heaanboot, notw that this is different from the later BT procedures proposed)


def _hb_mult_by_unit_poly(cipher, unit_coeffs):
    # Scale plaintext by the current top RNS prime so a following rescale keeps
    # ciphertext scale approximately/exactly stable.
    pscale = ZZ(q_list[cipher.level()])
    poly = _hb_scale_unit_poly(unit_coeffs, pscale)
    return hb_mult_by_poly(cipher, poly, pscale)


def heaan_coeff_to_slot(cipher, boot_context, boot_keys):
    slots = cipher.slots
    logSlots = ZZ(slots).nbits() - 1
    logk = logSlots // 2
    k = 1 << int(logk)

    rotvec = [None] * k
    rotvec[0] = cipher.copy()
    for j in range(1, k):
        rotvec[j] = hb_left_rotate_fast(rotvec[0], j, boot_keys)

    tmpvec = [None] * k
    for j in range(k):
        tmpvec[j] = _hb_mult_by_unit_poly(rotvec[j], boot_context.pvec_unit[j])
    for j in range(1, k):
        tmpvec[0] = hb_add(tmpvec[0], tmpvec[j])
    out = tmpvec[0]

    for ki in range(k, slots, k):
        for j in range(k):
            tmpvec[j] = _hb_mult_by_unit_poly(rotvec[j], boot_context.pvec_unit[j+ki])
        for j in range(1, k):
            tmpvec[0] = hb_add(tmpvec[0], tmpvec[j])
        tmp0 = hb_left_rotate_fast(tmpvec[0], ki, boot_keys)
        out = hb_add(out, tmp0)

    return hb_rescale_one(out)


def heaan_slot_to_coeff(cipher, boot_context, boot_keys):
    slots = cipher.slots
    logSlots = ZZ(slots).nbits() - 1
    logk = logSlots // 2
    k = 1 << int(logk)

    rotvec = [None] * k
    rotvec[0] = cipher.copy()
    for j in range(1, k):
        rotvec[j] = hb_left_rotate_fast(rotvec[0], j, boot_keys)

    tmpvec = [None] * k
    for j in range(k):
        tmpvec[j] = _hb_mult_by_unit_poly(rotvec[j], boot_context.pvecInv_unit[j])
    for j in range(1, k):
        tmpvec[0] = hb_add(tmpvec[0], tmpvec[j])
    out = tmpvec[0]

    for ki in range(k, slots, k):
        for j in range(k):
            tmpvec[j] = _hb_mult_by_unit_poly(rotvec[j], boot_context.pvecInv_unit[j+ki])
        for j in range(1, k):
            tmpvec[0] = hb_add(tmpvec[0], tmpvec[j])
        tmp0 = hb_left_rotate_fast(tmpvec[0], ki, boot_keys)
        out = hb_add(out, tmp0)

    return hb_rescale_one(out)


# exp2piAndEqual (degree-7 grouped polynomial)

"""def heaan_exp2pi(cipher, relin_key):
    #Direct RNS translation of Scheme::exp2piAndEqual().
    pi = HB_RR.pi()

    cipher2 = hb_square_rescale(cipher, relin_key)
    cipher4 = hb_square_rescale(cipher2, relin_key)

    cipher01 = hb_add_const(cipher, HB_RR(1)/(HB_RR(2)*pi))
    cipher01 = hb_mult_by_const_rescale(cipher01, HB_RR(2)*pi)

    cipher23 = hb_add_const(cipher, HB_RR(3)/(HB_RR(2)*pi))
    cipher23 = hb_mult_by_const_rescale(cipher23, HB_RR(4)*pi**3/HB_RR(3))
    cipher23 = hb_mult_rescale(cipher23, cipher2, relin_key)
    cipher23 = hb_add(cipher23, cipher01)

    cipher45 = hb_add_const(cipher, HB_RR(5)/(HB_RR(2)*pi))
    cipher45 = hb_mult_by_const_rescale(cipher45, HB_RR(4)*pi**5/HB_RR(15))

    out = hb_add_const(cipher, HB_RR(7)/(HB_RR(2)*pi))
    out = hb_mult_by_const_rescale(out, HB_RR(8)*pi**7/HB_RR(315))
    out = hb_mult_rescale(out, cipher2, relin_key)

    # C++ modDownBy(cipher45, logp) before adding at cipher's lower modulus.
    cipher45 = hb_mod_down_to_level(cipher45, out.level())
    out = hb_add(out, cipher45)

    out = hb_mult_rescale(out, cipher4, relin_key)
    cipher23 = hb_mod_down_to_level(cipher23, out.level())
    out = hb_add(out, cipher23)
    return out
"""

def heaan_exp2pi(cipher, relin_key, plain_scale):
    """
    Degree-7 grouped exponential approximation.    
    """
    pi = HB_RR.pi()

    cipher2 = hb_square_rescale(cipher, relin_key)
    cipher4 = hb_square_rescale(cipher2, relin_key)

    cipher01 = hb_add_const(
        cipher,
        HB_RR(1)/(HB_RR(2)*pi),
        const_scale=plain_scale
    )
    cipher01 = hb_mult_by_const_rescale(
        cipher01,
        HB_RR(2)*pi,
        plain_scale=plain_scale
    )

    cipher23 = hb_add_const(
        cipher,
        HB_RR(3)/(HB_RR(2)*pi),
        const_scale=plain_scale
    )
    cipher23 = hb_mult_by_const_rescale(
        cipher23,
        HB_RR(4)*pi**3/HB_RR(3),
        plain_scale=plain_scale
    )
    cipher23 = hb_mult_rescale(
        cipher23,
        cipher2,
        relin_key
    )
    cipher23 = hb_add(cipher23, cipher01)

    cipher45 = hb_add_const(
        cipher,
        HB_RR(5)/(HB_RR(2)*pi),
        const_scale=plain_scale
    )
    cipher45 = hb_mult_by_const_rescale(
        cipher45,
        HB_RR(4)*pi**5/HB_RR(15),
        plain_scale=plain_scale
    )

    out = hb_add_const(
        cipher,
        HB_RR(7)/(HB_RR(2)*pi),
        const_scale=plain_scale
    )
    out = hb_mult_by_const_rescale(
        out,
        HB_RR(8)*pi**7/HB_RR(315),
        plain_scale=plain_scale
    )
    out = hb_mult_rescale(
        out,
        cipher2,
        relin_key
    )

    cipher45 = hb_mod_down_to_level(
        cipher45,
        out.level()
    )
    out = hb_add(out, cipher45)

    out = hb_mult_rescale(
        out,
        cipher4,
        relin_key
    )

    cipher23 = hb_mod_down_to_level(
        cipher23,
        out.level()
    )
    out = hb_add(out, cipher23)

    return out

# EvalExp


def heaan_eval_exp(cipher, boot_context, boot_keys, relin_key, logT, logI=4):
    boot_plain_scale = ZZ(16) * ZZ(q_list[0])
    slots = cipher.slots
    Nh = POLDEG//2
    logSlots = ZZ(slots).nbits() - 1
    logNh = ZZ(Nh).nbits() - 1
    logT = int(logT)
    logI = int(logI)
    repeats = logI + logT

    if logSlots == 0 and not cipher.is_complex:
        out = hb_imult(cipher)
        out = hb_div_by_po2(out, logT)
        out = heaan_exp2pi(out, relin_key)
        for _ in range(repeats):
            out = hb_square_rescale(out, relin_key)
        tmp = hb_conjugate(out, boot_keys)
        out = hb_sub(out, tmp)
        out = hb_idiv(out)
        out = hb_mult_by_const(out, HB_RR(1)/(HB_RR(4)*HB_RR.pi()))

    elif logSlots < logNh:
        tmp = hb_conjugate(cipher, boot_keys)
        out = hb_sub(cipher, tmp)
        out = hb_div_by_po2(out, logT + 1)
        out = heaan_exp2pi(out, relin_key)
        for _ in range(repeats):
            out = hb_square_rescale(out, relin_key)
        tmp = hb_conjugate(out, boot_keys)
        out = hb_sub(out, tmp)

        # p1 branch
        pscale = ZZ(q_list[out.level()])
        p1poly = _hb_scale_unit_poly(boot_context.p1_unit, pscale)
        tmp = hb_mult_by_poly(out, p1poly, pscale)
        tmprot = hb_left_rotate_fast(tmp, slots, boot_keys)
        tmp = hb_add(tmp, tmprot)

        # p2 branch
        p2poly = _hb_scale_unit_poly(boot_context.p2_unit, pscale)
        out = hb_mult_by_poly(out, p2poly, pscale)
        tmprot = hb_left_rotate_fast(out, slots, boot_keys)
        out = hb_add(out, tmprot)
        out = hb_add(out, tmp)

    else:
        tmp = hb_conjugate(cipher, boot_keys)
        c2 = hb_sub(cipher, tmp)
        out = hb_add(cipher, tmp)
        out = hb_imult(out)
        out = hb_div_by_po2(out, logT + 1)

        # C++ uses reScaleBy(c2, logT+1), which divides both coefficients and
        # metadata scale by the same small power of two.  No RNS prime is spent.
        c2 = hb_virtual_rescale_by_po2(c2, logT + 1)

        out = heaan_exp2pi(out, relin_key, boot_plain_scale)
        c2 = heaan_exp2pi(c2, relin_key, boot_plain_scale)
        for _ in range(repeats):
            c2 = hb_square_rescale(c2, relin_key)
            out = hb_square_rescale(out, relin_key)

        tmp = hb_conjugate(c2, boot_keys)
        c2 = hb_sub(c2, tmp)
        tmp = hb_conjugate(out, boot_keys)
        out = hb_sub(out, tmp)
        out = hb_imult(out)

        # Ring2Utils::subAndEqual2(c2, out) means out <- c2 - out.
        #out = hb_sub(c2, out)
        out = hb_sub_into_second(c2, out)

        out = hb_mult_by_const(out, HB_RR(1)/(HB_RR(4)*HB_RR.pi()))

    # C++ reScaleBy(logp + logI): in RNS, the multiplication above contributes
    # one top-q plaintext scale.  Consume that q-prime, then emulate the extra
    # 2^logI factor as a small virtual rescale.
    out = hb_rescale_one(out)
    if logI > 0:
        out = hb_virtual_rescale_by_po2(out, logI)
    return out



# ModRaise and bootstrapp

def _hb_center_lift_q0_int(a):
    a = ZZ(a) % ZZ(q_list[0])
    if a > ZZ(q_list[0])//2:
        a -= ZZ(q_list[0])
    return a


def _hb_modraise_component_q0(poly_q0, active):
    coeffs = _dense_coeffs(poly_q0, POLDEG)
    z = R([_hb_center_lift_q0_int(c) for c in coeffs]) % firreducible
    return [poly_mod_prime(z, q_list[j]) for j in range(active)]


def heaan_mod_raise_q0(raw_q0, active=None):

    if active is None:
        active = L
    active = int(active)
    if len(raw_q0[0]) != 1:
        raw_q0 = _hb_raw_mod_down_to_level(raw_q0, 0)

    c0 = _hb_modraise_component_q0(raw_q0[0][0], active)
    c1 = []
    for comp in raw_q0[1]:
        c1.append(_hb_modraise_component_q0(comp[0], active))
    return [c0, c1]


def client_encrypt_q0_poly(m_poly, key):
    """q0-only ModHE encryption of an already encoded CKKS polynomial."""
    rank = key[0].ncols()
    e0 = random_poly_error()
    e1 = matrix(R, 1, rank)
    v = matrix(R, 1, rank)
    for k in range(rank):
        e1[0, k] = random_poly_error()
        v[0, k] = random_poly_sparse_binary()

    b_pk = key[1]
    A = key[2]
    temp0 = (v*b_pk + m_poly + e0) % firreducible
    c0 = poly_mod_prime(temp0[0, 0], q_list[0])

    temp1 = (v*A + e1) % firreducible
    c1 = [[poly_mod_prime(temp1[0, k], q_list[0])] for k in range(rank)]
    return [[c0], c1]


def client_encrypt_q0_ckks(values, key, scale=p_enc, method="fft"):
    mpoly = ckks_encode(values, scale, method)
    return client_encrypt_q0_poly(mpoly, key)


def heaan_subsum(cipher, boot_keys):
    """HEAANBOOT SubSum exactly as in bootstrapAndEqual()."""
    out = cipher.copy()
    Nh = POLDEG//2
    logSlots = ZZ(out.slots).nbits() - 1
    logNh = ZZ(Nh).nbits() - 1

    for i in range(int(logSlots), int(logNh)):
        rot = hb_left_rotate_fast(out, 1 << i, boot_keys)
        out = hb_add(out, rot)

    if logSlots == 0 and not out.is_complex:
        cconj = hb_conjugate(out, boot_keys)
        out = hb_add(out, cconj)
        out = hb_div_by_po2(out, ZZ(POLDEG).nbits()-1)
    else:
        out = hb_div_by_po2(out, int(logNh))
    return out


def heaan_bootstrap(cipher, boot_context, boot_keys, relin_key,
                    logT=0, logI=4, output_scale=None,
                    boot_logp=None, active_raise=None, verbose=True):
    """
    Full HEAANBOOT bootstrap adapted to ModHE's MLWE/RNS ciphertexts.
    """
    if not isinstance(cipher, ModHEBootCipher):
        raise TypeError("heaan_bootstrap expects ModHEBootCipher; use hb_wrap()")
    if output_scale is None:
        output_scale = cipher.scale
    if boot_logp is None:
        boot_logp = ZZ(q_list[0]).nbits() - 1 + 4
    if active_raise is None:
        active_raise = L

    # HEAAN: modDownTo(logq), normalize, then reinterpret under large Q.
    q0raw = _hb_raw_mod_down_to_level(cipher.raw, 0)
    raised = heaan_mod_raise_q0(q0raw, active=active_raise)
    # C++ uses q=2^logq and then sets logp=logq+4, i.e. temporary scale=16*q.
    # For the prime q0 used by ModHE, 16*q0 is the exact corresponding scale.
    boot_scale = HB_RR(16 * ZZ(q_list[0]))
    out = ModHEBootCipher(raised, boot_scale,
                          cipher.slots, cipher.is_complex)

    if verbose:
        print("[bootstrap] after ModRaise:", out)

    out = heaan_subsum(out, boot_keys)
    if verbose:
        print("[bootstrap] after SubSum:", out)

    logSlots = ZZ(out.slots).nbits() - 1
    if logSlots == 0 and not out.is_complex:
        out = heaan_eval_exp(out, boot_context, boot_keys, relin_key, logT, logI)
    else:
        out = heaan_coeff_to_slot(out, boot_context, boot_keys)
        if verbose:
            print("[bootstrap] after CoeffToSlot:", out)
        out = heaan_eval_exp(out, boot_context, boot_keys, relin_key, logT, logI)
        if verbose:
            print("[bootstrap] after EvalExp:", out)
        out = heaan_slot_to_coeff(out, boot_context, boot_keys)

    # Exact analogue of `cipher.logp = logp` at the end of C++ bootstrap.
    out.scale = HB_RR(output_scale)
    if verbose:
        print("[bootstrap] after SlotToCoeff/final scale reset:", out)
    return out


# Some checks

def heaan_bootstrap_depth_report(logT=0, logI=4, packed=True):
    if packed:
        need = int(logT) + int(logI) + 9
    else:
        need = int(logT) + int(logI) + 4
    print("HEAANBOOT RNS depth report")
    print("  active q-shares available =", L)
    print("  rescale primes required    =", need)
    print("  shares required incl q0    =", need + 1)
    print("  enough chain?              =", bool(L >= need+1))
    return need


def hb_decrypt_slots(cipher, key, method="fft"):
    """Decode using the wrapper's current scale, directly from q0 share."""
    return decryption_ckks(cipher.raw, key, scale=cipher.scale,
                           slots=cipher.slots, method=method)

