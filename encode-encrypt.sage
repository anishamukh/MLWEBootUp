# encode-encrypt.sage

def encryption_poly(m_poly, key):
    """Encrypt an already-encoded plaintext polynomial m_poly in the base ring."""
    rank = key[0].ncols()

    e0 = matrix(R, 1, rank)
    e1 = matrix(R, 1, rank)
    v = matrix(R, 1, rank)
    for k in range(rank):
        e0[0, k] = random_poly_error()
        e1[0, k] = random_poly_error()
        v[0, k] = random_poly_sparse_binary()

    b_pk = key[1]
    A = key[2]

    temp0 = (v*b_pk + m_poly + e0[0, 0]) % firreducible
    c0 = poly_mod_prime(temp0[0, 0], Q)

    temp1 = (v*A + e1) % firreducible
    c1 = [poly_mod_prime(temp1[0, k], Q) for k in range(rank)]

    c0_qshares = shares_in_mod(c0, L, q_list)
    c1_qshares = [shares_in_mod(c1[k], L, q_list) for k in range(rank)]

    return [c0_qshares, c1_qshares]


def encryption_ckks(slots, key, scale=p_enc, method="fft"):
    m_poly = ckks_encode(slots, scale, method)
    return encryption_poly(m_poly, key)


def decrypt_poly_q0(ct, key, level1=0, level2=0):
    c0_q = ct[0]
    c1_q = ct[1]
    s = key[0]
    rank = s.ncols()

    if level2 > level1:
        raise ValueError("level1 must be >= level2")

    if level1 > level2:
        c0_dec_q_shares = rescaling_ciphertext(c0_q, level1, level2)
        c1_dec_q_shares = [rescaling_ciphertext(c1_q[k], level1, level2) for k in range(rank)]
    else:
        c0_dec_q_shares = c0_q
        c1_dec_q_shares = c1_q

    temp = c0_dec_q_shares[0]
    for k in range(rank):
        temp = (temp + c1_dec_q_shares[k][0]*s[0, k]) % firreducible

    dec_q0 = poly_mod_prime(temp, q_list[0])
    return center_lift_poly_mod_q(dec_q0, q_list[0])


def decryption_ckks(ct, key, level1=0, level2=0, scale=p_enc, slots=None, method="fft"):
    """Decrypt a ModHE ciphertext and CKKS-decode the base-ring polynomial."""
    m_poly = decrypt_poly_q0(ct, key, level1, level2)
    return ckks_decode(m_poly, scale, slots, method)


def decryption_constant_compat(ct, key, level1=0, level2=0, scale=p_enc):
    m_poly = decrypt_poly_q0(ct, key, level1, level2)
    coeffs = _dense_coeffs(m_poly, POLDEG)
    return RRK(coeffs[0]) / RRK(scale)
