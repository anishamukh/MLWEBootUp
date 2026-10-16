# function for ENCRYPTION for arbitrary MLWE rank

def encryption(m_plain, key):
    e0 = matrix(R, 1, MLWE)
    for k in range(MLWE):
        e0[0, k] = random_poly_error()

    # Step2: generate error vector e1
    e1 = matrix(R, 1, MLWE)
    for k in range(MLWE):
        e1[0, k] = random_poly_error()

    # Step3: generate sparse binary vector v
    v = matrix(R, 1, MLWE)
    for k in range(MLWE):
        v[0, k] = random_poly_sparse_binary()

    # Step4: encode one real number as the constant coefficient
    m_enc = round(m_plain*p_enc)
    print("m=", m_plain)

    b_pk = key[1]
    A = key[2]

    temp0 = (v*b_pk + m_enc + e0[0, 0]) % firreducible
    c0 = poly_mod_prime(temp0[0, 0], Q)

    temp1 = (v*A + e1) % firreducible
    c1 = [poly_mod_prime(temp1[0, k], Q) for k in range(MLWE)]

    c0_qshares = shares_in_mod(c0, L, q_list)
    c1_qshares = [shares_in_mod(c1[k], L, q_list) for k in range(MLWE)]

    return [c0_qshares, c1_qshares]


# function for DECRYPTION for arbitrary MLWE rank

def decryption(ct, key, level1, level2):

    print("decryption:")
    c0_q = ct[0]
    c1_q = ct[1]
    s = key[0]

    if level2 > level1:
        print('ERROR: level1 must be >= level2')
        exit()

    if level1 > level2:
        c0_dec_q_shares = rescaling_ciphertext(c0_q, level1, level2)
        c1_dec_q_shares = [rescaling_ciphertext(c1_q[k], level1, level2) for k in range(MLWE)]
    else:
        c0_dec_q_shares = c0_q
        c1_dec_q_shares = c1_q

    temp = c0_dec_q_shares[0]
    for k in range(MLWE):
        temp = (temp + c1_dec_q_shares[k][0]*s[0, k]) % firreducible

    dec = poly_mod_prime(temp, q_list[0])
    dec_coeff = dec.coefficients(sparse=False)
    dec_m = float(dec_coeff[0]/p_enc)
    return dec_m
