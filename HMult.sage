# HOMOMORPHIC MULTIPLICATION for arbitrary MLWE rank.
# Output format: [d0, d1, d2]
#   d0: scalar component shares
#   d1[i]: linear component for s_i
#   d2[t]: quadratic component for monomial s_i*s_k, where t indexes
#          _quadratic_pairs(MLWE) = [(0,0),(0,1),...,(r-1,r-1)].

def _quadratic_pairs(rank):
    return [(i, k) for i in range(rank) for k in range(i, rank)]


def homomorphic_mult(ciphertext1, ciphertext2, level, key=None):

    c0_1 = ciphertext1[0]
    c1_1 = ciphertext1[1]

    c0_2 = ciphertext2[0]
    c1_2 = ciphertext2[1]

    rank = len(c1_1)
    active = level + 1

    d0_qshares = [None]*active
    d1_qshares = [[None]*active for _ in range(rank)]

    pairs = _quadratic_pairs(rank)
    d2_qshares = [[None]*active for _ in range(len(pairs))]

    for j in range(active):
        qj = q_list[j]

        d0_qshares[j] = poly_mod_prime((c0_1[j]*c0_2[j]) % firreducible, qj)

        for i in range(rank):
            temp = (c0_1[j]*c1_2[i][j] + c1_1[i][j]*c0_2[j]) % firreducible
            d1_qshares[i][j] = poly_mod_prime(temp, qj)

        for t, (i, k) in enumerate(pairs):
            if i == k:
                temp = (c1_1[i][j]*c1_2[k][j]) % firreducible
            else:
                temp = (c1_1[i][j]*c1_2[k][j] + c1_1[k][j]*c1_2[i][j]) % firreducible
            d2_qshares[t][j] = poly_mod_prime(temp, qj)

    return [d0_qshares, d1_qshares, d2_qshares]


# Auxiliary function: EXTENDED DECRYPTION of the non-linearized product.
# This is generalized mostly for debugging/validation.

def extended_decryption(extended_ciphertext, key, level1, level2):

    d0 = extended_ciphertext[0]
    d1 = extended_ciphertext[1]
    d2 = extended_ciphertext[2]
    s = key[0]
    rank = len(d1)
    pairs = _quadratic_pairs(rank)

    if level2 > level1:
        print('ERROR: level1 must be >= level2')
        exit()

    if level1 > level2:
        d0_q = rescaling_ciphertext(d0, level1, level2)
        d1_q = [rescaling_ciphertext(d1[i], level1, level2) for i in range(rank)]
        d2_q = [rescaling_ciphertext(d2[t], level1, level2) for t in range(len(pairs))]
    else:
        d0_q = d0
        d1_q = d1
        d2_q = d2

    temp = d0_q[0]
    for i in range(rank):
        temp = (temp + d1_q[i][0]*s[0, i]) % firreducible
    for t, (i, k) in enumerate(pairs):
        temp = (temp + d2_q[t][0]*s[0, i]*s[0, k]) % firreducible

    ext_dec = poly_mod_prime(temp, q_list[0])
    print("ext_dec=", ext_dec)
    extdec_coeff = ext_dec.coefficients(sparse=False)
    extdec_mm1 = float(extdec_coeff[0]/p_enc)
    print("mm1=", extdec_mm1)
