
def _center_lift_int_q0(a):
    a = int(a) % int(q_list[0])
    if a > int(q_list[0])//2:
        a -= int(q_list[0])
    return a


def center_lift_poly_from_q0(poly_q0):
    coeffs = _dense_coeffs(poly_q0, POLDEG)
    return R([_center_lift_int_q0(c) for c in coeffs])


def shares_from_center_lift_q0(poly_q0, active=L):
    z = center_lift_poly_from_q0(poly_q0)
    return [poly_mod_prime(z, q_list[j]) for j in range(active)]


def toy_modraise_q0_to_full(ciphertext_q0, active=L):
    if len(ciphertext_q0[0]) != 1:
        raise ValueError("toy_modraise_q0_to_full expects a q0-only ciphertext")

    c0_full = shares_from_center_lift_q0(ciphertext_q0[0][0], active)
    c1_full = []
    for comp in ciphertext_q0[1]:
        if len(comp) != 1:
            raise ValueError("toy_modraise_q0_to_full expects q0-only c1 components")
        c1_full.append(shares_from_center_lift_q0(comp[0], active))
    return [c0_full, c1_full]


def client_encrypt_q0(m_plain, key):
    """Actual client-side low-rank, q0-only encryption.

    This computes only the q0 residue of ModHE encryption.  It intentionally
    avoids constructing the full Q ciphertext, so it models the desired client
    optimization more faithfully than encrypting at Q and slicing q0.
    """
    rank = key[0].ncols()

    e0 = random_poly_error()
    e1 = matrix(R, 1, rank)
    v = matrix(R, 1, rank)
    for k in range(rank):
        e1[0, k] = random_poly_error()
        v[0, k] = random_poly_sparse_binary()

    m_enc = round(m_plain*p_enc)
    b_pk = key[1]
    A = key[2]

    temp0 = (v*b_pk + m_enc + e0) % firreducible
    c0_q0 = poly_mod_prime(temp0[0, 0], q_list[0])

    temp1 = (v*A + e1) % firreducible
    c1_q0 = [[poly_mod_prime(temp1[0, k], q_list[0])] for k in range(rank)]

    return [[c0_q0], c1_q0]