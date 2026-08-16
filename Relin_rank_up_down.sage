# Relinearisation for arbitrary MLWE rank.
# Evaluation key layout:
#   evk_key = [pairs, evkb_shares, Arow_shares]
#   pairs[t] = (i,k) represents monomial s_i*s_k
#   evkb_shares[t][j][u] is b-share for pair t, decomposition modulus q_j,
#                         output modulus moduli[u]
#   Arow_shares[t][j][k][u] is A-row share for output linear component k.

# function for EVALUATION KEY
def relin_keygen(key):
    s = key[0]
    rank = s.ncols()
    pairs = _quadratic_pairs(rank)

    evkb_shares = []
    Arow_shares = []

    for (u, v) in pairs:
        pair_b = []
        pair_A = []
        monomial = (s[0, u]*s[0, v]) % firreducible

        for j in range(L):
            Arow = matrix(R, 1, rank)
            for k in range(rank):
                Arow[0, k] = random_poly_uniform(p*Q)

            err = random_poly_error()
            b = (-(Arow*s.transpose())[0, 0] + p*B[j]*monomial + err) % firreducible
            b = poly_mod_prime(b, p*Q)

            pair_b.append(shares_in_mod(b, L+1, moduli))
            pair_A.append([shares_in_mod(Arow[0, k], L+1, moduli) for k in range(rank)])

        evkb_shares.append(pair_b)
        Arow_shares.append(pair_A)

    print('relinearization key generated for', len(pairs), 'quadratic monomials')
    return [pairs, evkb_shares, Arow_shares]


# function for changing modulus (MOD DOWN)
def ModDown(b_share_p, b_shares_q, active=None):
    if active is None:
        active = L
    result_q = [0]*active
    for i in range(active):
        a_share_q = b_share_p % q_list[i]
        result_q[i] = p_inv_q[i]*(b_shares_q[i] - a_share_q) % q_list[i]
    return result_q


def ModDown_poly(poly_shares_q_plus_p, active):
    """Apply ModDown to a polynomial represented in q_0...q_level plus p."""
    q_coeffs = [_dense_coeffs(poly_shares_q_plus_p[j], POLDEG) for j in range(active)]
    p_coeffs = _dense_coeffs(poly_shares_q_plus_p[active], POLDEG)
    out_coeffs = [[0]*POLDEG for _ in range(active)]

    for i in range(POLDEG):
        shares_q = [int(q_coeffs[j][i]) for j in range(active)]
        share_p = int(p_coeffs[i])
        scaled = ModDown(share_p, shares_q, active)
        for j in range(active):
            out_coeffs[j][i] = int(scaled[j])

    return [R(out_coeffs[j]) for j in range(active)]

def relinearisation(non_linearized_ciphertext, evk_key):

    d0_qshares = non_linearized_ciphertext[0]
    d1_qshares = non_linearized_ciphertext[1]
    d2_qshares = non_linearized_ciphertext[2]

    rank = len(d1_qshares)

    # Number of q_i currently present in the ciphertext.
    # E.g. level 21 -> active = 22 -> q_0,...,q_21.
    active = len(d0_qshares)

    pairs = evk_key[0]
    evkb_shares = evk_key[1]
    Arow_shares = evk_key[2]

    if len(d2_qshares) != len(pairs):
        raise ValueError(
            "ciphertext quadratic part and evaluation-key pairs do not match"
        )

    # ------------------------------------------------------------
    # IMPORTANT:
    #
    # Evaluation keys were generated over the FULL modulus list
    #
    #     q_0,...,q_{L-1},p
    #
    # and therefore the p-share is permanently stored at the final
    # key index, irrespective of the current ciphertext level.
    # ------------------------------------------------------------
    p_key_idx = len(moduli) - 1

    # Temporary p-extended result has layout
    #
    #     q_0,...,q_{active-1},p
    #
    # and therefore has active+1 entries.
    ct0_prime = [R(0) for _ in range(active + 1)]
    ct1_prime = [
        [R(0) for _ in range(active + 1)]
        for _ in range(rank)
    ]

    for out_pos in range(active + 1):

        # --------------------------------------------------------
        # Distinguish:
        #
        # out_pos:
        #     index in the SHORTENED temporary result
        #
        # key_idx:
        #     index in the FULL evaluation key.
        # --------------------------------------------------------
        if out_pos < active:

            # q_0,...,q_{active-1}
            mod = q_list[out_pos]
            key_idx = out_pos

        else:

            # Last temporary position is p.
            mod = p

            # But the evaluation-key p-share is always at the
            # FINAL position of the full key.
            key_idx = p_key_idx

        acc0 = R(0)
        acc1 = [R(0) for _ in range(rank)]

        for t in range(len(pairs)):

            for j in range(active):

                # j indexes the active RNS decomposition digit.
                d = d2_qshares[t][j]

                acc0 += poly_mod_prime(
                    (
                        d *
                        evkb_shares[t][j][key_idx]
                    ) % firreducible,
                    mod
                )

                for k in range(rank):

                    acc1[k] += poly_mod_prime(
                        (
                            d *
                            Arow_shares[t][j][k][key_idx]
                        ) % firreducible,
                        mod
                    )

        ct0_prime[out_pos] = poly_mod_prime(
            acc0,
            mod
        )

        for k in range(rank):
            ct1_prime[k][out_pos] = poly_mod_prime(
                acc1[k],
                mod
            )

    # ------------------------------------------------------------
    # ModDown:
    #
    #   (q_0,...,q_{active-1},p)
    #            ->
    #   (q_0,...,q_{active-1})
    # ------------------------------------------------------------
    ct0_scaled = ModDown_poly(
        ct0_prime,
        active
    )

    ct1_scaled = [
        ModDown_poly(ct1_prime[k], active)
        for k in range(rank)
    ]

    # Add the original linear portion of the product.
    cmult0_qshares = [None] * active

    cmult1_qshares = [
        [None] * active
        for _ in range(rank)
    ]

    for j in range(active):

        cmult0_qshares[j] = poly_mod_prime(
            d0_qshares[j] + ct0_scaled[j],
            q_list[j]
        )

        for k in range(rank):

            cmult1_qshares[k][j] = poly_mod_prime(
                d1_qshares[k][j] + ct1_scaled[k][j],
                q_list[j]
            )

    return [
        cmult0_qshares,
        cmult1_qshares
    ]

# function for RELINEARISATION
def relinearisation_old(non_linearized_ciphertext, evk_key):

    d0_qshares = non_linearized_ciphertext[0]
    d1_qshares = non_linearized_ciphertext[1]
    d2_qshares = non_linearized_ciphertext[2]

    rank = len(d1_qshares)
    active = len(d0_qshares)     # number of q-prime shares currently active
    # q_list is a Sage vector in constants.sage.  Convert the slice to a
    # Python list before appending the special prime p; otherwise Sage tries
    # to add a Vector_integer_dense and a Python list.
    out_moduli = list(q_list[:active]) + [p]

    pairs = evk_key[0]
    evkb_shares = evk_key[1]
    Arow_shares = evk_key[2]

    if len(d2_qshares) != len(pairs):
        print('ERROR: ciphertext quadratic part and evaluation-key pairs do not match')
        exit()

    # Accumulate the p-extended products before ModDown.
    ct0_prime = [R(0) for _ in range(active+1)]
    ct1_prime = [[R(0) for _ in range(active+1)] for _ in range(rank)]

    """for out_idx in range(active+1):

    # Output representation is:
    #   0,...,active-1 -> q_0,...,q_{active-1}
    #   active         -> special modulus p
        if out_idx < active:
            mod = q_list[out_idx]
            key_idx = out_idx
        else:
            mod = p

            # IMPORTANT:
            # Evaluation keys were generated over
            # q_0,...,q_{L-1},p,
            # so their p-share is permanently stored at index L.
            key_idx = L

        acc0 = R(0)
        acc1 = [R(0) for _ in range(rank)]

        for t in range(len(pairs)):
            for j in range(active):
                d = d2_qshares[t][j]

                acc0 = acc0 + poly_mod_prime(
                    (
                        d *
                        evkb_shares[t][j][key_idx]
                    ) % firreducible,
                    mod
                )

                for k in range(rank):
                    acc1[k] = acc1[k] + poly_mod_prime(
                        (
                            d *
                            Arow_shares[t][j][k][key_idx]
                        ) % firreducible,
                        mod
                    )

    ct0_prime[out_idx] = poly_mod_prime(
        acc0,
        mod
    )

    for k in range(rank):
        ct1_prime[k][out_idx] = poly_mod_prime(
            acc1[k],
            mod
        )"""
    for out_idx in range(active+1):
        mod = out_moduli[out_idx]

        acc0 = R(0)
        acc1 = [R(0) for _ in range(rank)]

        for t in range(len(pairs)):
            for j in range(active):
                d = d2_qshares[t][j]

                acc0 = acc0 + poly_mod_prime((d*evkb_shares[t][j][out_idx]) % firreducible, mod)

                for k in range(rank):
                    acc1[k] = acc1[k] + poly_mod_prime((d*Arow_shares[t][j][k][out_idx]) % firreducible, mod)

        ct0_prime[out_idx] = poly_mod_prime(acc0, mod)
        for k in range(rank):
            ct1_prime[k][out_idx] = poly_mod_prime(acc1[k], mod)

    # Scale the p-extended products back from p*q to q.
    ct0_scaled = ModDown_poly(ct0_prime, active)
    ct1_scaled = [ModDown_poly(ct1_prime[k], active) for k in range(rank)]

    # Add the original degree-0 and degree-1 parts.
    cmult0_qshares = [None]*active
    cmult1_qshares = [[None]*active for _ in range(rank)]

    for j in range(active):
        cmult0_qshares[j] = poly_mod_prime(d0_qshares[j] + ct0_scaled[j], q_list[j])
        for k in range(rank):
            cmult1_qshares[k][j] = poly_mod_prime(d1_qshares[k][j] + ct1_scaled[k][j], q_list[j])

    return [cmult0_qshares, cmult1_qshares]
