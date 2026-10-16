# -----------------------------------------------------------------------------
# Flexibility.sage
# General rank-up and rank-reduction/rank-down routines for the ModHE Sage
# prototype.
#
# rank_red is the ModHE rank-reduction operation: it keeps the first r_tgt
# secret/ciphertext components and key-switches the removed components
# s_{r_tgt},...,s_{r_src-1} under the retained secret.
# -----------------------------------------------------------------------------


def _secret_to_list(secret_or_key):
    try:
        if isinstance(secret_or_key, (list, tuple)) and len(secret_or_key) >= 3:
            s = secret_or_key[0]
            return [s[0, i] for i in range(s.ncols())]
    except Exception:
        pass

    try:
        return [secret_or_key[0, i] for i in range(secret_or_key.ncols())]
    except Exception:
        pass

    if isinstance(secret_or_key, (list, tuple)):
        return list(secret_or_key)

    return [secret_or_key]


def _poly_coeff(poly, i):
    coeffs = poly.coefficients(sparse=False)
    if i < len(coeffs):
        return coeffs[i]
    return 0


def _share_index(out_idx, active):
    if out_idx < active:
        return out_idx
    return L


def _mod_at_output_index(out_idx, active):
    if out_idx < active:
        return q_list[out_idx]
    return p


def _zero_poly_shares(size):
    return [R(0) for _ in range(size)]


def _moddown_poly_shares(poly_shares_q_plus_p, active):
    return ModDown_poly(poly_shares_q_plus_p, active)


def extend_key_with_extra_secret(key_prefix, target_rank):
    prefix = _secret_to_list(key_prefix)
    if target_rank < len(prefix):
        raise ValueError("extend_key_with_extra_secret: target_rank must be >= prefix rank")

    secret_polys = list(prefix)
    for _ in range(target_rank - len(prefix)):
        secret_polys.append(random_poly_binary(force_monic_last=True))

    return public_private_keygen(rank=target_rank, secret_polys=secret_polys)


# Rank-up / rank-lifting key switching

def rankup_keygen(s_src, s_tgt):
    s_src_list = _secret_to_list(s_src)
    print("\nsource secret=", s_src_list)

    s_tgt_list = _secret_to_list(s_tgt)
    print("\ntarget secret=", s_tgt_list)
    
    r_src = len(s_src_list)
    r_tgt = len(s_tgt_list)

    b_matrix = []
    a_matrix = []

    for u in range(r_src):
        b_for_u = []
        a_for_u_by_target = [[] for _ in range(r_tgt)]

        for j in range(L):
            a_vec = [random_poly_uniform(p*Q) for _ in range(r_tgt)]
            err = random_poly_error()

            inner = R(0)
            for t in range(r_tgt):
                inner = inner + a_vec[t]*s_tgt_list[t]

            b_poly = (-inner + err + p*B[j]*s_src_list[u]) % firreducible
            b_poly = poly_mod_prime(b_poly, p*Q)
            b_for_u.append(shares_in_mod(b_poly, L+1, moduli))

            for t in range(r_tgt):
                a_poly = poly_mod_prime(a_vec[t] % firreducible, p*Q)
                a_for_u_by_target[t].append(shares_in_mod(a_poly, L+1, moduli))

        b_matrix.append(b_for_u)
        a_matrix.append(a_for_u_by_target)

    print('rank-up key generated:', r_src, '->', r_tgt)
    return [b_matrix, a_matrix]


def rank_up(ciphertext, rankup_key):
    c0_qshares = ciphertext[0]
    c1_qshares = ciphertext[1]
    active = len(c0_qshares)

    b_matrix = rankup_key[0]
    a_matrix = rankup_key[1]
    r_src = len(b_matrix)
    r_tgt = len(a_matrix[0])

    if len(c1_qshares) != r_src:
        raise ValueError("rank_up: ciphertext rank and rank-up key source rank do not match")

    ks0 = _zero_poly_shares(active + 1)
    ks1 = [[R(0) for _ in range(active + 1)] for _ in range(r_tgt)]

    for out_idx in range(active + 1):
        mod = _mod_at_output_index(out_idx, active)
        key_idx = _share_index(out_idx, active)

        acc0 = R(0)
        acc1 = [R(0) for _ in range(r_tgt)]

        for u in range(r_src):
            for j in range(active):
                c = c1_qshares[u][j]
                acc0 = acc0 + poly_mod_prime((c*b_matrix[u][j][key_idx]) % firreducible, mod)
                for t in range(r_tgt):
                    acc1[t] = acc1[t] + poly_mod_prime((c*a_matrix[u][t][j][key_idx]) % firreducible, mod)

        ks0[out_idx] = poly_mod_prime(acc0, mod)
        for t in range(r_tgt):
            ks1[t][out_idx] = poly_mod_prime(acc1[t], mod)

    ks0_scaled = _moddown_poly_shares(ks0, active)
    ks1_scaled = [_moddown_poly_shares(ks1[t], active) for t in range(r_tgt)]

    out0 = []
    for j in range(active):
        out0.append(poly_mod_prime(c0_qshares[j] + ks0_scaled[j], q_list[j]))

    out1 = []
    for t in range(r_tgt):
        out1.append([poly_mod_prime(ks1_scaled[t][j], q_list[j]) for j in range(active)])

    return [out0, out1]

def rankup_keygen_from_client_secret(s_client, s_eval):
    return rankup_keygen([s_client], s_eval)

# Rank reduction / rank down

def rankred_keygen(source_key_or_secret, target_rank):
    s_src = _secret_to_list(source_key_or_secret)
    r_src = len(s_src)

    if target_rank >= r_src:
        raise ValueError("rankred_keygen: target_rank must be smaller than source rank")
    if target_rank < 1:
        raise ValueError("rankred_keygen: target_rank must be at least 1")

    s_tgt = s_src[:target_rank]
    removed_indices = list(range(target_rank, r_src))

    b_matrix = []
    a_matrix = []

    for rem in removed_indices:
        b_for_rem = []
        a_for_rem_by_target = [[] for _ in range(target_rank)]

        for j in range(L):
            a_vec = [random_poly_uniform(p*Q) for _ in range(target_rank)]
            err = random_poly_error()

            inner = R(0)
            for t in range(target_rank):
                inner = inner + a_vec[t]*s_tgt[t]

            b_poly = (-inner + err + p*B[j]*s_src[rem]) % firreducible
            b_poly = poly_mod_prime(b_poly, p*Q)
            b_for_rem.append(shares_in_mod(b_poly, L+1, moduli))

            for t in range(target_rank):
                a_poly = poly_mod_prime(a_vec[t] % firreducible, p*Q)
                a_for_rem_by_target[t].append(shares_in_mod(a_poly, L+1, moduli))

        b_matrix.append(b_for_rem)
        a_matrix.append(a_for_rem_by_target)

    print('rank-reduction key generated:', r_src, '->', target_rank)
    return [target_rank, removed_indices, b_matrix, a_matrix]

def rankred_key(key, target_rank=1):
    return rankred_keygen(key, target_rank)


def rank_red(ciphertext, rankred_key):
    c0_qshares = ciphertext[0]
    c1_qshares = ciphertext[1]
    active = len(c0_qshares)
    r_src = len(c1_qshares)

    target_rank = rankred_key[0]
    removed_indices = rankred_key[1]
    b_matrix = rankred_key[2]
    a_matrix = rankred_key[3]

    if target_rank >= r_src:
        raise ValueError("rank_red: target rank must be smaller than ciphertext rank")
    if removed_indices != list(range(target_rank, r_src)):
        raise ValueError("rank_red: key is not compatible with prefix rank reduction")

    ks0 = _zero_poly_shares(active + 1)
    ks1 = [[R(0) for _ in range(active + 1)] for _ in range(target_rank)]

    for out_idx in range(active + 1):
        mod = _mod_at_output_index(out_idx, active)
        key_idx = _share_index(out_idx, active)

        acc0 = R(0)
        acc1 = [R(0) for _ in range(target_rank)]

        for u, rem in enumerate(removed_indices):
            for j in range(active):
                c = c1_qshares[rem][j]
                acc0 = acc0 + poly_mod_prime((c*b_matrix[u][j][key_idx]) % firreducible, mod)
                for t in range(target_rank):
                    acc1[t] = acc1[t] + poly_mod_prime((c*a_matrix[u][t][j][key_idx]) % firreducible, mod)

        ks0[out_idx] = poly_mod_prime(acc0, mod)
        for t in range(target_rank):
            ks1[t][out_idx] = poly_mod_prime(acc1[t], mod)

    ks0_scaled = _moddown_poly_shares(ks0, active)
    ks1_scaled = [_moddown_poly_shares(ks1[t], active) for t in range(target_rank)]

    out0 = []
    for j in range(active):
        out0.append(poly_mod_prime(c0_qshares[j] + ks0_scaled[j], q_list[j]))

    out1 = []
    for t in range(target_rank):
        comp = []
        for j in range(active):
            comp.append(poly_mod_prime(c1_qshares[t][j] + ks1_scaled[t][j], q_list[j]))
        out1.append(comp)

    return [out0, out1]
