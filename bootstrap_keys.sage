# -----------------------------------------------------------------------------
# BootstrapKeys.sage
# Generic ModHE/ModRNS key-switching, rotation, and conjugation keys.
#
# Expected load order:
#   load("constants.sage")
#   load("auxiliary_functions.sage")
#   load("KeyGen.sage")
#   load("Encrypt-Decrypt.sage")
#   load("HMult.sage")
#   load("Relin.sage")
#   load("BootstrapKeys.sage")
#   load("BootstrapOps.sage")
# -----------------------------------------------------------------------------


def _bs_secret_to_list(s_or_key):
    """Return a Python list [s0,...,s_{r-1}] from a key, Sage row matrix, or list."""
    try:
        # full key = [s, b_pk, A]
        if isinstance(s_or_key, (list, tuple)) and len(s_or_key) >= 3 and hasattr(s_or_key[0], 'ncols'):
            s_or_key = s_or_key[0]
    except Exception:
        pass

    try:
        return [s_or_key[0, i] for i in range(s_or_key.ncols())]
    except Exception:
        pass

    if isinstance(s_or_key, (list, tuple)):
        return list(s_or_key)

    return [s_or_key]


def _bs_secret_matrix(secret_list):
    """Convert a list of secret polynomials to a 1 x r Sage matrix."""
    s = matrix(R, 1, len(secret_list))
    for i in range(len(secret_list)):
        s[0, i] = secret_list[i]
    return s


def _bs_uniform_poly_pQ():
    """Uniform-ish polynomial modulo p*Q for switching keys."""
    return random_poly_uniform(p*Q) % firreducible


def ks_keygen(s_from, s_to):
    """Generate a generic ModHE key-switching key.

    It switches ciphertexts decryptable under s_from to ciphertexts decryptable
    under s_to.  For every input component u and RNS digit j it stores an
    encryption of p*B[j]*s_from[u] under the target secret s_to:

        b[u,j] = - <a[u,j], s_to> + e + p*B[j]*s_from[u]  mod pQ.

    Output layout:
        ksk["b"][u][j][ell]
        ksk["a"][u][t][j][ell]

    where ell=0..L-1 is an output q_ell share and ell=L is the special p share.
    """
    sf = _bs_secret_to_list(s_from)
    st = _bs_secret_to_list(s_to)
    r_from = len(sf)
    r_to = len(st)

    b_shares = []
    a_shares = []

    for u in range(r_from):
        b_for_u = []
        a_for_u_by_target = [[] for _ in range(r_to)]

        for j in range(L):
            a_vec = [_bs_uniform_poly_pQ() for _ in range(r_to)]
            inner = R(0)
            for t in range(r_to):
                inner = inner + a_vec[t]*st[t]

            b_poly = (-inner + random_poly_error() + p*B[j]*sf[u]) % firreducible
            b_poly = poly_mod_prime(b_poly, p*Q)
            b_for_u.append(shares_in_mod(b_poly, L+1, moduli))

            for t in range(r_to):
                a_poly = poly_mod_prime(a_vec[t], p*Q)
                a_for_u_by_target[t].append(shares_in_mod(a_poly, L+1, moduli))

        b_shares.append(b_for_u)
        a_shares.append(a_for_u_by_target)

    return {
        "type": "ksk",
        "rank_from": r_from,
        "rank_to": r_to,
        "b": b_shares,
        "a": a_shares,
    }


def rotation_keygen(s_or_key, k):
    """Generate key for kappa_k followed by switching back to the original key."""
    s_list = _bs_secret_to_list(s_or_key)
    s = _bs_secret_matrix(s_list)
    s_auto = _bs_secret_matrix([automorphism_poly(si, k) for si in s_list])
    return {
        "type": "rotation",
        "k": Integer(k) % (2*POLDEG),
        "ksk": ks_keygen(s_auto, s),
    }


def slot_rotation_keygen(s_or_key, r):
    """Generate key for slot rotation by r, i.e. k = 5^r mod 2N."""
    k = power_mod(5, Integer(r), 2*POLDEG)
    return rotation_keygen(s_or_key, k)


def conjugation_keygen(s_or_key):
    """Generate key for CKKS conjugation, i.e. kappa_{-1}."""
    return rotation_keygen(s_or_key, -1)


def boot_keygen(key, rotation_exponents=None, slot_rotation_steps=None, include_relin=True):
    """Bundle the keys used by the traditional CKKS-style bootstrap scaffold.

    rotation_exponents are automorphism exponents k.  slot_rotation_steps are
    slot offsets r and are converted to k=5^r mod 2N.
    """
    if rotation_exponents is None:
        rotation_exponents = []
    if slot_rotation_steps is None:
        slot_rotation_steps = []

    s = key[0]
    rot_keys = {}

    # identity does not need a key
    for k in rotation_exponents:
        kk = Integer(k) % (2*POLDEG)
        if kk != 1:
            rot_keys[kk] = rotation_keygen(s, kk)

    for r in slot_rotation_steps:
        kk = power_mod(5, Integer(r), 2*POLDEG)
        if kk != 1:
            rot_keys[kk] = rotation_keygen(s, kk)

    out = {
        "rot_keys": rot_keys,
        "conj_key": conjugation_keygen(s),
    }

    if include_relin:
        out["relin_key"] = relin_keygen(key)
    else:
        out["relin_key"] = None

    return out
