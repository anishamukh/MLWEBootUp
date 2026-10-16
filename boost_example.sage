# -----------------------------------------------------------------------------
# client_rankup_bootstrap_test.sage
#
# Correctness checks for the hypothesis in ModHE/MLWE-CKKS:
#   1. client encrypts at low rank r_c (one in our case) and only one base modulus q_0;
#   2. server performs a proper key-switching rank-up r_c -> r_E;
# -----------------------------------------------------------------------------

_sage_const_2 = Integer(2); _sage_const_0p4 = RealNumber('0.4'); _sage_const_0p6 = RealNumber('0.6'); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_4 = Integer(4)
load("constants.sage")
load("auxiliary_functions.sage")
load("KeyGen_rank_up_down.sage")
load("Encrypt-Decrypt_rank_up_down.sage")
load("HMult.sage")
load("Relin_rank_up_down.sage")
load("Flexibility.sage")
load("boost_mlwe.sage")

def _center_lift_int(a, q):
    a = int(a) % int(q)
    if a > int(q)//2:
        a -= int(q)
    return a


def _center_lift_poly_from_q0(poly_q0):
    coeffs = _dense_coeffs(poly_q0, POLDEG)
    return R([_center_lift_int(c, q_list[0]) for c in coeffs])


def _shares_from_center_lift_q0(poly_q0, active=L):
    z = _center_lift_poly_from_q0(poly_q0)
    return [poly_mod_prime(z, q_list[j]) for j in range(active)]


def toy_modraise_q0_to_full(ct_q0, active=L):
    if len(ct_q0[0]) != 1:
        raise ValueError("toy_modraise_q0_to_full expects a q0-only ciphertext")

    c0_full = _shares_from_center_lift_q0(ct_q0[0][0], active)
    c1_full = []
    for comp in ct_q0[1]:
        if len(comp) != 1:
            raise ValueError("toy_modraise_q0_to_full expects q0-only c1 components")
        c1_full.append(_shares_from_center_lift_q0(comp[0], active))
    return [c0_full, c1_full]


def client_encrypt_q0(m_plain, key):
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


def ct_shape(ct):
    return (len(ct[1]), len(ct[0]))  # (rank, active q-shares)


def check_close(label, got, expected, tol=1e-5):
    diff = abs(float(got) - float(expected))
    print(label, "=", got, " expected=", expected, " diff=", diff)
    if diff > tol:
        print("WARNING:", label, "diff exceeds tolerance", tol)
    return diff


# -----------------------------------------------------------------------------
# Demo parameters.
# You can set r_client=1 and r_eval=2 or 4.  For fast debugging, reduce POLDEG in
# constants.sage.  POLDEG=2^14 is the intended base-ring setting but may be slow.
# -----------------------------------------------------------------------------
r_client = 1
r_eval = 4
m_plain = 0.4

print("\n=== Client q0-only encryption -> rank-up -> toy ModRaise correctness test ===")
print("POLDEG=", POLDEG, " L=", L, " r_client=", r_client, " r_eval=", r_eval)

# Independent client and evaluation-ready keys.  The rank-up key switches from the
# client secret to the evaluation secret; it is not zero-padding.
key_client = public_private_keygen(rank=r_client)
key_eval = public_private_keygen(rank=r_eval)

ct_client_q0 = client_encrypt_q0(m_plain, key_client)
print(ct_client_q0)
print("ct_client_q0 shape (rank, shares)=", ct_shape(ct_client_q0))
dec_client = decryption(ct_client_q0, key_client, 0, 0)
check_close("client q0 decrypt", dec_client, m_plain)

rup_key = rankup_keygen(key_client[0], key_eval[0])
ct_eval_q0 = rank_up(ct_client_q0, rup_key)
print("ct_eval_q0 shape after rank-up (rank, shares)=", ct_shape(ct_eval_q0))
dec_rankup = decryption(ct_eval_q0, key_eval, 0, 0)
check_close("after rank-up decrypt", dec_rankup, m_plain)

ct_eval_full = toy_modraise_q0_to_full(ct_eval_q0, active=L)
print("ct_eval_full shape after toy ModRaise (rank, shares)=", ct_shape(ct_eval_full))

# Decryption is still checked at q0, as in ModHE/CKKS decryption.  The extra
# shares are for the subsequent real bootstrapping/evaluation pipeline.
dec_full_q0 = decryption(ct_eval_full, key_eval, 0, 0)
check_close("after toy ModRaise decrypt at q0", dec_full_q0, m_plain)

# Sanity check that the q0 share is preserved by toy ModRaise.
print("q0 c0 share preserved:", ct_eval_full[0][0] == ct_eval_q0[0][0])
print("q0 c1 shares preserved:", all(ct_eval_full[1][k][0] == ct_eval_q0[1][k][0] for k in range(r_eval)))

print("\nNOTE: This validates correctness preservation only.  It does not validate")
print("MLWE bootstrapping security/noise or post-bootstrap refresh. They require")
print("a real MLWE-CKKS bootstrap implementation.")
