# function to generate PUBLIC KEY for arbitrary MLWE rank

def public_private_keygen():

    # Generation of secret polynomial vector s=[s0, s1, ...]
    s = matrix(R, 1, MLWE)
    for k in range(MLWE):
        s[0, k] = random_poly_binary(force_monic_last=True)

    # Generation of error polynomial vector e=[e0, e1, ...] for encryption
    e = matrix(R, 1, MLWE)
    for k in range(MLWE):
        e[0, k] = random_poly_error()

    # Generate a random matrix A in R_Q^{MLWE x MLWE}
    A = matrix_vector_random_uniform(MLWE, MLWE, Q)

    b_pk = -A*s.transpose() + e.transpose()

    # sk = {1, s}, pk = {b_pk, A}
    key = [s, b_pk, A]
    print('private-public keypair generated for MLWE rank', MLWE)

    return key
