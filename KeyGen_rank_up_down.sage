# function to generate PUBLIC KEY for arbitrary MLWE rank
# Optional arguments make it possible to create several compatible keys in the
# same example, e.g., a rank-2 key and a rank-4 key whose first two secret
# components are the rank-2 secret.

def public_private_keygen(rank=None, secret_polys=None):

    if rank is None:
        rank = MLWE

    if secret_polys is not None:
        if len(secret_polys) != rank:
            raise ValueError("public_private_keygen: len(secret_polys) must equal rank")

    # Generation of secret polynomial vector s=[s0, s1, ...]
    s = matrix(R, 1, rank)
    for k in range(rank):
        if secret_polys is None:
            s[0, k] = random_poly_binary(force_monic_last=True)
        else:
            s[0, k] = secret_polys[k]

    # Generation of error polynomial vector e=[e0, e1, ...] for encryption
    e = matrix(R, 1, rank)
    for k in range(rank):
        e[0, k] = random_poly_error()

    # Generate a random matrix A in R_Q^{rank x rank}
    A = matrix_vector_random_uniform(rank, rank, Q)

    b_pk = -A*s.transpose() + e.transpose()

    # sk = {1, s}, pk = {b_pk, A}
    key = [s, b_pk, A]
    print('private-public keypair generated for MLWE rank', rank)

    return key
