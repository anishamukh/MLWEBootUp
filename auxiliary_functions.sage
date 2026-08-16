# Auxiliary functions for ModHE / MLWE-CKKS Sage prototype.
# The routines below avoid coefficient-by-coefficient polynomial construction where
# possible.  This matters at POLDEG = 2^14, where repeated `+ coeff*x^i` loops are
# the main bottleneck of the original proof-of-concept code.

def random_between(j, k):
    return int(random()*(k-j+1)) + j


def _dense_coeffs(poly, n=POLDEG):
    """Return exactly n dense coefficients of a Sage polynomial."""
    try:
        coeffs = poly.coefficients(sparse=False)
    except AttributeError:
        coeffs = [poly]
    if len(coeffs) < n:
        coeffs = coeffs + [0]*(n-len(coeffs))
    return coeffs[:n]


def poly_mod_prime(pol_a, mod):
    """Reduce a polynomial modulo x^N+1 and coefficient-wise modulo mod."""
    try:
        if pol_a.degree() >= POLDEG:
            pol_a = pol_a % firreducible
    except AttributeError:
        pass
    coeffs = _dense_coeffs(pol_a, POLDEG)
    return R([int(c) % mod for c in coeffs])


def random_poly_binary(force_monic_last=False):
    coeffs = [random_between(0, 1) for _ in range(POLDEG)]
    if force_monic_last:
        coeffs[POLDEG-1] = 1
    return R(coeffs)


def random_poly_sparse_binary():
    return R([random_between(0, 1)*random_between(0, 1) for _ in range(POLDEG)])


def random_poly_error():
    return R([random_between(0, Berr) for _ in range(POLDEG)])


def random_poly_uniform(mod):
    # Uses Sage integer sampling instead of floating-point `int(mod*random())`.
    return R([ZZ.random_element(mod) for _ in range(POLDEG)])


# function to evaluate shares
def shares_in_mod(poly_mat, size, mod_list):
    return [poly_mod_prime(poly_mat, mod_list[i]) for i in range(size)]


def matrix_vector_random_error(rows, cols):
    M = matrix(R, rows, cols)
    for i in range(rows):
        for j in range(cols):
            M[i, j] = random_poly_error()
    return M


def matrix_vector_random_uniform(rows, cols, mod):
    M = matrix(R, rows, cols)
    for i in range(rows):
        for j in range(cols):
            M[i, j] = random_poly_uniform(mod)
    return M


# function for rescaling

def rescaling_shares(coeff_shares, level_in, level_out):
    coeff_shares_scaled = copy(coeff_shares)
    for l in range(level_in, level_out, -1):
        for j in range(l):
            temp = (q_list[l]^-1) % q_list[j]
            coeff_shares_scaled[j] = temp*(coeff_shares_scaled[j] - coeff_shares_scaled[l]) % q_list[j]
        coeff_shares_scaled[l] = 0
    return coeff_shares_scaled


def rescaling_ciphertext(ct_shares, level_in, level_out):
    """Fast coefficient-wise modulus switching for one polynomial component."""
    ct_coeffs = [_dense_coeffs(ct_shares[j], POLDEG) for j in range(level_in+1)]
    out_coeffs = [[0]*POLDEG for _ in range(level_out+1)]

    for i in range(POLDEG):
        in_shares = [int(ct_coeffs[j][i]) for j in range(level_in+1)]
        out_shares = rescaling_shares(in_shares, level_in, level_out)
        for j in range(level_out+1):
            out_coeffs[j][i] = int(out_shares[j])

    return [R(out_coeffs[j]) for j in range(level_out+1)]
