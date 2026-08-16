# -----------------------------------------------------------------------------
# bootstrap_chain.sage
#
# Parameter helper for the HEAANBOOT -> ModHE Sage port.
#
# IMPORTANT
# ---------
# Call configure_heaanboot_chain() *after* loading constants.sage and *before*
# generating any public/evaluation/rotation keys.  It appends RNS primes whose
# size is close to the HEAANBOOT temporary bootstrapping scale 2^(log(q0)+4).
# These extra primes are consumed by CoeffToSlot/EvalExp/SlotToCoeff.
#
# This is a correctness/prototyping parameter helper, NOT a security-certified
# parameter selection routine.  Re-run a concrete MLWE estimator for any real
# parameter set.
# -----------------------------------------------------------------------------


def _hb_floor_log2(n):
    n = ZZ(n)
    if n <= 0:
        raise ValueError("log2 input must be positive")
    return n.nbits() - 1


def _hb_ntt_prime_below(log_target, used=None, start_skip=0):
    """Return a prime q = 1 mod 2*POLDEG just below 2^log_target.

    log_target follows HEAAN's convention: the target magnitude is 2^log_target.
    """
    if used is None:
        used = set()
    step = ZZ(2 * POLDEG)
    target = ZZ(1) << ZZ(log_target)
    k = (target - 2) // step - ZZ(start_skip)
    while k > 1:
        cand = k * step + 1
        if cand not in used and cand.is_prime():
            return cand
        k -= 1
    raise RuntimeError("could not find NTT-friendly prime")


def _hb_ntt_prime_below_target(target, used=None, start_skip=0):
    """Return a prime q = 1 mod 2*POLDEG just below an integer target."""
    if used is None:
        used = set()
    step = ZZ(2 * POLDEG)
    target = ZZ(target)
    if target <= step + 1:
        raise ValueError("target is too small for an NTT-friendly prime")
    k = (target - 2) // step - ZZ(start_skip)
    while k > 1:
        cand = k * step + 1
        if cand not in used and cand.is_prime():
            return cand
        k -= 1
    raise RuntimeError("could not find NTT-friendly prime below target")


def heaanboot_required_rescale_primes(logT, logI=4, packed=True):
    """Number of large RNS primes consumed by the translated HEAANBOOT circuit.

    packed=True corresponds to the usual packed branch:
        CoeffToSlot: 1
        exp2pi:      3
        squarings:   logT + logI
        EvalExp end: 1
        SlotToCoeff: 1
      total = logT + logI + 6.

    The single-real branch omits the two linear transforms, hence +4 instead.
    """
    """if packed:
        return ZZ(logT) + ZZ(logI) + 6
    return ZZ(logT) + ZZ(logI) + 4"""

    if packed:
        # RNS adaptation consumes one additional prime for:
        #   1. SubSum divByPo2(logNh)
        #   2. EvalExp div/rescaleByPo2(logT+1)
        #   3. final reScaleBy(logI)
        return ZZ(logT) + ZZ(logI) + 9

def configure_heaanboot_chain(logT=0, logI=4, extra_guard_primes=0,
                               keep_existing_eval_primes=True,
                               boot_logp=None, verbose=True):
    """Extend the currently loaded ModHE RNS chain for HEAANBOOT.

    The original q0 and, by default, all existing q1,... are retained as the
    post-bootstrap/evaluation chain.  Large bootstrapping primes are appended at
    the *top* of the chain, so normal ModHE rescaling consumes them first.

    This function recomputes all global constants used by Relin/Flexibility.
    Therefore it MUST be called before key generation.
    """
    global q_list, p, p_enc, L, moduli, Q, p_inv_q, q_hat, q_hat_inv, B

    old_qs = [ZZ(q) for q in list(q_list)]
    if len(old_qs) < 1:
        raise ValueError("q_list must contain q0")

    q0 = old_qs[0]
    if boot_logp is None:
        boot_logp = _hb_floor_log2(q0) + 4
    boot_logp = ZZ(boot_logp)
    # In HEAANBOOT q=2^logq and the boot scale is 2^(logq+4)=16*q.
    # ModHE q0 is an odd RNS prime, so target the exact analogue 16*q0.
    boot_target = ZZ(16) * q0

    # Packed HEAANBOOT is the target here.  By default no guard prime is
    # retained, so after bootstrap the original evaluation chain is exposed again.
    need = int(heaanboot_required_rescale_primes(logT, logI, packed=True))
    nboot = need + int(extra_guard_primes)

    if keep_existing_eval_primes:
        base_qs = old_qs
    else:
        base_qs = [q0]

    used = set(base_qs)
    boot_qs = []
    skip = 0
    while len(boot_qs) < nboot:
        qi = _hb_ntt_prime_below_target(boot_target, used=used, start_skip=skip)
        boot_qs.append(qi)
        used.add(qi)
        # Move the next search farther away so repeated searches are quick.
        skip += 1

    new_qs = base_qs + boot_qs
    q_list = vector(ZZ, new_qs)
    L = len(new_qs)

    # One special prime p, also near the boot scale and distinct from q_i.
    p = _hb_ntt_prime_below_target(boot_target, used=set(new_qs), start_skip=nboot + 7)

    moduli = list(new_qs) + [p]
    Q = ZZ(1)
    for qi in new_qs:
        Q *= ZZ(qi)

    p_inv_q = [ZZ(p).inverse_mod(ZZ(qi)) for qi in new_qs]

    q_hat = [ZZ(0)] * L
    q_hat_inv = [ZZ(0)] * L
    B = [ZZ(0)] * L
    for j in range(L):
        q_hat[j] = Q // ZZ(q_list[j])
        q_hat_inv[j] = q_hat[j].inverse_mod(ZZ(q_list[j]))
        B[j] = q_hat_inv[j] * q_hat[j]

    if verbose:
        print("HEAANBOOT RNS chain configured")
        print("  POLDEG             =", POLDEG)
        print("  log2(q0) floor      =", _hb_floor_log2(q0))
        print("  boot_logp hint      =", boot_logp)
        print("  boot scale target   = 16*q0 =", boot_target)
        print("  existing eval primes=", len(base_qs))
        print("  appended boot primes=", len(boot_qs))
        print("  total L             =", L)
        print("  required rescalings =", need)
        print("  special p bits      =", ZZ(p).nbits())
        print("  WARNING: parameters are for Sage correctness experiments only;")
        print("           concrete MLWE security has not been certified here.")

    return {
        "boot_logp": int(boot_logp),
        "boot_target": boot_target,
        "boot_primes": boot_qs,
        "required_rescales": need,
        "L": L,
        "p": p,
    }
