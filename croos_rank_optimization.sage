# =============================================================================
# crossrank_c2s_poc.sage
#
# Proof-of-concept extension for the current ModHE/HEAANBOOT Sage prototype.
#
# Goal:
#   q0 rank-1 client ct
#       -> public ModRaise while rank 1
#       -> C2S baby branches: phi_j on rank-1 ct + direct KS phi_j(s_C)->s_E
#       -> ordinary rank-r giant rotations
#       -> ordinary rank-r EvalExp / S2C
#
# IMPORTANT:
#   * Load this file AFTER Flexibility.sage and bootstrap_heaanboot_modhe.sage.
#   * This first PoC keeps HEAANBOOT's existing BSGS baby size k.  That isolates
#     correctness of the fused cross-rank baby layer.  Changing k to the
#     rank-aware optimum requires regenerating the C2S plaintext diagonals.
#   * The current rankup_keygen() generates full-chain P*Q switching material.
#     This is convenient for correctness testing, but over-generates key data;
#     a later optimized implementation should generate only the C2S level keys.
# =============================================================================


class HEAANCrossRankBabyKeys(object):
    def __init__(self, keys, rotation_indices, source_rank, target_rank, k):
        self.keys = keys
        self.rotation_indices = list(rotation_indices)
        self.source_rank = int(source_rank)
        self.target_rank = int(target_rank)
        self.k = int(k)


def _heaan_c2s_default_k(slots):
    slots = int(slots)
    logSlots = ZZ(slots).nbits() - 1
    logk = logSlots // 2
    return 1 << int(logk)


def heaan_crossrank_baby_keygen(key_client, key_eval, slots, verbose=True):
    """Generate C2S baby keys phi_j(s_C) -> s_E for j=0,...,k-1.

    key_client must be rank 1.  key_eval may have arbitrary target rank r.
    The underlying project's rankup_keygen() is used as generalized KSGen.
    Since q_list/Q/p/B have already been extended to the bootstrap chain,
    these keys are generated at the full P*Q chain; rank_up() later consumes
    only the active q-digits at the C2S level.
    """
    if key_client[0].ncols() != 1:
        raise ValueError("heaan_crossrank_baby_keygen expects rank-1 client key")

    r_tgt = int(key_eval[0].ncols())
    k = _heaan_c2s_default_k(slots)
    keys = {}

    if verbose:
        print("Generating", k, "cross-rank C2S baby keys: 1 ->", r_tgt)

    for rot in range(k):
        if rot == 0:
            s_src_rot = _secret_to_list(key_client[0])
        else:
            expn = hb_rotation_exponent(rot)
            s_src_rot = _hb_secret_automorphism(key_client[0], expn)

        # Generalized key switch phi_rot(s_C) -> s_E.
        keys[int(rot)] = rankup_keygen(s_src_rot, key_eval[0])

        if verbose:
            print("  cross baby key rot=", rot,
                  "source=phi_%d(s_C)" % rot,
                  "target rank=", r_tgt)

    return HEAANCrossRankBabyKeys(
        keys=keys,
        rotation_indices=list(range(k)),
        source_rank=1,
        target_rank=r_tgt,
        k=k
    )


def hb_cross_rank_rotate_up(cipher, rot_slots, cross_baby_keys):
    """Compute phi_rot(cipher), then KS phi_rot(s_C) -> s_E.

    Every call starts from the SAME rank-1 input ciphertext supplied by the
    caller.  The baby branches are therefore parallel, not sequential.
    """
    rot_slots = int(rot_slots)

    if cipher.rank() != 1:
        raise ValueError("hb_cross_rank_rotate_up expects a rank-1 input")
    if cross_baby_keys.source_rank != 1:
        raise ValueError("cross baby key source rank must be 1")
    if rot_slots not in cross_baby_keys.keys:
        raise KeyError("missing cross-rank baby key for rotation %d" % rot_slots)

    if rot_slots == 0:
        raw_auto = _hb_raw_copy(cipher.raw)
    else:
        expn = hb_rotation_exponent(rot_slots)
        raw_auto = _hb_raw_automorphism(cipher.raw, expn)

    # raw_auto is under phi_rot(s_C).  Switch directly to s_E.
    raw_out = rank_up(raw_auto, cross_baby_keys.keys[rot_slots])

    out = cipher.copy()
    out.raw = raw_out

    if out.rank() != cross_baby_keys.target_rank:
        raise ValueError("cross-rank rotation produced unexpected target rank")
    return out


def heaan_coeff_to_slot_crossrank(cipher, boot_context,
                                  cross_baby_keys, eval_boot_keys,
                                  verbose=False):
    """C2S with parallel 1->r baby branches and ordinary r->r giants.

    This version intentionally uses the same k and the same pvec_unit[] as the
    current heaan_coeff_to_slot().  It is therefore a correctness PoC for the
    fused rank transition, not yet the rank-aware-k optimization.
    """
    if cipher.rank() != 1:
        raise ValueError("cross-rank C2S expects rank-1 input")

    slots = int(cipher.slots)
    k = _heaan_c2s_default_k(slots)
    if int(cross_baby_keys.k) != k:
        raise ValueError("cross baby key k does not match current C2S context")

    # ------------------------------------------------------------------
    # Baby layer: ALL branches start independently from 'cipher'.
    # Each branch does phi_j at rank 1, followed by 1->r key switching.
    # ------------------------------------------------------------------
    rotvec = [None] * k
    for j in range(k):
        rotvec[j] = hb_cross_rank_rotate_up(cipher, j, cross_baby_keys)

    if verbose:
        print("[cross-C2S] materialized", k,
              "parallel baby ciphertexts at rank", rotvec[0].rank())

    # From here onward all ciphertexts are under the same rank-r s_E.
    tmpvec = [None] * k
    for j in range(k):
        tmpvec[j] = _hb_mult_by_unit_poly(
            rotvec[j], boot_context.pvec_unit[j]
        )
    for j in range(1, k):
        tmpvec[0] = hb_add(tmpvec[0], tmpvec[j])
    out = tmpvec[0]

    # Giant steps remain ordinary same-rank rotations under s_E.
    for ki in range(k, slots, k):
        for j in range(k):
            tmpvec[j] = _hb_mult_by_unit_poly(
                rotvec[j], boot_context.pvec_unit[j + ki]
            )
        for j in range(1, k):
            tmpvec[0] = hb_add(tmpvec[0], tmpvec[j])

        tmp0 = hb_left_rotate_fast(tmpvec[0], ki, eval_boot_keys)
        out = hb_add(out, tmp0)

    return hb_rescale_one(out)


def heaan_bootstrap_crossrank_c2s(cipher, boot_context,
                                  cross_baby_keys, eval_boot_keys,
                                  relin_key, logT=0, logI=4,
                                  output_scale=None, active_raise=None,
                                  verbose=True):
    """Fully-packed bootstrap with ModRaise-before-rank-transition.

    Pipeline:
      rank-1/q0
        -> public ModRaise at rank 1
        -> full-packed public SubSum normalization
        -> C2S baby layer: phi_j(s_C) -> s_E, 1->r
        -> ordinary rank-r EvalExp
        -> ordinary rank-r S2C

    This PoC deliberately rejects sparse packing because sparse SubSum performs
    rotations before C2S and would need its own cross-rank treatment.
    """
    if not isinstance(cipher, ModHEBootCipher):
        raise TypeError("expected ModHEBootCipher")
    if cipher.rank() != 1:
        raise ValueError("optimized bootstrap expects rank-1 client ciphertext")
    if int(cipher.slots) != int(POLDEG // 2):
        raise NotImplementedError(
            "PoC supports fully packed ciphertexts only (slots = N/2)"
        )

    if output_scale is None:
        output_scale = cipher.scale
    if active_raise is None:
        active_raise = L

    # Same q0 normalization / public centered lift as current bootstrap.
    q0raw = _hb_raw_mod_down_to_level(cipher.raw, 0)
    raised = heaan_mod_raise_q0(q0raw, active=active_raise)

    boot_scale = HB_RR(16 * ZZ(q_list[0]))
    out = ModHEBootCipher(
        raised, boot_scale, cipher.slots, cipher.is_complex
    )

    if verbose:
        print("[cross-bootstrap] after rank-1 ModRaise:", out)

    # In the current implementation, for full packing the SubSum rotation loop
    # is empty and the remaining operation is the public 2^{-log(N/2)} scaling.
    # Call it directly to make explicit that no rank-1 high-Q rotation key is used.
    logNh = ZZ(POLDEG // 2).nbits() - 1
    out = hb_div_by_po2(out, int(logNh))

    if verbose:
        print("[cross-bootstrap] after public full-pack normalization:", out)

    # First secret-dependent high-modulus stage: C2S cross-rank baby layer.
    out = heaan_coeff_to_slot_crossrank(
        out,
        boot_context,
        cross_baby_keys,
        eval_boot_keys,
        verbose=verbose
    )

    if verbose:
        print("[cross-bootstrap] after cross-rank CoeffToSlot:", out)

    # From this point on the ciphertext is a standard rank-r MLWE ciphertext.
    out = heaan_eval_exp(
        out, boot_context, eval_boot_keys, relin_key, logT, logI
    )
    if verbose:
        print("[cross-bootstrap] after EvalExp:", out)

    out = heaan_slot_to_coeff(out, boot_context, eval_boot_keys)
    out.scale = HB_RR(output_scale)

    if verbose:
        print("[cross-bootstrap] after SlotToCoeff/final scale reset:", out)
    return out


# =============================================================================
# Correctness diagnostics
# =============================================================================

def _poc_centered_poly_diff_q0(a, b):
    diff = poly_mod_prime((a - b) % firreducible, q_list[0])
    return center_lift_poly_mod_q(diff, q_list[0])


def _poc_max_abs_coeff(poly):
    coeffs = _dense_coeffs(poly, POLDEG)
    if len(coeffs) == 0:
        return ZZ(0)
    return max([abs(ZZ(c)) for c in coeffs])


def test_cross_rank_rotation_primitive(ct_client, key_client, key_eval,
                                       cross_baby_keys,
                                       rotations=None, active_raise=None):
    """Check Dec_sE(XRotUp_j(Lift(ct))) ~= phi_j(Dec_sC(ct)) mod q0.

    This tests the cross-rank rotation primitive independently of C2S/EvalExp.
    The comparison is performed on centered q0 decryption polynomials, so it
    does not depend on CKKS scale metadata.
    """
    if active_raise is None:
        active_raise = L
    if rotations is None:
        rotations = list(cross_baby_keys.rotation_indices)

    q0raw = _hb_raw_mod_down_to_level(ct_client.raw, 0)
    raised_raw = heaan_mod_raise_q0(q0raw, active=active_raise)
    raised = ModHEBootCipher(
        raised_raw, ct_client.scale, ct_client.slots, ct_client.is_complex
    )

    src_dec_q0 = decrypt_poly_q0(ct_client.raw, key_client)

    print("\n================================================")
    print("CROSS-RANK ROTATION PRIMITIVE TEST")
    print("================================================")

    results = {}
    for rot in rotations:
        rot = int(rot)
        got_ct = hb_cross_rank_rotate_up(raised, rot, cross_baby_keys)
        got_dec_q0 = decrypt_poly_q0(got_ct.raw, key_eval)

        if rot == 0:
            expected_q0 = src_dec_q0
        else:
            expn = hb_rotation_exponent(rot)
            expected_q0 = center_lift_poly_mod_q(
                poly_mod_prime(
                    hb_automorphism_poly(src_dec_q0, expn), q_list[0]
                ),
                q_list[0]
            )

        diff = _poc_centered_poly_diff_q0(got_dec_q0, expected_q0)
        maxdiff = _poc_max_abs_coeff(diff)
        results[rot] = maxdiff

        print("rotation", rot,
              ": output rank=", got_ct.rank(),
              " max centered q0 coeff error=", maxdiff)

    return results


def report_crossrank_c2s_work(slots, target_rank, active_c2s=None):
    """Structural ring-product counts for the current generalized KS loops.

    Counts the products with switching-key b/a polynomials performed by
    rank_up(), suppressing NTT/modulus widths.  This matches the current Sage
    loop structure and is useful as a sanity check, not a runtime model.
    """
    slots = int(slots)
    r = int(target_rank)
    k = _heaan_c2s_default_k(slots)
    m = int(slots // k)

    if active_c2s is None:
        # Current full-packed bootstrap consumes one q-prime in SubSum before C2S.
        active_c2s = L - 1
    active_c2s = int(active_c2s)

    same_rank_one_rot = active_c2s * r * (r + 1)
    cross_rank_one_rot = active_c2s * (r + 1)

    base_baby = (k - 1) * same_rank_one_rot
    opt_baby = k * cross_rank_one_rot
    giants = (m - 1) * same_rank_one_rot

    # Existing q0-only standalone RankUp consumes only one q digit.
    baseline_q0_rankup = (r + 1)

    print("\n================================================")
    print("C2S STRUCTURAL KS WORK REPORT")
    print("================================================")
    print("slots S                       =", slots)
    print("target rank r                 =", r)
    print("current HEAAN BSGS k          =", k)
    print("giant groups m=S/k            =", m)
    print("active q digits at C2S        =", active_c2s)
    print("one r->r rotation products    =", same_rank_one_rot)
    print("one 1->r rotation products    =", cross_rank_one_rot)
    print("baseline C2S baby products    =", base_baby)
    print("optimized C2S baby products   =", opt_baby)
    print("giant products (same both)    =", giants)
    print("baseline q0 RankUp products   =", baseline_q0_rankup)
    print("baseline total (RU+C2S KS)    =",
          baseline_q0_rankup + base_baby + giants)
    print("optimized C2S KS total        =", opt_baby + giants)

    if opt_baby + giants > 0:
        print("structural ratio baseline/opt =",
              RR(baseline_q0_rankup + base_baby + giants) /
              RR(opt_baby + giants))

    return {
        "k": k,
        "m": m,
        "active_c2s": active_c2s,
        "baseline_q0_rankup": baseline_q0_rankup,
        "baseline_baby": base_baby,
        "optimized_baby": opt_baby,
        "giants": giants,
    }
