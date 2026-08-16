# -----------------------------------------------------------------------------
# boost_rankup_bootstrap_example.sage
#
# Combined client-Boost path for this project:
#   rank-1 q0-only client encryption
#       -> MLWE rank-up/key-switch to evaluation rank
#       -> HEAANBOOT-style ModHE bootstrap
#       -> refreshed evaluation ciphertext.
# -----------------------------------------------------------------------------

load("constants.sage")
load("bootstrap_chain.sage")

LOGT = 0
LOGI = 4
configure_heaanboot_chain(logT=LOGT, logI=LOGI,
                          extra_guard_primes=0,
                          keep_existing_eval_primes=True,
                          verbose=True)

load("auxiliary_functions.sage")
load("KeyGen_rank_up_down.sage")
load("Encrypt-Decrypt_rank_up_down.sage")
load("HMult.sage")
load("Relin_rank_up_down.sage")
load("Flexibility.sage")
load("encode.sage")
load("encode-encrypt.sage")
load("bootstrap_heaanboot_modhe.sage")

R_CLIENT = 1
R_EVAL = MLWE
SLOTS = POLDEG // 2

key_client = public_private_keygen(rank=R_CLIENT)
key_eval = public_private_keygen(rank=R_EVAL)

# All server-side bootstrap/multiplication keys are under the evaluation secret.
relin_key = relin_keygen(key_eval)
boot_context = heaan_add_boot_context(SLOTS, verbose=True)
boot_keys = heaan_boot_keygen(key_eval, SLOTS, verbose=True)

# Key-switch from the client MLWE secret to the evaluation MLWE secret.
rankup_key = rankup_keygen(key_client[0], key_eval[0])

base_vals = [
    0.035 + 0.010*I,
   -0.020 + 0.015*I,
    0.010 - 0.025*I,
   -0.030 - 0.005*I,
]
vals = [base_vals[i % len(base_vals)] for i in range(SLOTS)]

# Client: one RNS base, rank 1.
ct_client_raw = client_encrypt_q0_ckks(vals, key_client, scale=p_enc)
ct_client = hb_wrap(ct_client_raw, scale=p_enc, slots=SLOTS, is_complex=True)
print("\nclient ciphertext:", ct_client)
print("client decrypt:", hb_decrypt_slots(ct_client, key_client))

# Server: rank-up at q0.  This is a genuine same-message key switch, not zero padding.
ct_eval_q0_raw = rank_up(ct_client.raw, rankup_key)
ct_eval_q0 = hb_wrap(ct_eval_q0_raw, scale=p_enc, slots=SLOTS, is_complex=True)
print("\nafter rank-up:", ct_eval_q0)
print("rank-up decrypt:", hb_decrypt_slots(ct_eval_q0, key_eval))

# Server: full modulus refresh/bootstrapping under the evaluation key.
ct_boot = heaan_bootstrap(ct_eval_q0, boot_context, boot_keys, relin_key,
                          logT=LOGT, logI=LOGI,
                          output_scale=p_enc,
                          active_raise=L,
                          verbose=True)
print("\nafter bootstrap:", ct_boot)
print("bootstrap decrypt:", hb_decrypt_slots(ct_boot, key_eval))
