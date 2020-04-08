
def calc_size(digest_size, seed_size, num_rounds, lowmc_k, lowmc_r, lowmc_m, is_unruh):
    lowmc_n = lowmc_k

    # bytes required to store one input share
    input_size = (lowmc_k + 7) >> 3;
    # bytes required to store one output share
    output_size = (lowmc_n + 7) >> 3;
    # number of bits per view per LowMC round
    view_round_size = lowmc_m * 3;
    # bytes required to store communicated bits (i.e. views) of one round
    view_size = (view_round_size * lowmc_r + 7) >> 3;

    # bytes required to store collapsed challenge
    collapsed_challenge_size = (num_rounds + 3) >> 2;

    if is_unruh:
      unruh_without_input_bytes_size = seed_size + view_size;
      unruh_with_input_bytes_size    = unruh_without_input_bytes_size + input_size;
    else:
      unruh_without_input_bytes_size = unruh_with_input_bytes_size = 0;

    # we can use unruh_without_input_bytes_size here. In call cases where we need
    # to write more, we do not need to write the input share
    per_round_size = input_size + view_size + digest_size + 2 * seed_size + unruh_without_input_bytes_size;
    max_signature_size = collapsed_challenge_size + 32 + num_rounds * per_round_size;

    print(unruh_without_input_bytes_size, unruh_with_input_bytes_size, max_signature_size)


# Picnic with partial Sbox layer
calc_size(32, 16, 219, 128, 20, 10, False)
calc_size(32, 16, 219, 128, 20, 10, True)

calc_size(48, 24, 329, 192, 30, 10, False)
calc_size(48, 24, 329, 192, 30, 10, True)

calc_size(64, 32, 438, 256, 38, 10, False)
calc_size(64, 32, 438, 256, 38, 10, True)

# Picnic with full Sbox layer
calc_size(32, 16, 219, 129, 4, 43, False)
calc_size(48, 24, 329, 192, 4, 64, False)
calc_size(64, 32, 438, 255, 4, 85, False)
