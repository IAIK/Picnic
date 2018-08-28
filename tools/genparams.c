#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <m4ri/m4ri.h>
#include <openssl/rand.h>

typedef struct {
  mzd_t* k_matrix;
  mzd_t* l_matrix;
  mzd_t* constant;
} lowmc_round_t;

/**
 * Represents the LowMC parameters as in https://bitbucket.org/malb/lowmc-helib/src,
 * with the difference that key in a separate struct
 */
typedef struct {
  size_t m;
  size_t n;
  size_t r;
  size_t k;

  mzd_t* k0_matrix;
  lowmc_round_t* rounds;

  mzd_t* precomputed_non_linear_part_matrix;
  mzd_t* precomputed_constant_linear;
  mzd_t* precomputed_constant_non_linear;
} lowmc_t;

static void mzd_randomize_ssl(mzd_t* val) {
  for (rci_t i = 0; i < val->nrows; ++i) {
    RAND_bytes((unsigned char*)val->rows[i], val->width * sizeof(word));
  }
}

static mzd_t* mzd_sample_matrix_word(uint32_t n, uint32_t k, int rank, bool with_xor) {
  // use mzd_init for A since m4ri will work with it in mzd_echolonize
  // also, this function cannot be parallelized as mzd_echolonize will call
  // mzd_init and mzd_free at will causing various crashes.
  mzd_t* A = mzd_init(n, k);
  mzd_t* B = mzd_init(n, k);
  do {
    mzd_randomize_ssl(A);
    if (with_xor) {
      for (uint32_t i = 0; i < n; i++) {
        mzd_xor_bits(A, n - i - 1, (k + i + 1) % k, 1, 1);
      }
    }
    mzd_copy(B, A);
  } while (mzd_echelonize(A, 0) != rank);
  mzd_free(A);
  return B;
};

/**
 * Samples the L matrix for the LowMC instance
 *
 * \param n the blocksize
 */
static mzd_t* mzd_sample_lmatrix(uint32_t n) {
  return mzd_sample_matrix_word(n, n, n, false);
}

/**
 * Samples the K matrix for the LowMC instance
 * \param n the blocksize
 */
static mzd_t* mzd_sample_kmatrix(uint32_t n, uint32_t k) {
  return mzd_sample_matrix_word(n, k, MIN(n, k), true);
}

#if defined(REDUCED_LINEAR_LAYER)
static bool precompute_values_for_lowmc(lowmc_t* lowmc) {
  mzd_t* L_inverse[lowmc->r];
  mzd_t* Li_K[lowmc->r];
  mzd_t* Li_C[lowmc->r];

  mzd_t* identity = mzd_init(lowmc->n, lowmc->k);
  mzd_set_ui(identity, 1);

  lowmc->precomputed_constant_non_linear = mzd_init(1, (3 * lowmc->m + 2) * lowmc->r);
  lowmc->precomputed_non_linear_part_matrix = mzd_init(lowmc->n, (3 * lowmc->m + 2) * lowmc->r);
  const unsigned int bound = lowmc->n - 3 * lowmc->m;

  for (unsigned round_i = 0; round_i < lowmc->r; ++round_i) {
    L_inverse[round_i] = mzd_invert_naive(NULL, lowmc->rounds[round_i].l_matrix, identity);
    Li_K[round_i]      = mzd_mul_naive(NULL, lowmc->rounds[round_i].k_matrix, L_inverse[round_i]);
    Li_C[round_i]      = mzd_mul_naive(NULL, lowmc->rounds[round_i].constant, L_inverse[round_i]);

    for (unsigned int row = bound; row < lowmc->n; row++) {
      mzd_row_clear_offset(L_inverse[round_i], row, 0);
    }
  }

  mzd_free(identity);

  for (unsigned int round_i = 0; round_i < lowmc->r; ++round_i) {
    mzd_t* tmp = mzd_copy(NULL, Li_K[round_i]);
    mzd_t* tmpC = mzd_copy(NULL, Li_C[round_i]);
    for (unsigned int i = round_i + 1; i < lowmc->r; ++i) {
      mzd_t* x = mzd_copy(NULL, Li_K[i]);
      mzd_t* c = mzd_copy(NULL, Li_C[i]);

      for (unsigned int j = i; j > round_i; --j) {
        mzd_t* t = mzd_mul_naive(NULL, x, L_inverse[j - 1]);
        mzd_free(x);
        x = t;

        t = mzd_mul_naive(NULL, c, L_inverse[j - 1]);
        mzd_free(c);
        c = t;
      }
      mzd_add(tmp, tmp, x);
      mzd_add(tmpC, tmpC, c);

      mzd_free(x);
      mzd_free(c);
    }

    const unsigned int idx = round_i * (3 * lowmc->m + 2);
    for (unsigned int c = 0; c < 3 * lowmc->m; ++c) {
      mzd_write_bit(lowmc->precomputed_constant_non_linear, 0, idx + 2 + c, mzd_read_bit(tmpC, 0, bound + c));
      for (unsigned int row = 0; row < lowmc->n; ++row) {
        mzd_write_bit(lowmc->precomputed_non_linear_part_matrix, row, idx + 2 + c, mzd_read_bit(tmp, row, bound + c));
      }
    }

    if (!round_i) {
      for (unsigned int row = 0; row < lowmc->n; ++row) {
        mzd_row_clear_offset(tmp, row, bound);
      }
      mzd_add(lowmc->k0_matrix, lowmc->k0_matrix, tmp);

      mzd_row_clear_offset(tmpC, 0, bound);
      lowmc->precomputed_constant_linear = tmpC;
    } else {
      mzd_free(tmpC);
    }
    mzd_free(tmp);

    // mzd_local_copy(round->precomputed_non_linear_part_matrix, tmp);
  }

  for (unsigned round_i = 0; round_i < lowmc->r; ++round_i) {
    mzd_free(L_inverse[round_i]);
    mzd_free(Li_K[round_i]);
    mzd_free(Li_C[round_i]);
  }


  return true;
}
#endif

static void writeMZD_TStructToFile(mzd_t* matrix, FILE* file) {
  fwrite(&matrix->nrows, sizeof(uint32_t), 1, file);
  fwrite(&matrix->ncols, sizeof(uint32_t), 1, file);

  for (int i = 0; i < matrix->nrows; i++) {
    fwrite(matrix->rows[i], matrix->rowstride * sizeof(word), 1, file);
  }
}

static bool lowmc_write_file(lowmc_t* lowmc, const char* file_name) {
  FILE* file = fopen(file_name, "w");
  if (!file) {
    return false;
  }

  uint32_t supported_type = 0;
#if defined(REDUCED_LINEAR_LAYER)
  supported_type = 0x01;
#endif

  fwrite(&supported_type, sizeof(supported_type), 1, file);
  fwrite(&lowmc->n, sizeof(lowmc->n), 1, file);
  fwrite(&lowmc->k, sizeof(lowmc->k), 1, file);
  fwrite(&lowmc->m, sizeof(lowmc->m), 1, file);
  fwrite(&lowmc->r, sizeof(lowmc->r), 1, file);

  writeMZD_TStructToFile(lowmc->k0_matrix, file);
  for (size_t i = 0; i < lowmc->r; ++i) {
    writeMZD_TStructToFile(lowmc->rounds[i].l_matrix, file);
#if !defined(REDUCED_LINEAR_LAYER)
    writeMZD_TStructToFile(lowmc->rounds[i].k_matrix, file);
    writeMZD_TStructToFile(lowmc->rounds[i].constant, file);
#endif
  }
#if defined(REDUCED_LINEAR_LAYER)
  writeMZD_TStructToFile(lowmc->precomputed_non_linear_part_matrix, file);
  writeMZD_TStructToFile(lowmc->precomputed_constant_linear, file);
  writeMZD_TStructToFile(lowmc->precomputed_constant_non_linear, file);
#endif

  fclose(file);
  return true;
}

static bool lowmc_generate(lowmc_t* lowmc, size_t m, size_t n, size_t r, size_t k) {
  if (!lowmc) {
    return false;
  }

  if (n - 3 * m < 2) {
    return false;
  }

  lowmc->m = m;
  lowmc->n = n;
  lowmc->r = r;
  lowmc->k = k;

  lowmc->rounds    = calloc(sizeof(lowmc_round_t), r);
  lowmc->k0_matrix = mzd_sample_kmatrix(k, n);
  for (unsigned int i = 0; i < r; ++i) {
    lowmc->rounds[i].l_matrix = mzd_sample_lmatrix(n);
    lowmc->rounds[i].k_matrix = mzd_sample_kmatrix(k, n);
    lowmc->rounds[i].constant = mzd_init(1, n);
    mzd_randomize_ssl(lowmc->rounds[i].constant);
  }
#ifdef REDUCED_LINEAR_LAYER
  precompute_values_for_lowmc(lowmc);
#endif

  return true;
}

static void lowmc_clear(lowmc_t* lowmc) {
  for (unsigned i = 0; i < lowmc->r; ++i) {
    mzd_free(lowmc->rounds[i].constant);
    mzd_free(lowmc->rounds[i].k_matrix);
    mzd_free(lowmc->rounds[i].l_matrix);
  }
  mzd_free(lowmc->precomputed_non_linear_part_matrix);
  mzd_free(lowmc->precomputed_constant_linear);
  mzd_free(lowmc->precomputed_constant_non_linear);

  mzd_free(lowmc->k0_matrix);
  free(lowmc->rounds);
}

static bool parse_arg(long* value, const char* arg) {
  errno        = 0;
  const long v = strtol(arg, NULL, 10);

  if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN)) || (errno != 0 && v == 0)) {
    return false;
  }
  *value = v;

  return true;
}

static int parse_args(long params[4], int argc, char** argv) {
  if (argc != 6) {
    printf("usage: %s [Number of SBoxes] [Blocksize] [Rounds] [Keysize] [File]\n", argv[0]);
    return -1;
  }

  for (size_t s = 0; s < 4; ++s) {
    if (!parse_arg(&params[s], argv[s + 1])) {
      printf("Unable to parse '%s' as base 10 number.\n", argv[s + 1]);
      return -1;
    }

    if (params[s] <= 0) {
      printf("Expexted positive number, got: %li.\n", params[s]);
      return -1;
    }
  }

  if (params[0] * 3 > params[1]) {
    printf("Number of S-boxes * 3 exceeds block size!");
    return -1;
  }

  return 0;
}

int main(int argc, char** argv) {
  long args[4];
  int ret = parse_args(args, argc, argv);

  if (!ret) {
    lowmc_t lowmc;
    if (!lowmc_generate(&lowmc, args[0], args[1], args[2], args[3])) {
      printf("Failed to generate LowMC instance.\n");
      ret = 1;
    } else {
      if (!lowmc_write_file(&lowmc, argv[5])) {
        printf("Failed to write LowMC instance.\n");
        ret = 1;
      }
      lowmc_clear(&lowmc);
    }
  }

  return ret;
}
