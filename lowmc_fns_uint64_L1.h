/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include "lowmc_fns_undef.h"

#define ADDMUL mzd_addmul_v_uint64_128
#define MUL mzd_mul_v_uint64_128
#define XOR mzd_xor_uint64_128
#define COPY mzd_copy_uint64_128
#define MPC_MUL mpc_matrix_mul_uint64_128

#if defined(WITH_LOWMC_126_126_4)
#define LOWMC_INSTANCE lowmc_126_126_4
#endif
#define LOWMC_N LOWMC_L1_N
#define LOWMC_R LOWMC_L1_R
#define LOWMC_M LOWMC_L1_M
