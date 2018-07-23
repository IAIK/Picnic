#ifndef LOWMC_192_192_30_H
#define LOWMC_192_192_30_H

#include "lowmc_pars.h"

#if !defined(MUL_M4RI)
const lowmc_t* get_lowmc_192_192_30(void);
#else
lowmc_t* get_lowmc_192_192_30(void);
#endif

#endif
