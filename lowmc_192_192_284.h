#ifndef LOWMC_192_192_284_H
#define LOWMC_192_192_284_H

#include "lowmc_pars.h"

#if !defined(MUL_M4RI)
const lowmc_t* get_lowmc_192_192_284(void);
#else
lowmc_t* get_lowmc_192_192_284(void);
#endif

#endif
