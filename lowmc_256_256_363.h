#ifndef LOWMC_256_256_363_H
#define LOWMC_256_256_363_H

#include "lowmc_pars.h"

#if !defined(MUL_M4RI)
const lowmc_t* get_lowmc_256_256_363(void);
#else
lowmc_t* get_lowmc_256_256_363(void);
#endif

#endif
