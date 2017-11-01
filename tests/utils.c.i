
#include <openssl/rand.h>

void mzd_randomize_ssl(mzd_local_t* val) {
  for (unsigned int i = 0; i < val->nrows; ++i) {
    RAND_bytes((unsigned char*) ROW(val, i), val->width * sizeof(word));
  }
}

mzd_local_t* mzd_convert(mzd_t* v) {
  mzd_local_t* r = mzd_local_init(v->nrows, v->ncols);

  for (rci_t i = 0; i < v->nrows; ++i) {
    memcpy(ROW(r, i), v->rows[i], v->width * sizeof(word));
  }

  return r;
}
