#include <xtract/xtract_types.h>

/* FFT functions */
void cdft(int n, int isgn, real_t* a, int* ip, real_t* w);
void rdft(int n, int isgn, real_t *a, int *ip, real_t *w);
void ddct(int n, int isgn, real_t *a, int *ip, real_t *w);
void ddst(int n, int isgn, real_t *a, int *ip, real_t *w);
void dfct(int n, real_t *a, real_t *t, int *ip, real_t *w);
void dfst(int n, real_t *a, real_t *t, int *ip, real_t *w);

/* Auxiliary functions */
void makewt(int nw, int *ip, real_t *w);
void bitrv2(int n, int *ip, real_t *a);
void bitrv2conj(int n, int *ip, real_t *a);
void cftfsub(int n, real_t *a, real_t *w);
void cftbsub(int n, real_t *a, real_t *w);
void makect(int nc, int *ip, real_t *c);
void rftfsub(int n, real_t *a, int nc, real_t *c);
void rftbsub(int n, real_t *a, int nc, real_t *c);
void dctsub(int n, real_t *a, int nc, real_t *c);
void dstsub(int n, real_t *a, int nc, real_t *c);
void cft1st(int n, real_t *a, real_t *w);
void cftmdl(int n, int l, real_t *a, real_t *w);
