#ifndef _AVX_BLAS_
#define _AVX_BLAS_

extern "C" {
    void y_e_ax(long n, float *a, float *x, float *y) __attribute__ ((noinline));
    void y_pe_ax(long n, float *a, float *x, float *y) __attribute__ ((noinline));
    void z_e_axy(long n, float *a, float *x, float *y, float* z) __attribute__ ((noinline));
    void z_pe_axy(long n, float *a, float *x, float *y, float* z) __attribute__ ((noinline));
}

#endif



