#ifndef __DO_FFT_H__
#define __DO_FFT_H__
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#include "./NE10_fft.h"
typedef enum _TFFTFormat
{
    kHalfComplexInPlace = 0, // r[0],r[1],...,r[n/2],i[n/2 - 1]...i[1]
    kIntelPerm,// r[0],r[n/2],r[1],i[1],...,r[n/2 - 1],i[n/2 - 1]
    kIntelCCS // r[0],0,r[1],i[1],...,r[n/2 - 1],i[n/2 - 1],r[n/2],0
}TFFTFormat;

void Do_fftr(float* data_out, float* data_in, const int fft_len, TFFTFormat format);
void Do_ifftr(float* data_out, float* data_in, const int fft_len, TFFTFormat format);
#ifdef __cplusplus
}
#endif
#endif