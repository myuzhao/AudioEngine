/*
 *  Copyright 2013-16 ARM Limited and Contributors.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the name of ARM Limited nor the
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY ARM LIMITED AND CONTRIBUTORS "AS IS" AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL ARM LIMITED AND CONTRIBUTORS BE LIABLE FOR ANY
 *  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * NE10 Library : dsp/NE10_fft.h
 */
#include <stdint.h>

#ifndef NE10_FFT_H
#define NE10_FFT_H
#ifdef __cplusplus
extern "C" {
#endif
#define NE10_UNROLL_LEVEL 0
/////////////////////////////////////////////////////////
// constant values that are used across the library/////
// ////////////////////////////////////////////////////
#define NE10_OK 0
#define NE10_ERR -1

/////////////////////////////////////////////////////////
// some external definitions to be exposed to the users
/////////////////////////////////////////////////////////

    typedef int8_t   ne10_int8_t;
    typedef uint8_t  ne10_uint8_t;
    typedef int16_t  ne10_int16_t;
    typedef uint16_t ne10_uint16_t;
    typedef int32_t  ne10_int32_t;
    typedef uint32_t ne10_uint32_t;
    typedef int64_t  ne10_int64_t;
    typedef uint64_t ne10_uint64_t;
    typedef float    ne10_float32_t;
    typedef double   ne10_float64_t;
    typedef int      ne10_result_t;     // resulting [error-]code
    ///////////////////////////
// Internal macro define
///////////////////////////
#define NE10_FFT_BYTE_ALIGNMENT 8
#define NE10_INLINE inline static

/*
 * FFT Algorithm Flags
 *
 * These are used within Ne10 to decide, after factoring an FFT into stages, what
 * FFT algorithm should be used.
 *
 * - NE10_FFT_ALG_DEFAULT is a mixed radix 2/4 algorithm.
 * - NE10_FFT_ALG_ANY is designated specifically for non-power-of-two input sizes.
 */
#define NE10_FFT_ALG_DEFAULT  0
#define NE10_FFT_ALG_ANY      1

 /*
  * FFT Factor Flags
  *
  * These are used within Ne10 to decide how an input FFT size should be factored into
  * stages (i.e. what radices should be used).
  *
  * - NE10_FACTOR_DEFAULT factors into 2, 3, 4, 5.
  * - NE10_FACTOR_EIGHT_FIRST_STAGE is NE10_FACTOR_DEFAULT with the extended ability to
  *   have a radix-8 initial stage.
  * - NE10_FACTOR_EIGHT factors into 2, 3, 4, 5, 8.
  */
#define NE10_FACTOR_DEFAULT             0
#define NE10_FACTOR_EIGHT_FIRST_STAGE   1
#define NE10_FACTOR_EIGHT               2

  // Comment when do not want to scale output result
#define NE10_DSP_RFFT_SCALING
#define NE10_DSP_CFFT_SCALING

#define NE10_FFT_PARA_LEVEL 4
/////////////////////////////////////////////////////////
// definitions for fft
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
// constant values that are used across the library
/////////////////////////////////////////////////////////

#define NE10_PI (ne10_float32_t)(3.1415926535897932384626433832795)

/////////////////////////////////////////////////////////
// some external macro definitions to be exposed to the users
/////////////////////////////////////////////////////////

#define NE10_MALLOC malloc
#define NE10_FREE(p) \
    do { \
        free(p); \
        p = 0; \
    }while(0)

#define NE10_MIN(a,b) ((a)>(b)?(b):(a))
#define NE10_MAX(a,b) ((a)<(b)?(b):(a))

#define NE10_BYTE_ALIGNMENT(address, alignment) \
    do { \
        (address) = (((address) + ((alignment) - 1)) & ~ ((alignment) - 1)); \
    }while (0)

/**
 * @brief Structure for the floating point FFT function.
 */
#define NE10_MAXFACTORS             32
    typedef struct
    {
        ne10_float32_t r;
        ne10_float32_t i;
    } ne10_fft_cpx_float32_t;

    /**
     * @brief Structure for the floating point FFT state
     */
    typedef struct
    {
        ne10_int32_t nfft;
        ne10_int32_t* factors;
        ne10_fft_cpx_float32_t* twiddles;
        ne10_fft_cpx_float32_t* buffer;
        ne10_fft_cpx_float32_t* last_twiddles;
        /**
         *  @brief Flag to control scaling behaviour in forward floating point complex FFT.
         *  @note If is_forward_scaled is set 0, Ne10 will not scale output of forward floating
         *  point complex FFT. Otherwise, Ne10 will scale output of forward floating
         *  point complex FFT.
         *  @warning Only non-power-of-two FFTs are affected by this flag.
         */
        ne10_int32_t is_forward_scaled;
        /**
         *  @brief Flag to control scaling behaviour in backward floating point complex FFT.
         *  @note If is_backward_scaled is set 0, Ne10 will not scale output of backward floating
         *  point complex FFT. Otherwise, Ne10 will scale output of backward floating
         *  point complex FFT.
         *  @warning Only non-power-of-two FFTs are affected by this flag.
         */
        ne10_int32_t is_backward_scaled;
    } ne10_fft_state_float32_t;

    /**
     * @brief Configuration structure for floating point FFT.
     */
    typedef ne10_fft_state_float32_t* ne10_fft_cfg_float32_t;

    typedef struct
    {
        ne10_fft_cpx_float32_t* buffer;
#if (NE10_UNROLL_LEVEL == 0)
        ne10_int32_t ncfft;
        ne10_int32_t* factors;
        ne10_fft_cpx_float32_t* twiddles;
        ne10_fft_cpx_float32_t* super_twiddles;
#elif (NE10_UNROLL_LEVEL > 0)
        ne10_int32_t nfft;
        ne10_fft_cpx_float32_t* r_twiddles;
        ne10_int32_t* r_factors;
        ne10_fft_cpx_float32_t* r_twiddles_backward;
        ne10_fft_cpx_float32_t* r_twiddles_neon;
        ne10_fft_cpx_float32_t* r_twiddles_neon_backward;
        ne10_int32_t* r_factors_neon;
        ne10_fft_cpx_float32_t* r_super_twiddles_neon;
#endif
    } ne10_fft_r2c_state_float32_t;

    typedef ne10_fft_r2c_state_float32_t* ne10_fft_r2c_cfg_float32_t;
///////////////////////////
// function prototypes:
///////////////////////////
    ne10_fft_r2c_cfg_float32_t ne10_fft_alloc_r2c_float32(ne10_int32_t nfft);
    void ne10_fft_destory_r2c_float32(ne10_fft_r2c_cfg_float32_t cfg);
    void ne10_fft_r2c_1d_float32_c(ne10_fft_cpx_float32_t* fout,ne10_float32_t* fin,ne10_fft_r2c_cfg_float32_t cfg);
    void ne10_fft_c2r_1d_float32_c(ne10_float32_t* fout,ne10_fft_cpx_float32_t* fin,ne10_fft_r2c_cfg_float32_t cfg);
#ifdef __cplusplus
}
#endif

#endif
