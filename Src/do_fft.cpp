#include <stddef.h>
#include "../Include/do_fft.h"
#include "../Include/NE10_fft.h"

void Do_fftr(float* data_out, float* data_in, const int fft_len, TFFTFormat format)
{
	if (NULL == data_out || NULL == data_in){
		return;
	}
	ne10_fft_r2c_cfg_float32_t cfg = ne10_fft_alloc_r2c_float32(fft_len);
	ne10_fft_cpx_float32_t cx_out[4096 / 2 + 1]; //max 4096
	ne10_fft_r2c_1d_float32_c(cx_out, data_in, cfg);
	int idx = 0;

	switch (format)
	{
	case kHalfComplexInPlace:
		data_out[0] = cx_out[0].r;
		data_out[fft_len / 2] = cx_out[fft_len / 2].r;
		for (idx = 1; idx < fft_len / 2; idx++)
		{
			data_out[idx] = cx_out[idx].r;
			data_out[fft_len / 2 - idx] = cx_out[idx].i;
		}
		break;
	case kIntelPerm:
		data_out[0] = cx_out[0].r;
		data_out[1] = cx_out[fft_len / 2].r;
		data_out = data_out + 2;
		for (idx = 1; idx < fft_len / 2; idx++)
		{
			*data_out++ = cx_out[idx].r;
			*data_out++ = cx_out[idx].i;
		}
		break;
	case kIntelCCS:
		for (idx = 0; idx < fft_len / 2 + 1; idx++)
		{
			*data_out++ = cx_out[idx].r;
			*data_out++ = cx_out[idx].i;
		}
		break;
	default:
		break;
	}
	ne10_fft_destory_r2c_float32(cfg);
}

void Do_ifftr(float* data_out, float* data_in, const int fft_len, TFFTFormat format)
{
	if (NULL == data_out || NULL == data_in) {
		return;
	}
	ne10_fft_r2c_cfg_float32_t cfg = ne10_fft_alloc_r2c_float32(fft_len);
	ne10_fft_cpx_float32_t cx_in[4096 / 2 + 1]; //max 4096

	float* data_in_p;
	int idx = 0, i = 0;
	switch (format)
	{
	case kHalfComplexInPlace:
		data_in_p = data_in;
		cx_in[0].r = data_in_p[0];
		cx_in[0].i = 0;
		cx_in[fft_len / 2].r = data_in_p[fft_len / 2];
		cx_in[fft_len / 2].i = 0;

		for (idx = 1; idx < fft_len / 2; idx++)
		{
			cx_in[idx].r = data_in_p[idx];
			cx_in[idx].i = data_in_p[fft_len / 2 - idx];
		}
		break;
	case kIntelPerm:
		cx_in[0].r = data_in[0];
		cx_in[0].i = 0;
		cx_in[fft_len / 2].r = data_in[1];
		cx_in[fft_len / 2].i = 0;
		data_in_p = data_in + 2;
		for (idx = 1; idx < fft_len / 2; idx++)
		{
			cx_in[idx].r = *data_in_p++;
			cx_in[idx].i = *data_in_p++;
		}
		break;
	case kIntelCCS:
		data_in_p = data_in;
		for (idx = 0; idx < fft_len / 2 + 1; idx++)
		{
			cx_in[idx].r = *data_in_p++;
			cx_in[idx].i = *data_in_p++;
		}
		cx_in[0].i = 0;
		cx_in[fft_len / 2].i = 0;
		break;
	default:
		break;
	}

	if (cfg != NULL)
	{
		ne10_fft_c2r_1d_float32_c(data_out,cx_in,cfg);
	}
	ne10_fft_destory_r2c_float32(cfg);
}