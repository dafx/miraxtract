/*******************************************************
 * Copyright 2024 Li Zhou
 *
 * This software is the property of Li Zhou.
 * Unauthorized copying of this software, or any portion
 * of it, via any medium, is strictly prohibited. This
 * includes, but is not limited to, modification or
 * creation of derivative works based on the original
 * software.
 *
 * Any use or distribution of this software requires
 * explicit written permission from Li Zhou.
 *
 * Li Zhou reserves all rights not expressly
 * granted to the user under this agreement.
 *******************************************************/

#include "common.hpp"
#include "xtract/libxtract.h"

using namespace li;

spectrum::spectrum(size_t fft_size_) : fft_size(fft_size_)
{
    assert((fft_size & (fft_size - 1)) == 0);
    xtract_init_fft(fft_size, XTRACT_SPECTRUM);
}

void spectrum::calc_mag(const float* sig, float* spec, float fs)
{
    real_t argd[4];
    argd[0] = fs / fft_size;
    argd[1] = XTRACT_MAGNITUDE_SPECTRUM;
    argd[2] = 0.0;
    xtract_spectrum(sig, fft_size, argd, spec);
}

void spectrum::gen_hann_window(float *window, size_t sz)
{
    for (size_t i = 0; i < sz; ++i)
    {
        window[i] = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (sz - 1)));
    }
}
