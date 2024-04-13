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

#include "beat.hpp"
#include "xtract/libxtract.h"

using namespace li;

namespace li
{
    struct onset_func
    {
        ivec bands = ivec(XTRACT_BARK_BANDS);
        fvec bark = fvec(XTRACT_BARK_BANDS);
        fvec last = fvec(XTRACT_BARK_BANDS);
        size_t nfft;

        onset_func(int nfft_, float fs) : nfft(nfft_)
        {
            xtract_init_bark(nfft, fs, bands.data());
            std::fill(bark.begin(), bark.end(), __FLT_EPSILON__);
            std::fill(last.begin(), last.end(), __FLT_EPSILON__);
        }

        float process(const float* spec)
        {
            xtract_bark_coefficients(spec, nfft / 2, bands.data(), bark.data());
            float of = sub(bark, last);
            last = bark;
            return of;
        }

        fvec get_bark() const
        {
            return bark;
        }
    };
}

beat_det::beat_det(float fs_) : fs(fs_)
{
    onset = new onset_func(fft_size, fs);
    spectrum::gen_hann_window(window.data.data(), window_size);
    assert(onset);
}

beat_det::~beat_det()
{
    delete onset;
}

int beat_det::process(const fmat &in)
{
    assert(in.cols == frame_size);
    input.roll(in);
    input.weight(window);
    spec.calc_mag(input.row(0), temp.row(0), fs);
    float of = onset->process(temp.row(0));
    return process(of);
}

void beat_det::calc_spec(const fvec &in, fvec &out)
{
    assert(in.size() == fft_size);
    assert(out.size() == fft_size / 2);
}

float beat_det::calc_flux(const fvec &spec)
{
    return 0;
}

int beat_det::process(float flux)
{
    flux += __FLT_EPSILON__;

    history.roll(flux);

    update_score(flux);
    
    next_beat--;
    if (next_beat <= 0)
    {
        next_beat = update_prd();
    }

    return next_beat;
}

float beat_det::update_prd()
{
    current_prd = 0.0f;
    return 0.0f;
}

int beat_det::predict_next_beat()
{
    return 0;
}

void beat_det::update_score(float flux)
{
    score.roll(flux);
}
