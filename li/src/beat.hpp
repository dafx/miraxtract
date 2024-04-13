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

#pragma once

#include "common.hpp"

namespace li
{
    struct onset_func;

    struct beat_det
    {

        constexpr static size_t history_size = 512;         // ~ 11 s
        constexpr static size_t fft_size = 512;             // ~ 80 ms
        constexpr static size_t window_size = fft_size;
        constexpr static size_t frame_size = fft_size / 4;  // ~ 20 ms
        constexpr static size_t downsample = 4;             // 6000 hz

        beat_det(float fs);
        ~beat_det();

        /**
         * @brief Process a frame of audio data.
         * @return distance to the next beat.
         */        
        int process(const fmat &in);
        int process(float flux);

        void calc_spec(const fvec &in, fvec &out);
        float calc_flux(const fvec &spec);

    private:

        const range bpm_range = {70.0f, 160.0f};
        
        float update_prd();
        int predict_next_beat();
        void update_score(float flux);

        fmat window = {1, window_size};
        fmat input = {1, window_size};
        fmat temp = {1, window_size};
        fmat history = {1, history_size};
        fmat score = {1, history_size};

        onset_func *onset = nullptr;
        spectrum spec = {fft_size};
        float fs;
        float current_prd = 0.0f;
        int next_beat = 0;
        range prd_range = {0.0f, 0.0f};
    };
}
