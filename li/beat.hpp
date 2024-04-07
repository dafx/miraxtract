/*
 * Copyright (c) 2024, Yan Li
 * All rights reserved.
 */

#pragma once

#include "common.hpp"

namespace li
{
    struct onset_func;

    struct beat_det
    {

        constexpr static size_t history_size = 43;
        constexpr static size_t fft_size = 512;
        constexpr static size_t window_size = fft_size;
        constexpr static size_t frame_size = fft_size / 4;
        constexpr static size_t downsample = 4;

        beat_det(float fs);
        ~beat_det();

        /**
         * @brief Process a frame of audio data.
         * @return distance to the next beat.
         */        
        int process(const fmat &in);
        int process(float onsetfunc);

        void calc_spec(const fvec &in, fvec &out);
        float calc_onsetfunc(const fvec &spec);

        fmat window = {1, window_size};
        fmat input = {1, window_size};
        fmat temp = {1, window_size};
        fmat history = {1, history_size};

        onset_func *onset = nullptr;
        spectrum spec = {fft_size};
        float fs;
    };
}
