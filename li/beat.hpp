/*
 * Copyright (c) 2024, Yan Li
 * All rights reserved.
 */

#pragma once

#include "common.hpp"

namespace li
{
    struct beat_det {

        constexpr static size_t history_size = 43;
        constexpr static size_t window_size = 512;
        constexpr static size_t fft_size = 512;
        constexpr static size_t frame_size = 128;
        constexpr static size_t downsample = 4;

        beat_det()
        {
        }

        void process(const fvec &in)
        {
            assert(in.size() == frame_size);
        }

        int get_next_beat() {
            return 0;
        }

        fvec window = {window_size};
        fvec input = {fft_size};
        fvec temp = {fft_size};
        fmat history = {1, history_size};
    };
}
