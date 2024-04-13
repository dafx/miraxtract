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

#include <vector>
#include <cmath>
#include <cassert>
#include <cstring>

namespace li
{
    typedef std::vector<float> fvec;
    typedef std::vector<int> ivec;

    inline float sub(fvec &a, fvec &b)
    {
        assert(a.size() == b.size());
        float sum = 0.0f;
        for (size_t i = 0; i < a.size(); ++i)
        {
            sum += a[i] - b[i];
        }
        return sum;
    }

    struct fmat
    {
        fvec data;
        size_t rows, cols;

        fmat(size_t r = 1, size_t c = 1, float v = 0.0f) 
            : data(r * c), rows(r), cols(c) 
        {
            std::fill(data.begin(), data.end(), v);
        }

        void resize(size_t r, size_t c)
        {
            data.resize(r * c);
            rows = r;
            cols = c;
        }

        void fill(float v)
        {
            std::fill(data.begin(), data.end(), v);
        }

        float& operator()(size_t r, size_t c)
        {
            assert(r < rows && c < cols);
            return data[r * cols + c];
        }

        // quick access for first row
        float& operator[](size_t c)
        {
            assert(c < cols);
            return data[c];
        }

        float* row(size_t r)
        {
            assert(r < rows);
            return data.data() + r * cols;
        }

        void copy(float* src, size_t sz)
        {
            assert(sz <= data.size());
            std::copy(src, src + sz, data.begin());
        }

        void roll(const fmat& m)
        {
            assert(rows == 1 && m.rows == 1);
            assert(cols >= m.cols);
            std::memmove(data.data(), data.data() + m.data.size(), (cols - m.data.size()) * sizeof(float));
            std::copy(m.data.begin(), m.data.end(), data.data() + cols - m.data.size());
        }

        void roll(float v)
        {
            assert(rows == 1);
            std::memmove(data.data(), data.data() + 1, (cols - 1) * sizeof(float));
            data[cols - 1] = v;
        }

        void transpose()
        {
            fmat tmp(cols, rows);
            for (size_t r = 0; r < rows; ++r)
            {
                for (size_t c = 0; c < cols; ++c)
                {
                    tmp(c, r) = (*this)(r, c);
                }
            }
            *this = tmp;
        }

        void weight(const fmat& m) {
            assert(rows == 1 && m.rows == 1);
            assert(rows == m.rows && cols == m.cols);
            for (size_t c = 0; c < cols; ++c)
            {
                (*this)[c] *= m.data[c];
            }
        }
    };

    struct spectrum
    {
        static void gen_hann_window(float *window, size_t sz);

        spectrum(size_t fft_size_);
        void calc_mag(const float *sig, float *spec, float fs = 48000.0f);

        size_t fft_size;
    };
    
}
