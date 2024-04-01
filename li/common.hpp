/*
 * Copyright (c) 2024, Yan Li
 * All rights reserved.
 */

#pragma once

#include <vector>

namespace li
{
    typedef std::vector<float> fvec;

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

        void roll(fvec& v)
        {
            assert(rows == 1 && cols >= v.size());            
            std::memmove(data.data(), data.data() + v.size(), (cols - v.size()) * sizeof(float));
            std::copy(v.begin(), v.end(), data.data() + cols - v.size());
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
    };

}
