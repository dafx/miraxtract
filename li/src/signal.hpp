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

    static void iir2(const float *x, float *y, size_t n, const float *b, const float *a, float *d)
    {
        /*
         a[0]*y[n] = b[0] * x[n]               + d[0][n-1]
           d[0][n] = b[1] * x[n] - a[1] * y[n] + d[1][n-1]
           d[1][n] = b[2] * x[n] - a[2] * y[n]
         */
        const float one_over_a0 = 1.0f / a[0];
        for (size_t ix = 0; ix < n; ix++)
        {
            const float xx = x[ix];
            y[ix] = b[0] * xx + d[0];
            y[ix] *= one_over_a0;
            d[0] = b[1] * xx - a[1] * y[ix] + d[1];
            d[1] = b[2] * xx - a[2] * y[ix];
        }
    }

    struct sosfilt
    {
        const float *coeff; // 6 * num_sections coefficients
        float *zi;
        fvec zi_vec; // 2 * num_sections initial conditions
        size_t num_sections;

        sosfilt(const float *coeff_, const float *zi_, size_t num_sections_)
            : coeff(coeff_),
              zi_vec(zi_, zi_ + (num_sections_ * 2)),
              num_sections(num_sections_)
        {
            zi = zi_vec.data();
        }

        /**
         * @brief IIR filters in second-order sections.
         * This is the counterpart of scipy.signal.sosfilt .
         * @param input Input signal
         * @param output Output signal. Can be the same as input for in place
         * @param x_size Minimum size of input and output signal
         */
        void run(const float *input, const size_t size, float *output)
        {
            assert(num_sections > 0);

            iir2(input, output, size, coeff, coeff + 3, zi);

            for (size_t sect = 1; sect < num_sections; sect++)
            {
                iir2(
                    output,
                    output,
                    size,
                    coeff + sect * 6,
                    coeff + sect * 6 + 3,
                    zi + sect * 2);
            }
        }

        void init(float x0)
        {
            for (size_t sect = 0; sect < num_sections; sect++)
            {
                zi[sect * 2] *= x0;
                zi[sect * 2 + 1] *= x0;
            }
        }
    };

    struct deci_by_4
    {
        sosfilt sos;
        bool init_done = false;

        deci_by_4();
        void run(float *input, const size_t size, float *output);
    };
}
