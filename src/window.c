/*
 * Copyright (C) 2012 Jamie Bullock
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *
 */

/* window.c: defines window generation functions (formulae courtesy of Wikipedia (http://en.wikipedia.org/wiki/Window_function) */

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#include "xtract_window_private.h"

void gauss(real_t *window, const int N, const real_t sd)
{

    int n;
    const real_t M = N - 1;
    real_t num,
          den,
          exponent;

    for (n = 0; n < N; n++)
    {

        num = n - M / 2.0;
        den = sd * M / 2.0;

        exponent = -0.5 * pow(num / den, 2);

        window[n] = exp(exponent);

    }
}

void hamming(real_t *window, const int N)
{

    int n;
    const real_t M = N - 1;

    for (n = 0; n < N; n++)
        window[n] = 0.53836 - (0.46164 * cos(2.0 * M_PI * (real_t)n / M));

}

void hann(real_t *window, const int N)
{

    int n;
    const real_t M = N - 1;

    for (n = 0; n < N; n++)
        window[n] = 0.5 * (1.0 - cos(2.0 * M_PI * (real_t)n / M));

}

void bartlett(real_t *window, const int N)
{

    int n;
    const real_t M = N - 1;

    for (n = 0; n < N; n++)
        window[n] = 2.0 / M * (M / 2.0 - fabs(n - M / 2.0));

}

void triangular(real_t *window, const int N)
{

    int n;
    const real_t M = N - 1;

    for (n = 0; n < N; n++)
        window[n] = 2.0 / N * (N / 2.0 - fabs(n - M / 2.0));
}

void bartlett_hann(real_t *window, const int N)
{

    int n;
    const real_t M = N - 1,
                a0 = 0.62,
                a1 = 0.5,
                a2 = 0.38;
    real_t term1 = 0.0,
          term2 = 0.0;

    for (n = 0; n < N; n++)
    {

        term1 = a1 * fabs(n / M - 0.5);
        term2 = a2 * cos(2.0 * M_PI * (real_t)n / M);

        window[n] = a0 - term1 - term2;
    }
}

void blackman(real_t *window, const int N)
{

    int n;
    const real_t M = N - 1,
                a0 = 0.42,
                a1 = 0.5,
                a2 = 0.08;
    real_t term1 = 0.0,
          term2 = 0.0;

    for (n = 0; n < N; n++)
    {

        term1 = a1 * cos(2.0 * M_PI * (real_t)n / M);
        term2 = a2 * cos(4.0 * M_PI * (real_t)n / M);

        window[n] = a0 - term1 + term2;
    }
}

#define BIZ_EPSILON 1E-21 // Max error acceptable 

/* Based on code from mplayer window.c, and somewhat beyond me */
real_t besselI0(real_t x)
{

    real_t temp;
    real_t sum   = 1.0;
    real_t u     = 1.0;
    real_t halfx = x/2.0;
    int      n     = 1;

    do
    {

        temp = halfx/(real_t)n;
        u *=temp * temp;
        sum += u;
        n++;

    }
    while (u >= BIZ_EPSILON * sum);

    return(sum);

}

void kaiser(real_t *window, const int N, const real_t alpha)
{

    int n;
    const real_t M = N - 1;
    real_t num;

    for (n = 0; n < N; n++)
    {

        num = besselI0(alpha * sqrt(1.0 - pow((2.0 * n / M - 1), 2)));
        window[n] = num / besselI0(alpha);

    }
}

void blackman_harris(real_t *window, const int N)
{

    int n;
    const real_t M = N - 1,
                a0 = 0.35875,
                a1 = 0.48829,
                a2 = 0.14128,
                a3 = 0.01168;
    real_t term1 = 0.0,
          term2 = 0.0,
          term3 = 0.0;

    for (n = 0; n < N; n++)
    {

        term1 = a1 * cos(2.0 * M_PI * n / M);
        term2 = a2 * cos(4.0 * M_PI * n / M);
        term3 = a3 * cos(6.0 * M_PI * n / M);

        window[n] = a0 - term1 + term2 - term3;
    }
}
