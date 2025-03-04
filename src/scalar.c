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

/* scalar.c: defines functions that extract a feature as a single value from an input vector */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>

#ifndef DBL_MAX
#include <float.h> /* on Linux DBL_MAX is in float.h */
#endif

#include "dywapitchtrack/dywapitchtrack.h"

#include "xtract/libxtract.h"
#include "xtract/xtract_helper.h"
#include "xtract_macros_private.h"
#include "xtract_globals_private.h"

int xtract_mean(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N;

    *result = 0.0;

    while(n--)
        *result += data[n];

    *result /= N;

    return XTRACT_SUCCESS;
}

int xtract_variance(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N;
    const real_t arg0 = *(real_t *)argv;

    *result = 0.0;

    while(n--)
        *result += XTRACT_SQ(data[n] - arg0);

    *result /= (N - 1);

    return XTRACT_SUCCESS;
}

int xtract_standard_deviation(const real_t *data, const int N, const void *argv, real_t *result)
{

    *result = sqrt(*(real_t *)argv);

    return XTRACT_SUCCESS;
}

int xtract_average_deviation(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N;
    const real_t arg0 = *(real_t *)argv;

    *result = 0.0;

    while(n--)
        *result += fabs(data[n] - arg0);

    *result /= N;

    return XTRACT_SUCCESS;
}

int xtract_skewness(const real_t *data, const int N, const void *argv,  real_t *result)
{

    int n = N;

    real_t temp = 0.0;
    const real_t arg0 = ((real_t *)argv)[0];
    const real_t arg1 = ((real_t *)argv)[1];

    if (((real_t *)argv)[1] == 0)
    {
      *result = 0.0;
      return XTRACT_NO_RESULT;
    }

    *result = 0.0;

    if (arg1 == 0)
    {
      return XTRACT_NO_RESULT;
    }

    while(n--)
    {
        temp = (data[n] - arg0) / arg1;
        *result += XTRACT_POW3(temp);
    }

    *result /= N;


    return XTRACT_SUCCESS;
}

int xtract_kurtosis(const real_t *data, const int N, const void *argv,  real_t *result)
{

    int n = N;

    real_t temp = 0.0;
    const real_t arg0 = ((real_t *)argv)[0];
    const real_t arg1 = ((real_t *)argv)[1];

    if (arg1 == 0)
    {
        *result = 0.0;
        return XTRACT_NO_RESULT;
    }

    *result = 0.0;

    while(n--)
    {
        temp = (data[n] - arg0) / arg1;
        *result += XTRACT_POW4(temp);
    }

    *result /= N;
    *result -= 3.0;

    return XTRACT_SUCCESS;
}

int xtract_spectral_centroid(const real_t *data, const int N, const void *argv,  real_t *result)
{

    int n = (N >> 1);

    const real_t *freqs, *amps;
    real_t FA = 0.0, A = 0.0;

    amps = data;
    freqs = data + n;

    while(n--)
    {
        FA += freqs[n] * amps[n];
        A += amps[n];
    }

    if(A == 0.0)
        *result = 0.0;
    else
        *result = FA / A;

    return XTRACT_SUCCESS;
}

int xtract_spectral_mean(const real_t *data, const int N, const void *argv, real_t *result)
{

    return xtract_spectral_centroid(data, N, argv, result);

}

int xtract_spectral_variance(const real_t *data, const int N, const void *argv, real_t *result)
{

    int m;
    real_t A = 0.0;
    const real_t *freqs, *amps;
    const real_t arg0 = *(real_t *)argv;

    m = N >> 1;

    amps = data;
    freqs = data + m;

    *result = 0.0;

    while(m--)
    {
        A += amps[m];
        *result += XTRACT_SQ(freqs[m] - arg0) * amps[m];
    }

    if (A == 0.0)
    {
        *result = 0.0;
        return XTRACT_NO_RESULT;
    }
    *result = *result / A;

    return XTRACT_SUCCESS;
}

int xtract_spectral_standard_deviation(const real_t *data, const int N, const void *argv, real_t *result)
{

    *result = sqrt(*(real_t *)argv);

    return XTRACT_SUCCESS;
}

/*int xtract_spectral_average_deviation(const real_t *data, const int N, const void *argv, real_t *result){

    int m;
    real_t A = 0.0;
    const real_t *freqs, *amps;

    m = N >> 1;

    amps = data;
    freqs = data + m;

    *result = 0.0;

    while(m--){
        A += amps[m];
        *result += fabs((amps[m] * freqs[m]) - *(real_t *)argv);
    }

    *result /= A;

    return XTRACT_SUCCESS;
}*/

int xtract_spectral_skewness(const real_t *data, const int N, const void *argv,  real_t *result)
{

    int m;
    const real_t *freqs, *amps;
    const real_t arg0 = ((real_t *)argv)[0];
    const real_t arg1 = ((real_t *)argv)[1];

    *result = 0.0;

    if (arg1 == 0.0)
    {
        return XTRACT_NO_RESULT;
    }

    if (((real_t *)argv)[1] == 0.0)
    {
        *result = 0.0;
        return XTRACT_NO_RESULT;
    }

    m = N >> 1;

    amps = data;
    freqs = data + m;

    while(m--)
        *result += XTRACT_POW3(freqs[m] - arg0) * amps[m];

    *result /= XTRACT_POW3(arg1);

    return XTRACT_SUCCESS;
}

int xtract_spectral_kurtosis(const real_t *data, const int N, const void *argv,  real_t *result)
{

    int m;
    const real_t *freqs, *amps;
    const real_t arg0 = ((real_t *)argv)[0];
    const real_t arg1 = ((real_t *)argv)[1];

    if (((real_t *)argv)[1] == 0.0)
    {
        *result = 0.0;
        return XTRACT_NO_RESULT;
    }

    m = N >> 1;

    amps = data;
    freqs = data + m;

    *result = 0.0;

    while(m--)
        *result += XTRACT_POW4(freqs[m] - arg0) * amps[m];

    *result /= XTRACT_POW4(arg1);
    *result -= 3.0;

    return XTRACT_SUCCESS;
}

int xtract_irregularity_k(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n,
        M = N - 1;

    *result = 0.0;

    for(n = 1; n < M; n++)
        *result += fabs(data[n] - (data[n-1] + data[n] + data[n+1]) / 3.0);

    return XTRACT_SUCCESS;
}

int xtract_irregularity_j(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N - 1;

    real_t num = 0.0, den = 0.0;

    while(n--)
    {
        num += XTRACT_SQ(data[n] - data[n+1]);
        den += XTRACT_SQ(data[n]);
    }

    *result = num / den;

    return XTRACT_SUCCESS;
}

int xtract_tristimulus_1(const real_t *data, const int N, const void *argv, real_t *result)
{
    int n = N >> 1, i;
    real_t den = 0.0, p1 = 0.0, fund = 0.0, temp = 0.0, h = 0.0;
    const real_t *freqs;

    fund = *(real_t *)argv;
    freqs = data + n;

    for(i = 0; i < n; i++)
    {
        if((temp = data[i]))
        {
            den += temp;
            h = floor(freqs[i] / fund + 0.5);
            if(h > 0 && h < 2 && (int)h == 1)
                p1 += temp;
        }
    }

    if(den == 0.0 || p1 == 0.0)
    {
        *result = 0.0;
        return XTRACT_NO_RESULT;
    }
    else
    {
        *result = p1 / den;
        return XTRACT_SUCCESS;
    }
}

int xtract_tristimulus_2(const real_t *data, const int N, const void *argv, real_t *result)
{
    int n = N >> 1, i;
    real_t den, p2, p3, p4, ps, fund, temp, h;
    den = p2 = p3 = p4 = ps = fund = temp = h = 0.0;
    const real_t *freqs;

    fund = *(real_t *)argv;
    freqs = data + n;

    for(i = 0; i < n; i++)
    {
        if((temp = data[i]))
        {
            den += temp;
            h = floor(freqs[i] / fund + 0.5);
            if (h > 1 && h < 5)
            {
                switch ((int)h)
                {
                    case 2:
                        p2 += temp;
                    break;

                    case 3:
                        p3 += temp;
                    break;

                    case 4:
                        p4 += temp;

                    default:
                        break;
                }
            }
        }
    }

    ps = p2 + p3 + p4;

    if(den == 0.0 || ps == 0.0)
    {
        *result = 0.0;
        return XTRACT_NO_RESULT;
    }
    else
    {
        *result = ps / den;
        return XTRACT_SUCCESS;
    }

}

int xtract_tristimulus_3(const real_t *data, const int N, const void *argv, real_t *result)
{
    int n = N >> 1, i;
    real_t den = 0.0, num = 0.0, fund = 0.0, temp = 0.0, h = 0.0;
    const real_t *freqs;

    fund = *(real_t *)argv;
    freqs = data + n;

    for(i = 0; i < n; i++)
    {
        if((temp = data[i]))
        {
            den += temp;
            h = (int)floor(freqs[i] / fund + 0.5);
            if(h >= 5)
                num += temp;
        }
    }

    if(den == 0.0 || num == 0.0)
    {
        *result = 0.0;
        return XTRACT_NO_RESULT;
    }
    else
    {
        *result = num / den;
        return XTRACT_SUCCESS;
    }
}

int xtract_smoothness(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n; 
    int M = N - 1;
    real_t prev = 0.0;
    real_t current = 0.0;
    real_t next = 0.0;
    real_t temp = 0.0;

    for(n = 1; n < M; n++)
    {
        if(n == 1)
        {
            prev = data[n-1] <= 0 ? XTRACT_LOG_LIMIT_DB : log(data[n-1]);
            current = data[n] <= 0 ? XTRACT_LOG_LIMIT_DB : log(data[n]);
        }
        else
        {
            prev = current;
            current = next;
        }
        
        next = data[n+1] <= 0 ? XTRACT_LOG_LIMIT_DB : log(data[n+1]);
        
        temp += fabs(20.0 * current - (20.0 * prev +
                         20.0 * current + 20.0 * next) / 3.0);
    }

    *result = temp;

    return XTRACT_SUCCESS;
}

int xtract_spread(const real_t *data, const int N, const void *argv, real_t *result)
{

    return xtract_spectral_variance(data, N, argv, result);
}

int xtract_zcr(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N;
    int count = 0;

    for(n = 1; n < N; n++)
        if(data[n] * data[n-1] < 0) count++;

    *result = (real_t)count / N;

    return XTRACT_SUCCESS;
}

int xtract_rolloff(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N;
    real_t pivot, temp, percentile;

    pivot = temp = 0.0;
    percentile = ((real_t *)argv)[1];

    while(n--) pivot += data[n];

    pivot *= percentile / 100.0;

    for(n = 0; temp < pivot; n++)
        temp += data[n];

    *result = n * ((real_t *)argv)[0];
    /* *result = (n / (real_t)N) * (((real_t *)argv)[1] * .5); */

    return XTRACT_SUCCESS;
}

int xtract_loudness(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N, rv;

    *result = 0.0;

    if(n > XTRACT_BARK_BANDS)
    {
        n = XTRACT_BARK_BANDS;
        rv = XTRACT_BAD_VECTOR_SIZE;
    }
    else
        rv = XTRACT_SUCCESS;

    while(n--)
    {
        // The first bark coefficients is negative and makes the result N/A
        if (n > 0)
        {
            *result += pow(data[n], 0.23);
        }
    }

    return rv;
}

int xtract_flatness(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n, count, denormal_found;

    real_t num, den, temp;

    num = 1.0;
    den = temp = 0.0;

    denormal_found = 0;
    count = 0;

    for(n = 0; n < N; n++)
    {
        if((temp = data[n]) != 0.0)
        {
            if (xtract_is_denormal(num))
            {
                denormal_found = 1;
                break;
            }
            num *= temp;
            den += temp;
            count++;
        }
    }

    if(!count)
    {
        *result = 0.0;
        return XTRACT_NO_RESULT;
    }

    num = pow(num, 1.0 / (real_t)N);
    den /= (real_t)N;


    *result = (real_t) (num / den);

    if(denormal_found)
        return XTRACT_DENORMAL_FOUND;
    else
        return XTRACT_SUCCESS;

}

int xtract_flatness_db(const real_t *data, const int N, const void *argv, real_t *result)
{

    real_t flatness;

    flatness = *(real_t *)argv;

    if (flatness <= 0)
        flatness = XTRACT_LOG_LIMIT;

    *result = 10 * log10(flatness);

    return XTRACT_SUCCESS;

}

int xtract_tonality(const real_t *data, const int N, const void *argv, real_t *result)
{

    real_t sfmdb;

    sfmdb = *(real_t *)argv;

    *result = XTRACT_MIN(sfmdb / -60.0, 1);

    return XTRACT_SUCCESS;
}

int xtract_crest(const real_t *data, const int N, const void *argv, real_t *result)
{

    real_t max, mean;

    max = mean = 0.0;

    max = *(real_t *)argv;
    mean = *((real_t *)argv+1);

    *result = max / mean;

    return XTRACT_SUCCESS;

}

int xtract_noisiness(const real_t *data, const int N, const void *argv, real_t *result)
{

    real_t h, i, p; /*harmonics, inharmonics, partials */

    i = p = h = 0.0;

    h = *(real_t *)argv;
    p = *((real_t *)argv+1);

    if (p == 0)
    {
      *result = 0;
      return XTRACT_NO_RESULT;
    }

    i = p - h;

    *result = i / p;

    return XTRACT_SUCCESS;

}

int xtract_rms_amplitude(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N;

    *result = 0.0;

    while(n--) *result += XTRACT_SQ(data[n]);

    *result = sqrt(*result / (real_t)N);

    return XTRACT_SUCCESS;
}

int xtract_spectral_inharmonicity(const real_t *data, const int N, const void *argv, real_t *result)
{
    int n = N >> 1, h = 0;
    real_t num = 0.0, den = 0.0, fund;
    const real_t *freqs, *amps;

    fund = *(real_t *)argv;
    amps = data;
    freqs = data + n;

    if (fund == 0)
    {
      *result = 0;
      return XTRACT_NO_RESULT;
    }

    while(n--)
    {
        if(amps[n])
        {
            h = (int)floor(freqs[n] / fund + 0.5);
            num += fabs(freqs[n] - h * fund) * XTRACT_SQ(amps[n]);
            den += XTRACT_SQ(amps[n]);
        }
    }

    if (den == 0)
    {
      *result = 0;
      return XTRACT_NO_RESULT;
    }

    *result = (2 * num) / (fund * den);

    return XTRACT_SUCCESS;
}


int xtract_power(const real_t *data, const int N, const void *argv, real_t *result)
{

    return XTRACT_FEATURE_NOT_IMPLEMENTED;

}

int xtract_odd_even_ratio(const real_t *data, const int N, const void *argv, real_t *result)
{
    int n = N >> 1, h = 0;
    real_t odd = 0.0, even = 0.0, fund, temp;
    const real_t *freqs;

    fund = *(real_t *)argv;
    freqs = data + n;

    while(n--)
    {
        if((temp = data[n]))
        {
            h = (int)floor(freqs[n] / fund + 0.5);
            if(XTRACT_IS_ODD(h))
            {
                odd += temp;
            }
            else
            {
                even += temp;
            }
        }
    }

    if(odd == 0.0 || even == 0.0)
    {
        *result = 0.0;
        return XTRACT_NO_RESULT;
    }
    else
    {
        *result = odd / even;
        return XTRACT_SUCCESS;
    }
}

int xtract_sharpness(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N, rv;
    real_t sl, g; /* sl = specific loudness */
    real_t temp;

    sl = g = 0.0;
    temp = 0.0;

    if(n > XTRACT_BARK_BANDS)
        rv = XTRACT_BAD_VECTOR_SIZE;
    else
        rv = XTRACT_SUCCESS;


    while(n--)
    {
        sl = pow(data[n], 0.23);
        g = (n < 15 ? 1.0 : 0.066 * exp(0.171 * n));
        temp += n * g * sl;
    }

    temp = 0.11 * temp / (real_t)N;
    *result = (real_t)temp;

    return rv;

}

int xtract_spectral_slope(const real_t *data, const int N, const void *argv, real_t *result)
{

    const real_t *freqs, *amps;
    real_t f, a,
          F, A, FA, FXTRACT_SQ; /* sums of freqs, amps, freq * amps, freq squared */
    int n, M;

    F = A = FA = FXTRACT_SQ = 0.0;
    n = M = N >> 1;

    amps = data;
    freqs = data + n;

    while(n--)
    {
        f = freqs[n];
        a = amps[n];
        F += f;
        A += a;
        FA += f * a;
        FXTRACT_SQ += f * f;
    }

    real_t temp = (real_t)M * FXTRACT_SQ - F * F;

    if (A == 0 || temp == 0)
    {
        *result = 0.0;
        return XTRACT_NO_RESULT;
    }

    *result = (1.0 / A) * ((real_t)M * FA - F * A) / temp;

    return XTRACT_SUCCESS;

}

int xtract_lowest_value(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N;

    *result = DBL_MAX;

    while(n--)
    {
        if(data[n] > *(real_t *)argv)
            *result = XTRACT_MIN(*result, data[n]);
    }

    if (*result == DBL_MAX)
        return XTRACT_NO_RESULT;
        
    return XTRACT_SUCCESS;
}

int xtract_highest_value(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N;

    *result = data[--n];

    while(n--)
        *result = XTRACT_MAX(*result, data[n]);

    return XTRACT_SUCCESS;
}


int xtract_sum(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N;

    *result = 0.0;

    while(n--)
        *result += *data++;

    return XTRACT_SUCCESS;

}

int xtract_nonzero_count(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N;

    *result = 0.0;

    while(n--)
        *result += (*data++ ? 1 : 0);

    return XTRACT_SUCCESS;

}

int xtract_hps(const real_t *data, const int N, const void *argv, real_t *result)
{
    int n, M, i, peak_index, position1_lwr;
    real_t tempProduct, peak, largest1_lwr, ratio1;

    n = N / 2;

    M = (int)ceil(n / 3.0);

    if (M <= 1)
    {
        /* Input data is too short. */
        *result = 0;
        return XTRACT_NO_RESULT;
    }

    peak_index = 0;

    tempProduct = peak = 0;
    for (i = 0; i < M; ++i)
    {
        tempProduct = data [i] * data [i * 2] * data [i * 3];

        if (tempProduct > peak)
        {
            peak = tempProduct;
            peak_index = i;
        }
    }

    largest1_lwr = position1_lwr = 0;

    for(i = 0; i < N; ++i)
    {
        if(data[i] > largest1_lwr && i != peak_index)
        {
            largest1_lwr = data[i];
            position1_lwr = i;
        }
    }

    ratio1 = data[position1_lwr] / data[peak_index];

    if(position1_lwr > peak_index * 0.4 && position1_lwr <
            peak_index * 0.6 && ratio1 > 0.1)
        peak_index = position1_lwr;

    *result = data [n + peak_index];

    return XTRACT_SUCCESS;
}

int xtract_f0(const real_t *data, const int N, const void *argv, real_t *result)
{

    int M, tau, n;
    real_t sr;
    size_t bytes;
    real_t f0, err_tau_1, err_tau_x, array_max,
          threshold_peak, threshold_centre,
          *input;

    sr = *(real_t *)argv;
    if(sr == 0)
        sr = 44100.0;

    input = (real_t*)malloc(bytes = N * sizeof(real_t));
    input = (real_t*)memcpy(input, data, bytes);
    /*  threshold_peak = *((real_t *)argv+1);
    threshold_centre = *((real_t *)argv+2);
    printf("peak: %.2\tcentre: %.2\n", threshold_peak, threshold_centre);*/
    /* add temporary dynamic control over thresholds to test clipping effects */

    /* FIX: tweak and  make into macros */
    threshold_peak = .8;
    threshold_centre = .3;
    M = N >> 1;
    err_tau_1 = 0;
    array_max = 0;

    /* Find the array max */
    for(n = 0; n < N; n++)
    {
        if (input[n] > array_max)
            array_max = input[n];
    }

    threshold_peak *= array_max;

    /* peak clip */
    for(n = 0; n < N; n++)
    {
        if(input[n] > threshold_peak)
            input[n] = threshold_peak;
        else if(input[n] < -threshold_peak)
            input[n] = -threshold_peak;
    }

    threshold_centre *= array_max;

    /* Centre clip */
    for(n = 0; n < N; n++)
    {
        if (input[n] < threshold_centre)
            input[n] = 0;
        else
            input[n] -= threshold_centre;
    }

    /* Estimate fundamental freq */
    for (n = 1; n < M; n++)
        err_tau_1 = err_tau_1 + fabs(input[n] - input[n+1]);
    /* FIX: this doesn't pose too much load if it returns 'early', but if it can't find f0, load can be significant for larger block sizes M^2 iterations! */
    for (tau = 2; tau < M; tau++)
    {
        err_tau_x = 0;
        for (n = 1; n < M; n++)
        {
            err_tau_x = err_tau_x + fabs(input[n] - input[n+tau]);
        }
        if (err_tau_x < err_tau_1)
        {
            f0 = sr / (tau + (err_tau_x / err_tau_1));
            *result = f0;
            free(input);
            return XTRACT_SUCCESS;
        }
    }
    *result = -0;
    free(input);
    return XTRACT_NO_RESULT;
}

int xtract_failsafe_f0(const real_t *data, const int N, const void *argv, real_t *result)
{

    real_t *spectrum = NULL, argf[4], *peaks = NULL, return_code, sr;

    return_code = xtract_f0(data, N, argv, result);

    if(return_code == XTRACT_NO_RESULT)
    {
        sr = *(real_t *)argv;
        if(sr == 0)
            sr = 44100.0;
        spectrum = (real_t *)malloc(N * sizeof(real_t));
        memset(spectrum, 0, N * sizeof(real_t));
        peaks = (real_t *)malloc(N * sizeof(real_t));
        argf[0] = sr / N;
        argf[1] = XTRACT_MAGNITUDE_SPECTRUM;
        argf[2] = 0.0f;
        argf[3] = 0.0f;
        xtract_spectrum(data, N, argf, spectrum);
        argf[1] = 10.0;
        xtract_peak_spectrum(spectrum, N >> 1, argf, peaks);
        argf[0] = 0.0;
        xtract_lowest_value(peaks+(N >> 1), N >> 1, argf, result);

        free(spectrum);
        free(peaks);
    }

    return XTRACT_SUCCESS;

}

int xtract_wavelet_f0(const real_t *data, const int N, const void *argv, real_t *result)
{
    /* real_t sr = *(real_t *)argv; */

    *result = dywapitch_computepitch(&wavelet_f0_state, data, 0, N);

    if (*result == 0.0)
    {
        return XTRACT_NO_RESULT;
    }

    return XTRACT_SUCCESS;
}

int xtract_midicent(const real_t *data, const int N, const void *argv, real_t *result)
{
    real_t f0 = *(real_t *)argv;
    real_t note = 0.0;
      
    note = 69 + log(f0 / 440.f) * 17.31234;
    note *= 100;
    note = floor( 0.5f + note ); // replace -> round(note);

    *result = note;
    
    if (note > 12700 || note < 0)
    {
        return XTRACT_ARGUMENT_ERROR;
    }
    
    return XTRACT_SUCCESS;
}

int xtract_peak(const real_t *data, const int N, const void *argv, real_t *result)
{
    real_t threshold = *(real_t *)argv;
    real_t current = data[N - 1];
    real_t average = 0.0;
    real_t maximum = -DBL_MAX;
    
    for (uint32_t n = 0; n < (uint32_t)N; ++n)
    {
        average += data[n];
        if (data[n] > maximum)
        {
            maximum = data[n];
        }
    }
    
    average /= (real_t)N;
        
    if (current != maximum)
    {
        return XTRACT_NO_RESULT;
    }
    
    if (current < average + threshold)
    {
        return XTRACT_NO_RESULT;
    }
    
    *result = current;
    
    return XTRACT_SUCCESS;
    
}


