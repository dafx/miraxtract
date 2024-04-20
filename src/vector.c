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

/* xtract_vector.c: defines functions that extract a feature as a single value from an input vector */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

#include "fft.h"

#include "xtract/libxtract.h"
#include "xtract_macros_private.h"
#include "xtract_globals_private.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

thread_local real_t** dct_cos_table = NULL;
thread_local int dct_cos_table_dim = 0;

int xtract_spectrum(const real_t *data, const int N, const void *argv, real_t *result)
{

    int vector     = 0;
    int withDC     = 0;
    int normalise  = 0;
    real_t q        = 0.0;
    real_t temp     = 0.0;
    real_t max      = 0.0;
    real_t NxN      = XTRACT_SQ((real_t)N);
    real_t real = 0.0;
    real_t imag = 0.0;
    unsigned int n = 0;
    unsigned int m = 0;
    unsigned int M = N >> 1;
#ifdef USE_OOURA
    real_t *fft = NULL;
#else 
    DSPDoubleSplitComplex *fft = NULL;
#endif

    q = *(real_t *)argv;
    vector = (int)*((real_t *)argv+1);
    withDC = (int)*((real_t *)argv+2);
    normalise = (int)*((real_t *)argv+3);

    XTRACT_CHECK_q;
#ifdef USE_OOURA
    if(!ooura_data_spectrum.initialised)
#else
    if(!vdsp_data_spectrum.initialised)
#endif
    {
        fprintf(stderr,
                "libxtract: error: xtract_spectrum() failed, "
                "fft data unitialised.\n");
        return XTRACT_NO_RESULT;
    }

#ifdef USE_OOURA
    /* ooura is in-place
     * the output format is
     * a[0] - DC, a[1] - nyquist, a[2...N-1] - remaining bins
     */
    fft = (real_t*)malloc(N * sizeof(real_t));
    assert(fft != NULL);
    memcpy(fft, data, N * sizeof(real_t));

    rdft(N, 1, fft, ooura_data_spectrum.ooura_ip, 
            ooura_data_spectrum.ooura_w);
#else
    fft = &vdsp_data_spectrum.fft;
    vDSP_ctozD((DSPDoubleComplex *)data, 2, fft, 1, N >> 1);
    vDSP_fft_zripD(vdsp_data_spectrum.setup, fft, 1, 
            vdsp_data_spectrum.log2N, FFT_FORWARD);
#endif

    switch(vector)
    {

        case XTRACT_LOG_MAGNITUDE_SPECTRUM:
        for(n = 0, m = 0; m < M; ++n, ++m)
        {
            if(n==0 && !withDC) /* discard DC and keep Nyquist */
            {
                ++n;
            }
#ifdef USE_OOURA
			/*
            if(n==1 && withDC) // discard Nyquist
            {
                ++n;
            }
			*/
			// OOURA discards the always 0 imaginary of DC and Nyquists
			if (n == M && !withDC)
			{
				real = fft[1];
				imag = 0.0;
			}
			else if (n == 0 && withDC) {
				real = fft[0];
				imag = 0.0;
			}
			else {
				real = fft[n * 2];
				imag = fft[n * 2 + 1];
			}
#else
			if (n == M && !withDC)
			{
				real = fft->imagp[0];
				imag = 0.0;
		}
			else if (n == 0 && withDC) {
				real = fft->realp[0];
				imag = 0.0;
			}
			else {
				real = fft->realp[n];
				imag = fft->imagp[n];
			}
#endif

            temp = XTRACT_SQ(real) + XTRACT_SQ(imag);
            if (temp > XTRACT_LOG_LIMIT)
            {
                temp = log(sqrt(temp) / (real_t)N);
            }
            else
            {
                temp = XTRACT_LOG_LIMIT_DB;
            }
            result[m] =
                /* Scaling */
                (temp + XTRACT_DB_SCALE_OFFSET) /
                XTRACT_DB_SCALE_OFFSET;

            XTRACT_SET_FREQUENCY;
            XTRACT_GET_MAX;
        }
        break;

    case XTRACT_POWER_SPECTRUM:
        for(n = 0, m = 0; m < M; ++n, ++m)
        {
            if(n==0 && !withDC) /* discard DC and keep Nyquist */
            {
                ++n;
            }
#ifdef USE_OOURA
			/*
			if(n==1 && withDC) // discard Nyquist
			{
				++n;
			}
			*/
			// OOURA discards the always 0 imaginary of DC and Nyquists
			if (n == M && !withDC)
			{
				real = fft[1];
				imag = 0.0;
			}
			else if (n == 0 && withDC) {
				real = fft[0];
				imag = 0.0;
			}
			else {
				real = fft[n * 2];
				imag = fft[n * 2 + 1];
			}
#else
			if (n == M && !withDC)
			{
				real = fft->imagp[0];
				imag = 0.0;
		}
			else if (n == 0 && withDC) {
				real = fft->realp[0];
				imag = 0.0;
			}
			else {
				real = fft->realp[n];
				imag = fft->imagp[n];
			}
#endif

            result[m] = (XTRACT_SQ(real) + XTRACT_SQ(imag)) / NxN;
            XTRACT_SET_FREQUENCY;
            XTRACT_GET_MAX;
        }
        break;

    case XTRACT_LOG_POWER_SPECTRUM:
        for(n = 0, m = 0; m < M; ++n, ++m)
        {
            if(n==0 && !withDC) /* discard DC and keep Nyquist */
            {
                ++n;
            }
#ifdef USE_OOURA
			/*
			if(n==1 && withDC) // discard Nyquist
			{
				++n;
			}
			*/
			// OOURA discards the always 0 imaginary of DC and Nyquists
			if (n == M && !withDC)
			{
				real = fft[1];
				imag = 0.0;
			}
			else if (n == 0 && withDC) {
				real = fft[0];
				imag = 0.0;
			}
			else {
				real = fft[n * 2];
				imag = fft[n * 2 + 1];
			}
#else
			if (n == M && !withDC)
			{
				real = fft->imagp[0];
				imag = 0.0;
		}
			else if (n == 0 && withDC) {
				real = fft->realp[0];
				imag = 0.0;
			}
			else {
				real = fft->realp[n];
				imag = fft->imagp[n];
			}
#endif

            if ((temp = XTRACT_SQ(real) + XTRACT_SQ(imag)) >
                    XTRACT_LOG_LIMIT)
                temp = log(temp / NxN);
            else
                temp = XTRACT_LOG_LIMIT_DB;

            result[m] = (temp + XTRACT_DB_SCALE_OFFSET) /
                        XTRACT_DB_SCALE_OFFSET;
            XTRACT_SET_FREQUENCY;
            XTRACT_GET_MAX;
        }
        break;

    case XTRACT_SPECTRUM_COEFFICIENTS:
        for(n = 0, m = 0; m < M; ++n, ++m)
        {
            if(n==0 && !withDC) /* discard DC and keep Nyquist */
            {
                ++n;
            }
#ifdef USE_OOURA
			/*
			if(n==1 && withDC) // discard Nyquist
			{
				++n;
			}
			*/
			// OOURA discards the always 0 imaginary of DC and Nyquists
			if (n == M && !withDC)
			{
				real = fft[1];
				imag = 0.0;
			}
			else if (n == 0 && withDC) {
				real = fft[0];
				imag = 0.0;
			}
			else {
				real = fft[n * 2];
				imag = fft[n * 2 + 1];
			}
#else
			if (n == M && !withDC)
			{
				real = fft->imagp[0];
				imag = 0.0;
		    }
			else if (n == 0 && withDC) {
				real = fft->realp[0];
				imag = 0.0;
			}
			else {
				real = fft->realp[n];
				imag = fft->imagp[n];
			}
#endif
            result[m*2] = real;
            result[m*2+1] = imag;
            XTRACT_GET_MAX;
            }
        break;

    default:
        /* MAGNITUDE_SPECTRUM */
        for(n = 0, m = 0; m < M; ++n, ++m)
        {
            if(n==0 && !withDC) /* discard DC and keep Nyquist */
            {
                ++n;
            }
#ifdef USE_OOURA
			/*
			if(n==1 && withDC) // discard Nyquist
			{
				++n;
			}
			*/
			// OOURA discards the always 0 imaginary of DC and Nyquists
			if (n == M && !withDC)
			{
				real = fft[1];
				imag = 0.0;
			}
			else if (n == 0 && withDC) {
				real = fft[0];
				imag = 0.0;
			}
			else {
				real = fft[n * 2];
				imag = fft[n * 2 + 1];
			}
#else
            if (n == M && !withDC)
            {
                real = fft->imagp[0];
                imag = 0.0;
            }
            else if (n == 0 && withDC) {
                real = fft->realp[0];
                imag = 0.0;
            }
            else {
                real = fft->realp[n];
                imag = fft->imagp[n];
            }
#endif
            result[m*2] = real;
            result[m*2+1] = imag;
            XTRACT_SET_FREQUENCY;
            XTRACT_GET_MAX;
        }
        break;
    }

    if(normalise)
    {
        max += FLT_EPSILON;
        for (n = 0; n < M; n++)
            result[n] /= max;
    }

#ifdef USE_OOURA
    free(fft);
#endif

    return XTRACT_SUCCESS;
}

int xtract_autocorrelation_fft(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n        = 0;
    int M        = N << 1;

#ifdef USE_OOURA
    real_t *rfft = NULL;
#else
    DSPDoubleSplitComplex *fft = NULL;
    real_t M_double = 0.0;
#endif


#ifdef USE_OOURA
    /* Zero pad the input vector */
    rfft = (real_t *)calloc(M, sizeof(real_t));
    memcpy(rfft, data, N * sizeof(real_t));
    
    rdft(M, 1, rfft, ooura_data_autocorrelation_fft.ooura_ip, 
            ooura_data_autocorrelation_fft.ooura_w);

    for(n = 2; n < M; ++n)
    {
        rfft[n*2] = XTRACT_SQ(rfft[n*2]) + XTRACT_SQ(rfft[n*2+1]);
        rfft[n*2+1] = 0.0;
    }

    rfft[0] = XTRACT_SQ(rfft[0]);
    rfft[1] = XTRACT_SQ(rfft[1]);

    rdft(M, -1, rfft, ooura_data_autocorrelation_fft.ooura_ip,
            ooura_data_autocorrelation_fft.ooura_w);

#else
    /* vDSP has its own autocorrelation function, but it doesn't fit the 
     * LibXtract model, e.g. we can't guarantee it's going to use
     * an FFT for all values of N */
    fft = &vdsp_data_autocorrelation_fft.fft;
    vDSP_ctozD((DSPDoubleComplex *)data, 2, fft, 1, N);
    vDSP_fft_zripD(vdsp_data_autocorrelation_fft.setup, fft, 1, 
            vdsp_data_autocorrelation_fft.log2N, FFT_FORWARD);

    for(n = 0; n < N; ++n)
    {
        fft->realp[n] = XTRACT_SQ(fft->realp[n]) + XTRACT_SQ(fft->imagp[n]);
        fft->imagp[n] = 0.0;
    }

    vDSP_fft_zripD(vdsp_data_autocorrelation_fft.setup, fft, 1, 
            vdsp_data_autocorrelation_fft.log2N, FFT_INVERSE);
#endif

    /* Normalisation factor */
    M = M * N;

#ifdef USE_OOURA
    for(n = 0; n < N; n++)
        result[n] = rfft[n] / (real_t)M;
    free(rfft);
#else
    M_double = (real_t)M;
    vDSP_ztocD(fft, 1, (DOUBLE_COMPLEX *)result, 2, N);
    vDSP_vsdivD(result, 1, &M_double, result, 1, N);
#endif

    return XTRACT_SUCCESS;
}

int xtract_mfcc(const real_t *data, const int N, const void *argv, real_t *result)
{

    xtract_mel_filter *f;
    int n, filter;
    real_t* temp;

    f = (xtract_mel_filter *)argv;
    temp = calloc(f->n_filters, sizeof(real_t));
    for(filter = 0; filter < f->n_filters; filter++)
    {
        for(n = 0; n < N; n++)
        {
            if (f->filters[filter][n] != 0)
                temp[filter] += data[n] * f->filters[filter][n];
        }
        if (temp[filter] < XTRACT_LOG_LIMIT)
            temp[filter] = XTRACT_LOG_LIMIT_DB;
        else
            temp[filter] = log(temp[filter]);
    }

    xtract_dct(temp, f->n_filters, NULL, result);
    free(temp);

    return XTRACT_SUCCESS;
}

int xtract_mmbses(const real_t *data, const int N, const void *argv, real_t *result)
{
    xtract_mel_filter *f;
    int n, filter;
    real_t* real = (real_t*)malloc(sizeof(real_t)*N);
	real_t* imag = (real_t*)malloc(sizeof(real_t)*N);

    f = (xtract_mel_filter *)argv;

    for (filter = 0; filter < f->n_filters; filter++)
    {
        int count = 0;
        real_t realMean = 0, realVariance = 0;
        real_t imagMean = 0, imagVariance = 0;
        real_t covariance = 0;
        real_t energy = 0;

        result[filter] = 0.0;
        for(n = 0; n < N; n++)
        {
          real_t tempReal = data[n*2]*f->filters[filter][n];
          real_t tempImag = data[n*2+1]*f->filters[filter][n];

            if (f->filters[filter][n] != 0)
            {
                real[count] = tempReal;
                imag[count] = tempImag;
                count++;
            }
            energy += sqrt(XTRACT_SQ(tempReal) + XTRACT_SQ(tempImag)) / (real_t)N;
        }
        if (count == 0)
          continue;
        if (count == 1)
        {
          energy = (real_t)2*M_PI*energy;

          if (energy < XTRACT_LOG_LIMIT)
              result[filter] = XTRACT_LOG_LIMIT_DB;
          else
              result[filter] = log(energy);

          continue;
        }
        // Calculate the arithmetic means of real and imaginary parts
        for(n = 0; n < count; n++)
        {
            realMean += real[n] / count;
            imagMean += imag[n] / count;
        }
        // Calculate the variances of real and imaginary parts
        for(n = 0; n < count; n++)
        {
            realVariance += XTRACT_SQ(real[n]-realMean) / count;
            imagVariance += XTRACT_SQ(imag[n]-imagMean) / count;
        }
        // Calculate the covariance between real and imaginary parts
        for(n = 0; n < count; n++)
        {
            covariance += (real[n]-realMean)*(imag[n]-imagMean);
        }
        covariance /= (count-1);
        // Calculate the final Mel based Multi-Band Spectral Entropy Signature coefficients
        real_t temp = realVariance*imagVariance-XTRACT_SQ(covariance);

        if (temp < XTRACT_LOG_LIMIT)
            temp = XTRACT_LOG_LIMIT_DB;
        else
            temp = log(temp);

        energy = (real_t)2*M_PI*energy;
        if (energy < XTRACT_LOG_LIMIT)
            energy = XTRACT_LOG_LIMIT_DB;
        else
            energy = log(energy);

        result[filter] = energy+temp / 2;
    }
    free(real);
    free(imag);
    return XTRACT_SUCCESS;
}

int xtract_spectral_subband_centroids(const real_t *data, const int N, const void *argv, real_t *result)
{
    xtract_mel_filter *f = (xtract_mel_filter *)argv;
    int n, filter;
    const real_t *freqs = data;
    const real_t *amps = data+N;

    for (filter = 0; filter < f->n_filters; filter++)
    {
        real_t FA = 0.0, A = 0.0;

        for(n = 0; n < N; n++)
        {
            real_t Multiplier = amps[n]*f->filters[filter][n];

            FA += freqs[n]*f->filters[filter][n]*Multiplier;
            A += Multiplier;
        }
        if (FA == 0.0 || A == 0.0)
            result[filter] = 0;
        else
            result[filter] = FA / A;
    }
    return XTRACT_SUCCESS;
}

int xtract_dct(const real_t *data, const int N, const void *argv, real_t *result)
{
    int n, m;
    // Extra variable to hold a reference for the dct lookup table since
    // accessing the thread local storage is expensive.
    real_t** temp_dct_table;

    // Free the dct table if the cached dimension is different from the new dimension
    if (dct_cos_table != NULL && dct_cos_table_dim != N)
    {
        for (n = 0; n < N; ++n)
        {
          free(dct_cos_table[n]);
        }
        free(dct_cos_table);
        dct_cos_table = NULL;
    }
    // Allocate the dct cache table
    if (dct_cos_table == NULL)
    {
        dct_cos_table_dim = N;
        dct_cos_table = calloc(N, sizeof(real_t*));
        for (n = 0; n < N; ++n)
        {
            dct_cos_table[n] = calloc(N, sizeof(real_t));
            for (m = 1; m <= N; ++m)
            {
                dct_cos_table[n][m-1] = cos(M_PI * (n / (real_t)N)*(m - 0.5));
            }
        }
    }
    // Calculate the dct transformation
    temp_dct_table = dct_cos_table;
    memset(result, 0, N * sizeof(real_t));
    for (n = 0; n < N; ++n)
    {
        for (m = 0; m < N; ++m)
            result[n] += data[m]*temp_dct_table[n][m];
    }

    return XTRACT_SUCCESS;
}

int xtract_autocorrelation(const real_t *data, const int N, const void *argv, real_t *result)
{

    /* Naive time domain implementation  */

    int n = N, i;

    real_t corr;

    while(n--)
    {
        corr = 0;
        for(i = 0; i < N - n; i++)
        {
            corr += data[i] * data[i + n];
        }
        result[n] = corr / N;
    }

    return XTRACT_SUCCESS;
}

int xtract_amdf(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N, i;

    real_t md, temp;

    while(n--)
    {
        md = 0.0;
        for(i = 0; i < N - n; i++)
        {
            temp = data[i] - data[i + n];
            temp = (temp < 0 ? -temp : temp);
            md += temp;
        }
        result[n] = md / (real_t)N;
    }

    return XTRACT_SUCCESS;
}

int xtract_asdf(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = N, i;

    real_t sd;

    while(n--)
    {
        sd = 0.0;
        for(i = 0; i < N - n; i++)
        {
            /*sd = 1;*/
            sd += XTRACT_SQ(data[i] - data[i + n]);
        }
        result[n] = sd / (real_t)N;
    }

    return XTRACT_SUCCESS;
}

int xtract_bark_coefficients(const real_t *data, const int N, const void *argv, real_t *result)
{

    int *limits, band, n;

    limits = (int *)argv;

    for(band = 0; band < XTRACT_BARK_BANDS - 1; band++)
    {
        result[band] = 0.0;
        for(n = limits[band]; n < limits[band + 1] && n < N; n++)
            result[band] += data[n];
    }

    return XTRACT_SUCCESS;
}

int xtract_peak_spectrum(const real_t *data, const int N, const void *argv, real_t *result)
{

    real_t threshold, max, y, y2, y3, p, q;
    int n = N, rv = XTRACT_SUCCESS;

    threshold = max = y = y2 = y3 = p = q = 0.0;

    if(argv != NULL)
    {
        q = ((real_t *)argv)[0];
        threshold = ((real_t *)argv)[1];
    }
    else
        rv = XTRACT_BAD_ARGV;

    if(threshold < 0 || threshold > 100)
    {
        threshold = 0;
        rv = XTRACT_BAD_ARGV;
    }

    XTRACT_CHECK_q;

    threshold *= .01 * max;

    result[0] = 0;
    result[N] = 0;

    for(n = 1; n < N; n++)
    {
        if(data[n] >= threshold)
        {
            if(data[n] > data[n - 1] && n + 1 < N && data[n] > data[n + 1])
            {
                result[N + n] = q * (n + 1 + (p = .5 * ((y = data[n-1]) -
                                                    (y3 = data[n+1])) / (data[n - 1] - 2 *
                                                            (y2 = data[n]) + data[n + 1])));
                result[n] = y2 - .25 * (y - y3) * p;
            }
            else
            {
                result[n] = 0;
                result[N + n] = 0;
            }
        }
        else
        {
            result[n] = 0;
            result[N + n] = 0;
        }
    }

    return (rv ? rv : XTRACT_SUCCESS);
}

int xtract_harmonic_spectrum(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n = (N >> 1), M = n;

    const real_t *freqs, *amps;
    real_t f0, threshold, ratio, nearest, distance;

    amps = data;
    freqs = data + n;
    f0 = *((real_t *)argv);
    threshold = *((real_t *)argv+1);

    ratio = nearest = distance = 0.0;

    while(n--)
    {
        if(freqs[n])
        {
            ratio = freqs[n] / f0;
			nearest = floor( 0.5f + ratio);				// replace -> nearest = round(ratio);
			distance = fabs(nearest - ratio);
            if(distance > threshold)
                result[n] = result[M + n] = 0.0;
            else
            {
                result[n] = amps[n];
                result[M + n] = freqs[n];
            }
        }
        else
            result[n] = result[M + n] = 0.0;
    }
    return XTRACT_SUCCESS;
}

int xtract_lpc(const real_t *data, const int N, const void *argv, real_t *result)
{

    int i, j, M, L;
    real_t r = 0.0,
          error = 0.0;

    real_t *ref = NULL,
           *lpc = NULL ;

    error = data[0];
    L = N - 1; /* The number of LPC coefficients */
    M = L * 2; /* The length of *result */
    ref = result;
    lpc = result+L;

    if(error == 0.0)
    {
        memset(result, 0, M * sizeof(real_t));
        return XTRACT_NO_RESULT;
    }

    memset(result, 0, M * sizeof(real_t));

    for (i = 0; i < L; i++)
    {

        /* Sum up this iteration's reflection coefficient. */
        r = -data[i + 1];
        for (j = 0; j < i; j++)
            r -= lpc[j] * data[i - j];
        ref[i] = r /= error;

        /* Update LPC coefficients and total error. */
        lpc[i] = r;
        for (j = 0; j < i / 2; j++)
        {
            real_t tmp      = lpc[j];
            lpc[j]          = r * lpc[i - 1 - j];
            lpc[i - 1 - j] += r * tmp;
        }
        if (i % 2) lpc[j] += lpc[j] * r;

        error *= 1 - r * r;
    }

    return XTRACT_SUCCESS;
}

int xtract_lpcc(const real_t *data, const int N, const void *argv, real_t *result)
{

    /* Given N lpc coefficients extract an LPC cepstrum of size argv[0] */
    /* Based on an an algorithm by rabiner and Juang */

    int n, k;
    real_t sum;
    int order = N - 1; /* Eventually change this to Q = 3/2 p as suggested in Rabiner */
    int cep_length;

    if(argv == NULL)
        cep_length = N - 1; /* FIX: if we're going to have default values, they should come from the descriptor */
    else
        cep_length = *(int *)argv;
    //cep_length = (int)((real_t *)argv)[0];

    memset(result, 0, cep_length * sizeof(real_t));

    for (n = 1; n <= order && n <= cep_length; n++)
    {
        sum = 0.0;
        for (k = 1; k < n; k++)
            sum += k * result[k-1] * data[n - k];
        result[n-1] = data[n] + sum / n;
    }

    /* be wary of these interpolated values */
    for(n = order + 1; n <= cep_length; n++)
    {
        sum = 0.0;
        for (k = n - (order - 1); k < n; k++)
            sum += k * result[k-1] * data[n - k];
        result[n-1] = sum / n;
    }

    return XTRACT_SUCCESS;

}
//int xtract_lpcc_s(const real_t *data, const int N, const void *argv, real_t *result){
//    return XTRACT_SUCCESS;
//}

int xtract_subbands(const real_t *data, const int N, const void *argv, real_t *result)
{

    int n, bw, xtract_func, nbands, scale, start, lower, *argi, rv;

    argi = (int *)argv;

    xtract_func = argi[0];
    nbands = argi[1];
    scale = argi[2];
    start = argi[3];

    if(scale == XTRACT_LINEAR_SUBBANDS)
        bw = floorf((N - start) / nbands);
    else
        bw = start;

    lower = start;
    rv = XTRACT_SUCCESS;

    for(n = 0; n < nbands; n++)
    {

        /* Bounds sanity check */
        if(lower >= N || lower + bw >= N)
        {
            //   printf("n: %d\n", n);
            result[n] = 0.0;
            continue;
        }

        rv = xtract[xtract_func](data+lower, bw, NULL, &result[n]);

        if(rv != XTRACT_SUCCESS)
            return rv;

        switch(scale)
        {
        case XTRACT_OCTAVE_SUBBANDS:
            lower += bw;
            bw = lower;
            break;
        case XTRACT_LINEAR_SUBBANDS:
            lower += bw;
            break;
        }

    }

    return rv;

}



