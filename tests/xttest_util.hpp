

#include <stdint.h>

#include "xtract/xtract_types.h"

// Fill table with sine wave at given frequency and amplitude
void xttest_gen_sine(real_t *table, uint32_t tablesize, real_t samplerate, real_t frequency, real_t amplitude);

// Fill table with sawtooth wave at given frequency and amplitude
void xttest_gen_sawtooth(real_t *table, uint32_t tablesize, real_t samplerate, real_t frequency, real_t amplitude);

// Fill table with noise at given frequency and amplitude
// N.B. The implementation actually provides "fake" noise from a table for reproducible testing
void xttest_gen_noise(real_t *table, uint32_t tablesize, real_t amplitude);

// Add table1 and table2 sample-by-sample leaving the result in table1
void xttest_add(real_t *table1, real_t *table2, uint32_t tablesize);

// Multiply table by a constant leavint the result in table
void xttest_mul(real_t *table, uint32_t tablesize, real_t constant);

// Return MIDI cent value for frequency
uint16_t xttest_ftom(real_t frequency);
