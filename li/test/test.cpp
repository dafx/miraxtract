#include "../examples/simpletest/WaveFile.h"

#include "../src/beat.hpp"

int main()
{
    WaveFile wavFile("test.wav");
    if (!wavFile.IsLoaded())
    {
        return EXIT_FAILURE;
    }

    float *wavData = (float *)wavFile.GetData(); // assume 32-bit float
    std::size_t wavBytes = wavFile.GetDataSize();
    uint64_t wavSamples = wavBytes / sizeof(float);
    li::fvec data(wavSamples);

    li::beat_det bd(wavFile.GetSampleRate());
    li::fmat in(1, bd.frame_size_full_rate);
    std::vector<float> flux(wavSamples / bd.frame_size_full_rate);
    for (size_t i = 0; i < wavSamples; i += bd.frame_size_full_rate)
    {
        in.copy(wavData + i, bd.frame_size_full_rate);
        bd.process(in);

        flux[i / bd.frame_size_full_rate] = bd.get_flux();
        printf("BPM: %.2f\n", bd.get_bpm());
    }
    return 0;
}
