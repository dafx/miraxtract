/*
 * Copyright (c) 2024, Yan Li
 * All rights reserved.
 */

#pragma once

#include "common.hpp"

namespace li
{
    struct sosfilt {
    };

    struct sosfilt_zi {
    };

    struct butter {
        static void lowpass(sosfilt &sos, sosfilt_zi &zi, double Wn, int N, double fs);
    };

    struct iir{
    };
}
