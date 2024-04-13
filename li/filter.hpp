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
