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
    struct phase_vocover {
        phase_vocover()
        {
        }

        void process(const fvec &in)
        {
            assert(in.size() == 512);
        }

        fvec window = {512};
        fvec input = {512};
        fvec temp = {512};
        fmat history = {1, 43};
    };

}
