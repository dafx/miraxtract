/*
 * Copyright (c) 2024, Yan Li
 * All rights reserved.
 */

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
