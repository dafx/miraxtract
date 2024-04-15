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

#include "signal.hpp"

using namespace li;

const float deci_4_filter_coefs[] = {
    /*00*/ 2.395964e-05f,
    /*01*/ 4.791929e-05f,
    /*02*/ 2.395964e-05f,
    /*03*/ 1.000000e+00f,
    /*04*/ -1.026351e+00f,
    /*05*/ 2.686402e-01f,
    /*06*/ 1.000000e+00f,
    /*07*/ 2.000000e+00f,
    /*08*/ 1.000000e+00f,
    /*09*/ 1.000000e+00f,
    /*10*/ -1.086858e+00f,
    /*11*/ 3.434309e-01f,
    /*12*/ 1.000000e+00f,
    /*13*/ 2.000000e+00f,
    /*14*/ 1.000000e+00f,
    /*15*/ 1.000000e+00f,
    /*16*/ -1.219725e+00f,
    /*17*/ 5.076635e-01f,
    /*18*/ 1.000000e+00f,
    /*19*/ 2.000000e+00f,
    /*20*/ 1.000000e+00f,
    /*21*/ 1.000000e+00f,
    /*22*/ -1.451580e+00f,
    /*23*/ 7.942511e-01f};

const float deci_4_filter_zi[] = {
    /*00*/ 3.715956e-04f,
    /*01*/ -8.230240e-05f,
    /*02*/ 5.771205e-03f,
    /*03*/ -1.722301e-03f,
    /*04*/ 7.950110e-02f,
    /*05*/ -3.732368e-02f,
    /*06*/ 9.143321e-01f,
    /*07*/ -7.085832e-01f};

deci_by_4::deci_by_4()
    : sos(deci_4_filter_coefs, deci_4_filter_zi, 4)
{
}

void deci_by_4::run(float *input, const size_t size, float *output)
{
    if(init_done == false) {
        sos.init(input[0]);
        init_done = true;
    }

    sos.run(input, size, input);

    for (size_t i = 0; i < size; i += 4)
    {
        output[i / 4] = input[i];
    }
}
