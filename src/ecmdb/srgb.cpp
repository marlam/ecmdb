/*
 * Copyright (C) 2011, 2012
 * Computer Graphics Group, University of Siegen, Germany.
 * Written by Martin Lambers <martin.lambers@uni-siegen.de>.
 * See http://www.cg.informatik.uni-siegen.de/ for contact information.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"

#include <cmath>

#include "srgb.h"

#include "dbg.h"


uint16_t srgb::lut_nonlinear_to_linear[256];
uint8_t srgb::lut_linear_to_nonlinear[2048];

void srgb::initialize_luts()
{
    // See e.g. GL_ARB_framebuffer_sRGB extension for the formulas
    for (int i = 0; i < 2048; i++) {
        float x = i / 2047.0f;
        float n = (x <= 0.04045f ? x / 12.92f : std::pow((x + 0.055f) / 1.055f, 2.4f));
        lut_linear_to_nonlinear[i] = std::round(n * 255.0);
    }
    for (int i = 0; i < 256; i++) {
        float x = i / 255.0f;
        float l = (x <= 0.0031308f ? (x * 12.92f) : (1.055f * std::pow(x, 1.0f / 2.4f) - 0.055f));
        lut_nonlinear_to_linear[i] = std::round(l * 2047.0f);
        assert(srgb::lut_linear_to_nonlinear[lut_nonlinear_to_linear[i]] == i);
    }
}
