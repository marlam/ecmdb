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

#ifndef SRGB_H
#define SRGB_H

#include <cstdint>

namespace srgb
{
    extern uint16_t lut_nonlinear_to_linear[256];
    extern uint8_t lut_linear_to_nonlinear[2048];

    void initialize_luts();

    inline uint16_t nonlinear_to_linear(uint8_t linear)
    {
        return lut_nonlinear_to_linear[linear];
    }

    inline uint8_t linear_to_nonlinear(uint16_t nonlinear)
    {
        // keep value in range!!
        return lut_linear_to_nonlinear[nonlinear];
    }
}

#endif
