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

#include <limits>
#include <cmath>
#include <cstdint>

#include "metadata.h"

#include "dbg.h"
#include "str.h"
#include "fio.h"


void metadata::set_quad_meta(const ecmdb& database,
        const blob& data, const blob& mask,
        ecmdb::quad_metadata* meta)
{
    meta->create(database.category());
    if (database.category() == ecmdb::category_elevation) {
        assert(database.channels() == 1);
        assert(database.type() == ecmdb::type_float32 || database.type() == ecmdb::type_int16);
        float min_elev = +std::numeric_limits<float>::max();
        float max_elev = -std::numeric_limits<float>::max();
        for (int i = 0; i < database.total_quad_size() * database.total_quad_size(); i++) {
            if (mask.ptr<uint8_t>()[i]) {
                float x;
                if (database.type() == ecmdb::type_float32) {
                    x = data.ptr<float>()[i];
                } else {
                    x = data.ptr<int16_t>()[i];
                }
                x = database.data_offset() + database.data_factor() * x;
                if (x < min_elev)
                    min_elev = x;
                if (x > max_elev)
                    max_elev = x;
            }
        }
        meta->elevation.min = min_elev;
        meta->elevation.max = max_elev;
    }
}

void metadata::update_global(const ecmdb& database,
        const ecmdb::global_metadata& quad_global_metadata,
        ecmdb::global_metadata& global_metadata)
{
    if (database.category() == ecmdb::category_elevation) {
        if (quad_global_metadata.elevation.min < global_metadata.elevation.min)
            global_metadata.elevation.min = quad_global_metadata.elevation.min;
        if (quad_global_metadata.elevation.max > global_metadata.elevation.max)
            global_metadata.elevation.max = quad_global_metadata.elevation.max;
    }
}
