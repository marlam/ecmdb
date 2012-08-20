/*
 * Copyright (C) 2012
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

#include "fio.h"
#include "str.h"

#include "compression-info.h"


void set_compression_info(const std::string& dbdir, bool lossy_compression, int lossy_compression_quality)
{
    if (lossy_compression) {
        std::string lcname = dbdir + '/' + "lossy-compression.txt";
        FILE* lc = fio::open(lcname, "w");
        fprintf(lc, "%d\n", lossy_compression_quality);
        fio::close(lc, lcname);
    }
}

void get_compression_info(const std::string& dbdir, bool* lossy_compression, int* lossy_compression_quality)
{
    std::string lcname = dbdir + '/' + "lossy-compression.txt";
    if (fio::test_f(lcname)) {
        FILE* lc = fio::open(lcname, "r");
        *lossy_compression = true;
        *lossy_compression_quality = str::to<int>(fio::readline(lc, lcname));
        fio::close(lc, lcname);
    } else {
        *lossy_compression = false;
        *lossy_compression_quality = 100;
    }
}

void remove_compression_info(const std::string& dbdir)
{
    std::string lcname = dbdir + '/' + "lossy-compression.txt";
    if (fio::test_f(lcname)) {
        fio::unlink(lcname);
    }
}
