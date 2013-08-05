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
#include <cinttypes>

#include <ecmdb/ecmdb.h>

#include "exc.h"
#include "msg.h"
#include "opt.h"
#include "dbg.h"
#include "str.h"


extern "C" void ecmdb_filename_help(void)
{
    msg::req("filename QUAD_SIDE QUAD_LEVEL QUAD_X QUAD_Y\n"
             "filename -r|--reverse FILENAME\n"
            "\n"
            "Print the file name corresponding to the given quad. The file name is relative to a data base URL.\n"
            "When -r is given, the reverse mapping is performed: for a given file name, the "
            "corresponding quad is printed.");
}

extern "C" int ecmdb_filename(int argc, char* argv[])
{
    std::vector<opt::option *> options;
    opt::info help("help", '\0', opt::optional);
    options.push_back(&help);
    opt::flag reverse("reverse", 'r', opt::optional);
    options.push_back(&reverse);
    std::vector<std::string> arguments;
    if (!opt::parse(argc, argv, options, 1, 4, arguments)) {
        return 1;
    }
    if (help.value()) {
        ecmdb_filename_help();
        return 0;
    }
    if (reverse.value() && arguments.size() != 1) {
        msg::err("Need exactly one argument for -r");
        return 1;
    }
    if (!reverse.value() && arguments.size() != 4) {
        msg::err("Need exactly four arguments for -r");
        return 1;
    }

    try {
        if (reverse.value()) {
            std::string filename = arguments[0];
            int qs, ql, qx, qy;
            bool ok = false;
            size_t first_slash = filename.find('/');
            if (first_slash != std::string::npos
                    && std::sscanf(filename.substr(0, first_slash).c_str(), "%d-%d", &qs, &ql) == 2
                    && qs >= 0 && qs <= 5 && ql >= 0 && ql < ecmdb::max_levels)
            {
                std::string xy = filename.substr(first_slash + 1);
                if (xy.size() >= 4 && xy.substr(xy.size() - 4) == ".gta") {
                    std::string istr = xy.substr(0, xy.size() - 4);
                    str::replace(istr, "/", "");
                    uint_least64_t i;
                    if (std::sscanf(istr.c_str(), "%" SCNxLEAST64, &i) == 1) {
                        qx = i & ((1 << ql) - 1);
                        qy = i >> ql;
                        if (qx >= 0 && qx < (1 << ql) && qy >= 0 && qy < (1 << ql))
                            ok = true;
                    }
                }
            }
            if (ok) {
                msg::req("%d %d %d %d", qs, ql, qx, qy);
            } else {
                throw exc(std::string("invalid filename: ") + filename);
            }
        } else {
            int qs = str::to<int>(arguments[0]);
            int ql = str::to<int>(arguments[1]);
            int qx = str::to<int>(arguments[2]);
            int qy = str::to<int>(arguments[3]);
            if (qs < 0 || qs > 5
                    || ql < 0 || ql >= ecmdb::max_levels
                    || qx < 0 || qx >= (1 << ql)
                    || qy < 0 || qy >= (1 << ql)) {
                throw exc(str::asprintf("invalid quad %d %d %d %d", qs, ql, qx, qy));
            }
            std::string filename = ecmdb::quad_filename(qs, ql, qx, qy);
            msg::req(filename);
        }
    }
    catch (std::exception& e) {
        msg::err_txt("%s", e.what());
        return 1;
    }

    return 0;
}
