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

#include <ecmdb/ecmdb.h>

#include "exc.h"
#include "msg.h"
#include "opt.h"
#include "dbg.h"


extern "C" void ecmdb_resolutions_help(void)
{
    msg::req("resolutions [OPTIONS...]\n"
            "[-E|--ellipsoid-name=NAME] Name of a well-known ellipsoid to use.\n"
            "                           Currently earth-wgs84, moon-nasa, or mars-nasa.\n"
            "[-e|--ellipsoid=a,b]       Semi-major and semi-minor axis of the ellipsoid.\n"
            "                           Overrides -E. Default: earth-wgs84.\n"
            "[-q|--quad-size=q]         Quad size (q x q). Default 512.\n"
            "\n"
            "List the ground resolution of one quad sample for all levels (0-%d),\n"
            "for the given ellipsoid and quad configuration.\n"
            "This list can be used to choose a suitable number of levels for\n"
            "a new database.", ecmdb::max_levels);
}

extern "C" int ecmdb_resolutions(int argc, char* argv[])
{
    std::vector<opt::option *> options;
    opt::info help("help", '\0', opt::optional);
    options.push_back(&help);
    std::vector<std::string> ellipsoid_names;
    ellipsoid_names.push_back("earth-wgs84");
    ellipsoid_names.push_back("moon-nasa");
    ellipsoid_names.push_back("mars-nasa");
    opt::val<std::string> ellipsoid_name("ellipsoid-name", 'E', opt::optional, ellipsoid_names, "");
    options.push_back(&ellipsoid_name);
    opt::tuple<double> ellipsoid("ellipsoid", 'e', opt::optional, 0.0, false, std::numeric_limits<double>::max(), true, std::vector<double>(), 2);
    options.push_back(&ellipsoid);
    opt::val<int> quad_size("quad-size", 'q', opt::optional, 1, ecmdb::max_quad_size, 512);
    options.push_back(&quad_size);
    std::vector<std::string> arguments;
    if (!opt::parse(argc, argv, options, 0, 0, arguments)) {
        return 1;
    }
    if (help.value()) {
        ecmdb_resolutions_help();
        return 0;
    }

    double semi_major_axis, semi_minor_axis;
    if (!ellipsoid.values().empty()) {
        semi_major_axis = ellipsoid.value()[0];
        semi_minor_axis = ellipsoid.value()[1];
    } else if (!ellipsoid_name.values().empty()) {
        if (ellipsoid_name.value() == "earth-wgs84") {
            semi_major_axis = ecm::semi_major_axis_earth_wgs84;
            semi_minor_axis = ecm::semi_minor_axis_earth_wgs84;
        } else if (ellipsoid_name.value() == "moon-nasa") {
            semi_major_axis = ecm::radius_moon_nasa;
            semi_minor_axis = ecm::radius_moon_nasa;
        } else {
            assert(ellipsoid_name.value() == "mars-nasa");
            semi_major_axis = ecm::semi_major_axis_mars_nasa;
            semi_minor_axis = ecm::semi_minor_axis_mars_nasa;
        }
    } else {
        semi_major_axis = ecm::semi_major_axis_earth_wgs84;
        semi_minor_axis = ecm::semi_minor_axis_earth_wgs84;
    }

    (void)semi_minor_axis;         // Unused; avoid compiler warning.
    for (int l = 0; l < ecmdb::max_levels; l++) {
        double quads = (1 << l);
        quads *= 4.0;  // front, right, back, left
        double pixels = quads * quad_size.value();
        double r = 2.0 * M_PI * semi_major_axis / pixels;
        double p = pixels / 360.0;
        if (r >= 10000.0) {
            r /= 1000.0;
            msg::req("Level %2d: ca. %6.4g km pixel spacing; ca. %6.4f pixels per degree at equator", l, r, p);
        } else if (r < 0.01) {
            r *= 1000.0;
            msg::req("Level %2d: ca. %6.4g mm pixel spacing; ca. %6.4f pixels per degree at equator", l, r, p);
        } else if (r < 1.0) {
            r *= 100.0;
            msg::req("Level %2d: ca. %6.4g cm pixel spacing; ca. %6.4f pixels per degree at equator", l, r, p);
        } else {
            msg::req("Level %2d: ca. %6.4g  m pixel spacing; ca. %6.4f pixels per degree at equator", l, r, p);
        }
    }

    return 0;
}
