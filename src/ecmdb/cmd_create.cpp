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

#include <ecmdb/ecmdb.h>

#include "exc.h"
#include "dbg.h"
#include "fio.h"
#include "msg.h"
#include "opt.h"

#include "compression-info.h"


extern "C" void ecmdb_create_help(void)
{
    msg::req("create OPTIONS... DIR\n"
            "[-E|--ellipsoid-name=NAME] Name of a well-known ellipsoid to use.\n"
            "                           Currently earth-wgs84, moon-nasa, or mars-nasa.\n"
            "[-e|--ellipsoid=a,b]       Semi-major and semi-minor axis of the ellipsoid.\n"
            "                           Overrides -E. Default: earth-wgs84.\n"
            "[-q|--quad-size=q]         Quad size (q x q). Default 512.\n"
            "-l|--levels=l              Number of levels in the new database.\n"
            "-C|--category=C            Category: elevation, texture, sar-amplitude, or data (= undefined).\n"
            "[-t|--type=t]              Data type: uint8, int16, or float32. Default depends on category.\n"
            "[-c|--channels=c]          Number of channels. Default depends on category.\n"
            "[-o|--overlap=o]           Quad overlap. Default 2 for elevation, 1 otherwise.\n"
            "[-O|--offset=O]            Data offset (v = F * v + O). Default 0.\n"
            "[-F|--factor=F]            Data factor (v = F * v + O). Default 1.\n"
            "[-l|--lossy[=on|off]]      Lossy compression. Enabled by default.\n"
            "[--lossy-quality=Q]        Quality of lossy compression, from 1 to 100, default 80.\n"
            "[-S|--short-description=S] Short description (one line). Default: empty.\n"
            "[-D|--description=D]       Description text. Default: short description.\n"
            "\n"
            "Create a new database in the given directory DIR, which will be created.");
}

extern "C" int ecmdb_create(int argc, char* argv[])
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
    opt::val<int> levels("levels", 'l', opt::required, 1, ecmdb::max_levels);
    options.push_back(&levels);
    std::vector<std::string> category_names;
    category_names.push_back("elevation");
    category_names.push_back("texture");
    category_names.push_back("data");
    category_names.push_back("sar-amplitude");
    opt::string category("category", 'C', opt::required, category_names);
    options.push_back(&category);
    std::vector<std::string> type_names;
    type_names.push_back("uint8");
    type_names.push_back("int16");
    type_names.push_back("float32");
    opt::string type("type", 't', opt::optional, type_names);
    options.push_back(&type);
    opt::val<int> channels("channels", 'c', opt::optional, 1, ecmdb::max_channels, 0);
    options.push_back(&channels);
    opt::val<int> overlap("overlap", 'o', opt::optional, 1, ecmdb::max_overlap, 1);
    options.push_back(&overlap);
    opt::val<float> offset("offset", 'O', opt::optional, 0.0f);
    options.push_back(&offset);
    opt::val<float> factor("factor", 'F', opt::optional, 1.0f);
    options.push_back(&factor);
    opt::flag lossy("lossy", 'l', opt::optional, true, true);
    options.push_back(&lossy);
    opt::val<int> lossy_quality("lossy-quality", '\0', opt::optional, 1, 100, 80);
    options.push_back(&lossy_quality);
    opt::string short_description("short-description", 'S', opt::optional);
    options.push_back(&short_description);
    opt::string description("description", 'D', opt::optional, "\n", "");
    options.push_back(&description);
    std::vector<std::string> arguments;
    if (!opt::parse(argc, argv, options, 1, 1, arguments)) {
        return 1;
    }
    if (help.value()) {
        ecmdb_create_help();
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

    // Fill in dependent defaults
    ecmdb::category_t db_category =
        (category.value().compare("elevation") == 0 ? ecmdb::category_elevation
         : category.value().compare("texture") == 0 ? ecmdb::category_texture
         : category.value().compare("sar-amplitude") == 0 ? ecmdb::category_sar_amplitude
         : ecmdb::category_data);
    ecmdb::type_t db_type;
    if (type.values().empty()) {
        switch (db_category) {
        case ecmdb::category_elevation:
            db_type = ecmdb::type_float32;
            break;
        case ecmdb::category_texture:
            db_type = ecmdb::type_uint8;
            break;
        case ecmdb::category_data:
            db_type = ecmdb::type_float32;
            break;
        case ecmdb::category_sar_amplitude:
            db_type = ecmdb::type_float32;
            break;
        default:
            db_type = ecmdb::type_float32;
            break;
        }
    } else {
        db_type = 
            (type.value().compare("uint8") == 0 ? ecmdb::type_uint8
             : type.value().compare("int16") == 0 ? ecmdb::type_int16
             : ecmdb::type_float32);
    }
    int db_channels;
    if (channels.values().empty()) {
        switch (db_category) {
        case ecmdb::category_elevation:
            db_channels = 1;
            break;
        case ecmdb::category_texture:
            db_channels = 3;
            break;
        case ecmdb::category_data:
            db_channels = 1;
            break;
        case ecmdb::category_sar_amplitude:
            db_channels = 1;
            break;
        default:
            db_channels = 1;
            break;
        }
    } else {
        db_channels = channels.value();
    }
    int db_overlap;
    if (overlap.values().empty()) {
        switch (db_category) {
        case ecmdb::category_elevation:
            db_overlap = 2;
            break;
        case ecmdb::category_sar_amplitude:
            db_overlap = 9;
            break;
        default:
            db_overlap = 1;
            break;
        }
    } else {
        db_overlap = overlap.value();
    }
    std::vector<std::string> db_description;
    if (description.values().empty()) {
        if (!short_description.value().empty()) {
            db_description.push_back(short_description.value());
        }
    } else {
        size_t last_line_start = 0;
        for (size_t i = 0; i < description.value().size(); i++) {
            if (description.value()[i] == '\n') {
                db_description.push_back(description.value().substr(last_line_start, i - last_line_start));
                last_line_start = i;
            }
        }
        if (description.value().size() > 0 && description.value()[description.value().size() - 1] != '\n') {
            db_description.push_back(description.value().substr(last_line_start,
                        description.value().size() - last_line_start));
        }
    }

    // Create database
    try {
        fio::mkdir(arguments[0]);
        ecmdb database;
        database.create(
                semi_major_axis, semi_minor_axis, quad_size.value(), levels.value(),
                db_category, db_type, db_channels, db_overlap, offset.value(), factor.value(),
                short_description.value(), db_description);
        database.write(arguments[0]);
        ecmdb::metadata metadata;
        metadata.create(db_category);
        metadata.write(arguments[0]);
        set_compression_info(arguments[0], lossy.value(), lossy_quality.value());
    }
    catch (std::exception& e) {
        msg::err_txt("%s", e.what());
        return 1;
    }

    return 0;
}
