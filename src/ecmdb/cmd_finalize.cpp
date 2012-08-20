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

#include <string>

#include <ecm/ecm.h>
#include <ecm/ecmdb.h>

#include "msg.h"
#include "str.h"
#include "opt.h"
#include "fio.h"
#include "exc.h"
#include "thread.h"

#include "quadlist.h"
#include "compression-info.h"


extern "C" void ecmdb_finalize_help(void)
{
    msg::req_txt("finalize [--force] DIR\n"
            "\n"
            "Remove all unnecessary files from previous 'add' and 'commit' commands from DIR.\n"
            "Once you run this command, you cannot add or commit more quads to the database in DIR.");
}

extern "C" int ecmdb_finalize(int argc, char* argv[])
{
    std::vector<opt::option *> options;
    opt::info help("help", '\0', opt::optional);
    options.push_back(&help);
    std::vector<std::string> arguments;
    if (!opt::parse(argc, argv, options, 1, 1, arguments)) {
        return 1;
    }
    if (help.value()) {
        ecmdb_finalize_help();
        return 0;
    }

    std::string dir = arguments[0];
    std::string infofilename = dir + "/ecmdb.txt";
    FILE *infofile = NULL;
    try {
        exc e;
        bool abort = false;
        ecmdb database;
        database.open(dir);
        FILE *infofile = fio::open(infofilename, "r+");
        if (!fio::writelock(infofile, infofilename)) {
            throw exc(dir + ": database is locked, probably by a commit or finalize command.");
        }
        std::vector<quadlist::entry> added_quads = quadlist::read_addfiles(dir);
        std::vector<quadlist::entry> committed_quads = quadlist::read_commitfile(dir);
        std::vector<quadlist::entry> uncommitted_quads = quadlist::get_uncommitted(added_quads, committed_quads);
        if (uncommitted_quads.size() > 0) {
            throw exc(dir + ": contains uncommitted quads");
        } else {
            msg::inf(dir + ": removing " + str::from(added_quads.size()) + " quad files");
            size_t added_quads_index = 0;
            #pragma omp parallel for
            for (size_t i = 0; i < added_quads.size(); i++) {
                #pragma omp flush (abort)
                if (!abort) {
                    try {
                        size_t j = atomic::fetch_and_inc(&added_quads_index);
                        std::string filename = dir + '/'
                            + ecmdb::quad_filename(added_quads[j].side, database.levels() - 1, added_quads[j].qx, added_quads[j].qy);
                        filename = filename.substr(0, filename.length() - 4) // remove ".gta"
                            + '.' + added_quads[j].datafile_id + ".gta";
                        fio::unlink(filename);
                    }
                    catch (std::exception& _e) {
                        #pragma omp critical
                        e = _e;
                        abort = true;
                        #pragma omp flush (abort)
                    }
                }
            }
            if (!e.empty()) {
                throw e;
            }
            std::vector<std::string> rm_filenames = fio::readdir(dir, "added.*.txt");
            if (fio::test_f(dir + "/committed.txt")) {
                rm_filenames.push_back("committed.txt");
            }
            msg::inf(dir + ": removing " + str::from(rm_filenames.size()) + " quad list files");
            #pragma omp parallel for
            for (size_t i = 0; i < rm_filenames.size(); i++) {
                #pragma omp flush (abort)
                if (!abort) {
                    try {
                        fio::unlink(dir + '/' + rm_filenames[i]);
                    }
                    catch (std::exception& _e) {
                        #pragma omp critical
                        e = _e;
                        abort = true;
                        #pragma omp flush (abort)
                    }
                }
            }
            if (!e.empty()) {
                throw e;
            }
            remove_compression_info(dir);
        }
        try { fio::close(infofile, infofilename); } catch (...) {}
        infofile = NULL;
    }
    catch (std::exception& e) {
        if (infofile) {
            try { fio::close(infofile, infofilename); } catch (...) {}
        }
        msg::err_txt("%s", e.what());
        return 1;
    }

    return 0;
}
