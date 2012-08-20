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

#include <cstdio>
#include <algorithm>

#include "fio.h"
#include "str.h"
#include "intcheck.h"

#include "quadlist.h"


/* An entry represented as a string:
 * side, qx, qy, id; space-separated, fixed widths */
static const size_t record_size = 1 + 1 + 10 + 1 + 10 + 1 + 36 + 1;

static void entry_to_record(const quadlist::entry& e, char r[record_size + 1])
{
    std::snprintf(r, record_size + 1, "%d %010d %010d %36s\n", e.side, e.qx, e.qy, e.datafile_id.c_str());
}

static quadlist::entry entry_from_record(const char* r)
{
    quadlist::entry e;
    char id[37];
    if (std::sscanf(r, "%d %d %d %36s", &e.side, &e.qx, &e.qy, id) != 4) {
        // fill in obviously invalid values
        e.side = -1;
        e.qx = -1;
        e.qy = -1;
        e.datafile_id = "INVALID-INVALID-INVALID-INVALID-INVA";
    } else {
        e.datafile_id = id;
    }
    return e;
}


std::string quadlist::get_addfilename(const std::string& dir, const std::string& process_id)
{
    return dir + "/added." + process_id + ".txt";
}

FILE* quadlist::open_addfile(const std::string& addfilename)
{
    FILE* f = fio::open(addfilename, "a");
    fio::disable_buffering(f, addfilename);
    return f;
}

void quadlist::add_entry(const std::string& addfilename, FILE* addfile, const entry& e)
{
    char r[record_size + 1];
    entry_to_record(e, r);
    fio::write(r, record_size, 1, addfile, addfilename);
}

static void read_file(const std::string& filename, std::vector<quadlist::entry>& entries)
{
    FILE* f = fio::open(filename, "r");
    fio::seek(f, 0, SEEK_END, filename);
    off_t len = fio::tell(f, filename);
    fio::rewind(f, filename);
    char r[record_size + 1];
    r[record_size] = '\0';
    for (size_t i = 0; i < checked_cast<size_t>(len / record_size); i++) {
        fio::read(r, record_size, 1, f, filename);
        entries.push_back(entry_from_record(r));
    }
    fio::close(f);
}

std::vector<quadlist::entry> quadlist::read_addfiles(const std::string& dir)
{
    std::vector<entry> entries;
    std::vector<std::string> addfilenames = fio::readdir(dir, "added.*.txt");
    for (size_t i = 0; i < addfilenames.size(); i++) {
        read_file(dir + '/' + addfilenames[i], entries);
    }
    std::sort(entries.begin(), entries.end());
    return entries;
}

std::vector<quadlist::entry> quadlist::read_commitfile(const std::string& dir)
{
    std::vector<entry> entries;
    if (fio::test_f(dir + "/committed.txt")) {
        read_file(dir + "/committed.txt", entries);
    }
    std::sort(entries.begin(), entries.end());
    return entries;
}

std::vector<quadlist::entry> quadlist::get_uncommitted(const std::vector<entry>& added, const std::vector<entry>& committed)
{
    std::vector<entry> entries;
    for (size_t i = 0, j = 0; i < added.size(); i++) {
        if (j < committed.size() && added[i] == committed[j]) {
            j++;
            continue;
        } else {
            entries.push_back(added[i]);
        }
    }
    return entries;
}

void quadlist::write_commitfile(const std::string& dir, const std::vector<entry>& entries)
{
    std::string filename = dir + "/committed.txt";
    FILE* f = fio::open(filename, "w");
    char r[record_size + 1];
    for (size_t i = 0; i < entries.size(); i++) {
        entry_to_record(entries[i], r);
        fio::write(r, record_size, 1, f, filename);
    }
    fio::close(f, filename);
}
