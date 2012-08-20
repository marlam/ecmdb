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

#ifndef QUADLIST_H

#include <cstdio>
#include <string>
#include <vector>

namespace quadlist
{
    class entry
    {
    public:
        std::string datafile_id;
        int side;
        int qx;
        int qy;
        
        entry()
        {
        }
        entry(const std::string& id, int s, int x, int y) :
            datafile_id(id), side(s), qx(x), qy(y)
        {
        }

        bool operator<(const entry& e) const {
            return (side < e.side
                    || (side == e.side && qx < e.qx)
                    || (side == e.side && qx == e.qx && qy < e.qy)
                    || (side == e.side && qx == e.qx && qy == e.qy && datafile_id.compare(e.datafile_id) < 0));
        }

        bool operator==(const entry& e) const {
            return (side == e.side
                    && qx == e.qx
                    && qy == e.qy
                    && (datafile_id.compare(e.datafile_id) == 0));
        }
    };

    /* For the 'add' command:
     * - Open DIR/added.process_id.txt in append mode
     * - Add entries for added quads
     */
    std::string get_addfilename(const std::string& dir, const std::string& process_id);
    FILE* open_addfile(const std::string& addfilename);
    void add_entry(const std::string& addfilename, FILE* addfile, const entry& e);

    /* For the 'commit' command:
     * - Read the current entries from all files DIR/added.*.txt (sorted)
     * - Read the current entries from DIR/committed.txt (sorted)
     * - Write new entries to DIR/committed.txt
     */
    std::vector<entry> read_addfiles(const std::string& dir);
    std::vector<entry> read_commitfile(const std::string& dir);
    std::vector<entry> get_uncommitted(const std::vector<entry>& added, const std::vector<entry>& committed);
    void write_commitfile(const std::string& dir, const std::vector<entry>& entries);
    
    /* For the 'finalize' command:
     * ... TODO
     */
}

#endif
