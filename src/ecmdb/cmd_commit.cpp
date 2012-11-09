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

#include <vector>
#include <set>
#include <string>
#include <memory>
#include <algorithm>
#include <cstring>
#include <cmath>

#include <ecmdb/ecmdb.h>

#include "str.h"
#include "msg.h"
#include "opt.h"
#include "fio.h"
#include "exc.h"
#include "intcheck.h"
#include "blob.h"
#include "sys.h"
#include "lru-cache.h"

#include "quadlist.h"
#include "srgb.h"
#include "metadata.h"
#include "compression-info.h"


extern "C" void ecmdb_commit_help(void)
{
    msg::req_txt("commit DIR\n"
            "\n"
            "Commit all changes that were made with 'add' commands to the database stored in DIR.\n"
            "You can still run more 'add' commands, even in parallel to "
            "the 'commit' command, but only one 'commit' command can run "
            "for a database at any given time.\n"
            "Once you know that you will not require additional 'commit' "
            "commands, you can run the 'finalize' command to remove "
            "unnecessary data files.");
}

static void commit(const std::string& dir, bool lossy_compression, int lossy_compression_quality,
        ecmdb& database, ecmdb::metadata& global_metadata,
        const std::vector<quadlist::entry>& added_quads,
        const std::vector<quadlist::entry>& uncommitted_quads,
        std::vector<quadlist::entry>& committed_quads);

extern "C" int ecmdb_commit(int argc, char* argv[])
{
    std::vector<opt::option *> options;
    opt::info help("help", '\0', opt::optional);
    options.push_back(&help);
    std::vector<std::string> arguments;
    if (!opt::parse(argc, argv, options, 1, 1, arguments)) {
        return 1;
    }
    if (help.value()) {
        ecmdb_commit_help();
        return 0;
    }

    std::string dir = arguments[0];
    std::string infofilename = dir + "/ecmdb.txt";
    FILE *infofile = NULL;
    try {
        ecmdb database;
        database.open(dir);
        FILE *infofile = fio::open(infofilename, "r+");
        if (!fio::writelock(infofile, infofilename)) {
            throw exc(dir + ": database is locked, probably by another commit command.");
        }
        bool lossy_compression;
        int lossy_compression_quality;
        get_compression_info(arguments[0], &lossy_compression, &lossy_compression_quality);
        ecmdb::metadata global_metadata;
        global_metadata.open(database.category(), dir);
        std::vector<quadlist::entry> added_quads = quadlist::read_addfiles(dir);
        std::vector<quadlist::entry> committed_quads = quadlist::read_commitfile(dir);
        std::vector<quadlist::entry> uncommitted_quads = quadlist::get_uncommitted(added_quads, committed_quads);
        if (uncommitted_quads.size() == 0) {
            msg::inf(dir + ": nothing to commit");
        } else {
            if (database.category() == ecmdb::category_texture) {
                srgb::initialize_luts();
            }
            commit(dir, lossy_compression, lossy_compression_quality,
                    database, global_metadata, added_quads, uncommitted_quads, committed_quads);
        }
        quadlist::write_commitfile(dir, committed_quads);
        database.write(dir);
        global_metadata.write(dir);
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

// A quad cache for quads from the same database and with constant level
// (i.e. only side, qx, and qy vary). This is used to build lower level
// quads from higher level quads, since there is a lot of overlap between
// the used quads.

class quad_location
{
public:
    int side, qx, qy;
    quad_location(int side, int qx, int qy) : side(side), qx(qx), qy(qy) {}
    bool operator<(const quad_location& ql) const {
        return (side < ql.side || (side == ql.side && (qx < ql.qx || (qx == ql.qx && qy < ql.qy))));
    }
};

class quad_data
{
public:
    std::unique_ptr<blob> data;
    std::unique_ptr<blob> mask;
    quad_data(blob* d, blob* m) : data(d), mask(m) {}
};

class quad_cache : public lru_cache<quad_data, quad_location>
{
public:
    quad_cache() : lru_cache<quad_data, quad_location>(std::min(sys::total_ram() / 8,
                static_cast<uintmax_t>(std::numeric_limits<size_t>::max() / 2)))
    {
    }
};

class quad_load_mutexes
{
private:
    mutex _mutex;
    std::map<quad_location, mutex*> _mutex_map;

public:
    quad_load_mutexes()
    {
    }

    ~quad_load_mutexes()
    {
        for (auto it = _mutex_map.begin(); it != _mutex_map.end(); it++)
            delete it->second;
    }

    mutex* get(const quad_location& ql)
    {
        _mutex.lock();
        mutex* m;
        auto it = _mutex_map.find(ql);
        if (it == _mutex_map.end()) {
            m = new mutex;
            _mutex_map.insert(std::pair<quad_location, mutex*>(ql, m));
        } else {
            m = it->second;
        }
        _mutex.unlock();
        return m;
    }
};

void load_quad_data(const std::string& dir, const ecmdb& database, int qs, int ql, int qx, int qy, quad_data** element, size_t* size)
{
    std::string filename = dir + "/" + ecmdb::quad_filename(qs, ql, qx, qy);
    if (fio::test_f(filename)) {
        std::unique_ptr<blob> data(new blob());
        data->resize(database.data_size());
        std::unique_ptr<blob> mask(new blob());
        mask->resize(database.mask_size());
        bool all_valid;
        ecmdb::metadata meta;
        database.load_quad(filename, data->ptr(), mask->ptr<uint8_t>(), &all_valid, &meta);
        if (all_valid) {
            std::memset(mask->ptr(), 0xff, database.mask_size());
        }
        *size = sizeof(quad_data) + database.data_size() + database.mask_size();
        *element = new quad_data(data.release(), mask.release());
    } else {
        *size = sizeof(quad_data);
        *element = new quad_data(new blob(), new blob());
    }
}


// Commit one quad of the highest level:
static ecmdb::metadata commit_quad(const std::string& dir, const ecmdb& database,
        bool lossy_compression, int lossy_compression_quality,
        int side, int qx, int qy, const std::vector<std::string>& ids);

// Commit one quad of a lower level:
static void commit_ll_quad(const std::string& dir, const ecmdb& database,
        bool lossy_compression, int lossy_compression_quality,
        quad_cache& qc, quad_load_mutexes& qlm,
        int side, int level, int qx, int qy);

void commit(const std::string& dir, bool lossy_compression, int lossy_compression_quality,
        ecmdb& database, ecmdb::metadata& global_metadata,
        const std::vector<quadlist::entry>& added_quads,
        const std::vector<quadlist::entry>& uncommitted_quads,
        std::vector<quadlist::entry>& committed_quads)
{
    exc e;
    bool abort = false;
    std::vector<size_t> uncommitted_indices;
    std::vector<size_t> added_indices;

    // for each unique side,qx,qy triple in uncommitted_quads,
    // add an entry to uncommitted_index pointing to the first matching entry
    // in the uncommitted_quads list, and add an entry to added_index pointing
    // to the first matching entry in the added_quads list.
    {
        size_t u_index = 0;
        size_t a_index = 0;
        do
        {
            int u_side = uncommitted_quads[u_index].side;
            int u_qx = uncommitted_quads[u_index].qx;
            int u_qy = uncommitted_quads[u_index].qy;
            uncommitted_indices.push_back(u_index);
            while (added_quads[a_index].side < u_side
                    || (added_quads[a_index].side == u_side && added_quads[a_index].qx < u_qx)
                    || (added_quads[a_index].side == u_side && added_quads[a_index].qx == u_qx && added_quads[a_index].qy < u_qy)) {
                a_index++;
            }
            assert(added_quads[a_index].side == u_side
                    && added_quads[a_index].qx == u_qx
                    && added_quads[a_index].qy == u_qy);
            added_indices.push_back(a_index);
            while (u_index < uncommitted_quads.size()
                    && uncommitted_quads[u_index].side == u_side
                    && uncommitted_quads[u_index].qx == u_qx
                    && uncommitted_quads[u_index].qy == u_qy) {
                u_index++;
            }
        }
        while (u_index < uncommitted_quads.size());
    }

    std::vector<quadlist::entry> new_committed_quads;
    const int level = database.levels() - 1;
    msg::inf(dir + " level " + str::from(level) + ": committing " + str::from(uncommitted_indices.size()) + " quads");
    size_t uncommitted_indices_index = 0;
    #pragma omp parallel for
    for (size_t i = 0; i < uncommitted_indices.size(); i++) {
        size_t j = atomic::fetch_and_inc(&uncommitted_indices_index);
        #pragma omp flush (abort)
        if (!abort) {
            try {
                const size_t u_index = uncommitted_indices[j];
                const size_t a_index = added_indices[j];
                const int side = uncommitted_quads[u_index].side;
                const int qx = uncommitted_quads[u_index].qx;
                const int qy = uncommitted_quads[u_index].qy;
                std::vector<std::string> ids;
                for (size_t a = a_index;
                        a < added_quads.size()
                        && added_quads[a].side == side
                        && added_quads[a].qx == qx
                        && added_quads[a].qy == qy; a++) {
                    ids.push_back(added_quads[a].datafile_id);
                }
                const ecmdb::metadata quad_metadata = commit_quad(dir, database,
                        lossy_compression, lossy_compression_quality, side, qx, qy, ids);
                #pragma omp critical
                {
                    database.add_quad(side, qx, qy);
                    metadata::update_global(database, quad_metadata, global_metadata);
                    for (size_t u = u_index;
                            u < uncommitted_quads.size()
                            && uncommitted_quads[u].side == side
                            && uncommitted_quads[u].qx == qx
                            && uncommitted_quads[u].qy == qy; u++) {
                        new_committed_quads.push_back(uncommitted_quads[u]);
                    }
                }
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
    std::sort(new_committed_quads.begin(), new_committed_quads.end());
    committed_quads.insert(committed_quads.end(), new_committed_quads.begin(), new_committed_quads.end());

    // lower levels
    if (database.levels() > 1) {
        std::set<quad_location> ll_quads_set;
        for (size_t i = 0; i < new_committed_quads.size(); i++) {
            ll_quads_set.insert(quad_location(
                        new_committed_quads[i].side,
                        new_committed_quads[i].qx / 2,
                        new_committed_quads[i].qy / 2));
        }
        std::vector<quad_location> ll_quads(ll_quads_set.begin(), ll_quads_set.end());
        ll_quads_set.clear();
        for (int l = database.levels() - 2; l >= 0; l--) {
            msg::inf(dir + " level " + str::from(l) + ": committing " + str::from(ll_quads.size()) + " quads");
            quad_cache qc;
            quad_load_mutexes qlm;
            // The quads in ll_quads are sorted according to side,x,y. Each quad requires 16 higher level
            // quads, and must get them from the cache or load them from disk. We commit the quads in a
            // special order to 1) increase cache coherency and 2) reduce races to load the same quads from disk.
            size_t ll_quads_index = 0;
            #pragma omp parallel for
            for (size_t i = 0; i < ll_quads.size(); i++) {
                #pragma omp flush (abort)
                if (!abort) {
                    try {
                        size_t j = atomic::fetch_and_inc(&ll_quads_index);
                        if (j % 2 == 1) {
                            if (j + 1 < ll_quads.size())
                                j++;
                        } else if (j > 0) {
                            j--;
                        }
                        commit_ll_quad(dir, database, lossy_compression, lossy_compression_quality,
                                qc, qlm, ll_quads[j].side, l, ll_quads[j].qx, ll_quads[j].qy);
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
            // compute quads for next level
            if (l > 0) {
                for (size_t i = 0; i < ll_quads.size(); i++) {
                    ll_quads_set.insert(quad_location(
                                ll_quads[i].side,
                                ll_quads[i].qx / 2,
                                ll_quads[i].qy / 2));
                }
                ll_quads.assign(ll_quads_set.begin(), ll_quads_set.end());
                ll_quads_set.clear();
            }
        }
    }
}

ecmdb::metadata commit_quad(const std::string& dir, const ecmdb& database,
        bool lossy_compression, int lossy_compression_quality,
        int qs, int qx, int qy, const std::vector<std::string>& ids)
{
    int ql = database.levels() - 1;
    const std::string quad_filename = dir + '/' + ecmdb::quad_filename(qs, ql, qx, qy);
    const std::string base_filename = quad_filename.substr(0, quad_filename.length() - 4);      // remove ".gta"

    if (ids.size() == 1) {
        // Common case optimization: only one quad: simply make a hard link
        std::string added_filename = base_filename + '.' + ids[0] + ".gta";
        ecmdb::metadata quad_metadata;
        database.load_quad_meta(added_filename, &quad_metadata);
        try {
            fio::link(added_filename, quad_filename);
        } catch (exc &e) {
            // delete a preexisting file if necessary
            if (e.sys_errno() == EEXIST) {
                fio::unlink(quad_filename);
                fio::link(added_filename, quad_filename);
            } else {
                throw e;
            }
        }
        return quad_metadata;
    }

    blob data(database.data_size());
    std::memset(data.ptr(), 0, database.data_size());
    blob mask(database.mask_size());
    std::memset(mask.ptr(), 0, database.mask_size());

    blob added_data[ids.size()];
    blob added_mask[ids.size()];
    for (size_t i = 0; i < ids.size(); i++) {
        std::string added_filename = base_filename + '.' + ids[i] + ".gta";
        added_data[i].resize(database.data_size());
        added_mask[i].resize(database.mask_size());
        bool all_valid;
        ecmdb::metadata quad_metadata;
        database.load_quad(added_filename, added_data[i].ptr(), added_mask[i].ptr<uint8_t>(), &all_valid, &quad_metadata);
        if (all_valid)
            std::memset(added_mask[i].ptr(), 0xff, database.mask_size());
    }

    bool all_valid = true;
    bool none_valid = true;
    float v[database.channels()];
    for (int e = 0; e < database.total_quad_size() * database.total_quad_size(); e++) {
        for (int c = 0; c < database.channels(); c++) {
            v[c] = 0.0f;
        }
        int valid_values = 0;
        for (size_t i = 0; i < ids.size(); i++) {
            uint8_t m = added_mask[i].ptr<uint8_t>()[e];
            if (m > 0x00) {
                valid_values++;
                void *element = added_data[i].ptr(e * database.element_size());
                for (int c = 0; c < database.channels(); c++) {
                    if (database.category() == ecmdb::category_texture) {
                        v[c] += srgb::nonlinear_to_linear(static_cast<uint8_t *>(element)[c]);
                    } else if (database.type() == ecmdb::type_uint8) {
                        v[c] += static_cast<uint8_t *>(element)[c];
                    } else if (database.type() == ecmdb::type_int16) {
                        v[c] += static_cast<int16_t *>(element)[c];
                    } else {
                        v[c] += static_cast<float *>(element)[c];
                    }
                }
            }
        }
        if (valid_values > 0) {
            mask.ptr<uint8_t>()[e] = 0xff;
            void *element = data.ptr(e * database.element_size());
            for (int c = 0; c < database.channels(); c++) {
                v[c] /= valid_values;
                if (database.category() == ecmdb::category_texture) {
                    static_cast<uint8_t *>(element)[c] = srgb::linear_to_nonlinear(v[c]);
                } else if (database.type() == ecmdb::type_uint8) {
                    static_cast<uint8_t *>(element)[c] = std::round(v[c]);
                } else if (database.type() == ecmdb::type_int16) {
                    static_cast<int16_t *>(element)[c] = std::round(v[c]);
                } else {
                    static_cast<float *>(element)[c] = v[c];
                }
            }
        }
        if (mask.ptr<uint8_t>()[e]) {
            none_valid = false;
        } else {
            all_valid = false;
        }
    }

    ecmdb::metadata quad_metadata;
    metadata::set_quad_meta(database, data, mask, &quad_metadata);

    if (!none_valid) {
        database.save_quad(quad_filename, data.ptr(), mask.ptr<const uint8_t>(), all_valid, &quad_metadata,
                lossy_compression ? 2 : 1, lossy_compression_quality);
    }

    return quad_metadata;
}

/* Helper functions to copy and clear quad areas */

static void copy_box(
        size_t element_size,
        size_t dst_quad_size,
        blob& dst_data,
        int dst_x, int dst_y, int dst_w, int dst_h,
        size_t src_quad_size,
        const blob& src_data,
        int src_x, int src_y,
        bool reverse_x, bool reverse_y)
{
    for (int y = 0; y < dst_h; y++) {
        int sy = (reverse_y ? src_y + dst_h - 1 - y : src_y + y);
        if (reverse_x) {
            for (int x = 0; x < dst_w; x++) {
                int sx = src_x + dst_w - 1 - x;
                const void* s = src_data.ptr((sy * src_quad_size + sx) * element_size);
                void* d = dst_data.ptr(((dst_y + y) * dst_quad_size + (dst_x + x)) * element_size);
                std::memcpy(d, s, element_size);
            }
        } else {
            const void* s = src_data.ptr((sy * src_quad_size + src_x) * element_size);
            void* d = dst_data.ptr(((dst_y + y) * dst_quad_size + dst_x) * element_size);
            std::memcpy(d, s, dst_w * element_size);
        }
    }
}

static void copy_box(
        size_t element_size,
        size_t dst_quad_size,
        blob& dst_data, blob& dst_mask,
        int dst_x, int dst_y, int dst_w, int dst_h,
        size_t src_quad_size,
        const blob& src_data, const blob& src_mask,
        int src_x, int src_y,
        bool reverse_x, bool reverse_y)
{
    copy_box(element_size,    dst_quad_size, dst_data, dst_x, dst_y, dst_w, dst_h, src_quad_size, src_data, src_x, src_y, reverse_x, reverse_y);
    copy_box(sizeof(uint8_t), dst_quad_size, dst_mask, dst_x, dst_y, dst_w, dst_h, src_quad_size, src_mask, src_x, src_y, reverse_x, reverse_y);
}

static void clear_box(
        size_t element_size, size_t dst_quad_size, blob& dst_data,
        int dst_x, int dst_y, int dst_w, int dst_h)
{
    for (int y = 0; y < dst_h; y++) {
        std::memset(dst_data.ptr(((dst_y + y) * dst_quad_size + dst_x) * element_size), 0, dst_w * element_size);
    }
}

static void clear_box(
        size_t element_size, size_t dst_quad_size,
        blob& dst_data, blob& dst_mask,
        int dst_x, int dst_y, int dst_w, int dst_h)
{
    clear_box(element_size,    dst_quad_size, dst_data, dst_x, dst_y, dst_w, dst_h);
    clear_box(sizeof(uint8_t), dst_quad_size, dst_mask, dst_x, dst_y, dst_w, dst_h);
}

static void set_quad_area(
        const std::string& dir, const ecmdb& database,
        quad_cache& qc, quad_load_mutexes& qlm,
        int side, int src_level,
        blob& dst_data, blob& dst_mask,
        int dst_x, int dst_y, int dst_w, int dst_h,
        int src_qx, int src_qy, int src_x, int src_y,
        int fb_src_qx, int fb_src_qy, int fb_src_x, int fb_src_y, bool fb_reverse_x, bool fb_reverse_y)
{
    size_t src_quad_size = database.total_quad_size();
    size_t dst_quad_size = 2 * src_quad_size;

    const quad_data* qd = NULL;
    if (src_qx >= 0 && src_qy >= 0 && src_qx < (1 << src_level) && src_qy < (1 << src_level)
            && database.has_quad(side, src_level, src_qx, src_qy)) {
        quad_location ql(side, src_qx, src_qy);
        qd = qc.locked_get(ql);
        if (!qd) {
            mutex* load_mutex = qlm.get(ql);
            if (load_mutex->trylock()) {
                // this thread wins the race to load this quad
                quad_data* new_qd = NULL;
                size_t new_qd_size;
                load_quad_data(dir, database, side, src_level, src_qx, src_qy, &new_qd, &new_qd_size);
                assert(new_qd);
                qc.locked_put(ql, new_qd, new_qd_size);
                qd = new_qd;
                load_mutex->unlock();
            } else {
                // wait for concurrent thread to finish loading and return result
                msg::dbg("waiting for concurrent thread to finish loading quad...");
                load_mutex->lock();
                load_mutex->unlock();
                qd = qc.locked_get(ql);
                assert(qd);
            }
        }
        if (qd->data.get()->ptr()) {
            copy_box(database.element_size(), dst_quad_size, dst_data, dst_mask, dst_x, dst_y, dst_w, dst_h,
                    src_quad_size, *(qd->data.get()), *(qd->mask.get()), src_x, src_y, false, false);
            return;
        }
    }
    if (fb_src_qx >= 0 && fb_src_qy >= 0 && fb_src_qx < (1 << src_level) && fb_src_qy < (1 << src_level)
            && database.has_quad(side, src_level, fb_src_qx, fb_src_qy)) {
        quad_location ql(side, fb_src_qx, fb_src_qy);
        qd = qc.locked_get(ql);
        if (!qd) {
            mutex* load_mutex = qlm.get(ql);
            if (load_mutex->trylock()) {
                // this thread wins the race to load this quad
                quad_data* new_qd = NULL;
                size_t new_qd_size;
                load_quad_data(dir, database, side, src_level, fb_src_qx, fb_src_qy, &new_qd, &new_qd_size);
                assert(new_qd);
                qc.locked_put(ql, new_qd, new_qd_size);
                qd = new_qd;
                load_mutex->unlock();
            } else {
                // wait for concurrent thread to finish loading and return result
                msg::dbg("waiting for concurrent thread to finish loading quad...");
                load_mutex->lock();
                load_mutex->unlock();
                qd = qc.locked_get(ql);
                assert(qd);
            }
        }
        if (qd->data.get()->ptr()) {
            copy_box(database.element_size(), dst_quad_size, dst_data, dst_mask, dst_x, dst_y, dst_w, dst_h,
                    src_quad_size, *(qd->data.get()), *(qd->mask.get()), fb_src_x, fb_src_y, fb_reverse_x, fb_reverse_y);
            return;
        }
    }
    clear_box(database.element_size(), dst_quad_size, dst_data, dst_mask, dst_x, dst_y, dst_w, dst_h);
}

/* Commit a lower level quad */

static void commit_ll_quad(const std::string& dir, const ecmdb& database,
        bool lossy_compression, int lossy_compression_quality,
        quad_cache& qc, quad_load_mutexes& qlm,
        int qs, int ql, int qx, int qy)
{
    msg::dbg("committing quad %d-%d-%d-%d", qs, ql, qx, qy);
    const std::string quad_basename = ecmdb::quad_filename(qs, ql, qx, qy);
    const std::string base_filename = quad_basename.substr(0, quad_basename.length() - 4);      // remove ".gta"
    const std::string quad_filename = dir + '/' + quad_basename;

    // create quad
    int quad_size = database.quad_size();
    int overlap = database.overlap();
    int total_quad_size = database.total_quad_size();
    blob data(database.data_size());
    blob mask(database.mask_size());

    // gather the required data (from 16 source quads) into one big quad
    blob big_data(4 * database.data_size());
    blob big_mask(4 * database.mask_size());
    // src quad 0
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            0, 0, 2 * overlap, 2 * overlap,
            2 * qx - 1, 2 * qy - 1, total_quad_size - 3 * overlap, total_quad_size - 3 * overlap,
            2 * qx, 2 * qy, overlap, overlap, true, true);
    // src quad 1
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap, 0, quad_size, 2 * overlap,
            2 * qx, 2 * qy - 1, overlap, total_quad_size - 3 * overlap,
            2 * qx, 2 * qy, overlap, overlap, false, true);
    // src quad 2
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap + quad_size, 0, quad_size, 2 * overlap,
            2 * qx + 1, 2 * qy - 1, overlap, total_quad_size - 3 * overlap,
            2 * qx + 1, 2 * qy, overlap, overlap, false, true);
    // src quad 3
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap + 2 * quad_size, 0, 2 * overlap, 2 * overlap,
            2 * qx + 2, 2 * qy - 1, overlap, total_quad_size - 3 * overlap,
            2 * qx + 1, 2 * qy, total_quad_size - 3 * overlap, overlap, true, true);
    // src quad 4
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            0, 2 * overlap, 2 * overlap, quad_size,
            2 * qx - 1, 2 * qy, total_quad_size - 3 * overlap, overlap,
            2 * qx, 2 * qy, overlap, overlap, true, false);
    // src quad 5
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap, 2 * overlap, quad_size, quad_size,
            2 * qx, 2 * qy, overlap, overlap,
            -1, -1, -1, -1, false, false);
    // src quad 6
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap + quad_size, 2 * overlap, quad_size, quad_size,
            2 * qx + 1, 2 * qy, overlap, overlap,
            -1, -1, -1, -1, false, false);
    // src quad 7
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap + 2 * quad_size, 2 * overlap, 2 * overlap, quad_size,
            2 * qx + 2, 2 * qy, overlap, overlap,
            2 * qx + 1, 2 * qy, total_quad_size - 3 * overlap, overlap, true, false);
    // src quad 8
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            0, 2 * overlap + quad_size, 2 * overlap, quad_size,
            2 * qx - 1, 2 * qy + 1, total_quad_size - 3 * overlap, overlap,
            2 * qx, 2 * qy + 1, overlap, overlap, true, false);
    // src quad 9
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap, 2 * overlap + quad_size, quad_size, quad_size,
            2 * qx, 2 * qy + 1, overlap, overlap,
            -1, -1, -1, -1, false, false);
    // src quad 10
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap + quad_size, 2 * overlap + quad_size, quad_size, quad_size,
            2 * qx + 1, 2 * qy + 1, overlap, overlap,
            -1, -1, -1, -1, false, false);
    // src quad 11
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap + 2 * quad_size, 2 * overlap + quad_size, 2 * overlap, quad_size,
            2 * qx + 2, 2 * qy + 1, overlap, overlap,
            2 * qx + 1, 2 * qy + 1, total_quad_size - 3 * overlap, overlap, true, false);
    // src quad 12
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            0, 2 * overlap + 2 * quad_size, 2 * overlap, 2 * overlap,
            2 * qx - 1, 2 * qy + 2, total_quad_size - 3 * overlap, overlap,
            2 * qx, 2 * qy + 1, overlap, total_quad_size - 3 * overlap, true, true);
    // src quad 13
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap, 2 * overlap + 2 * quad_size, quad_size, 2 * overlap,
            2 * qx, 2 * qy + 2, overlap, overlap,
            2 * qx, 2 * qy + 1, overlap, total_quad_size - 3 * overlap, false, true);
    // src quad 14
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap + quad_size, 2 * overlap + 2 * quad_size, quad_size, 2 * overlap,
            2 * qx + 1, 2 * qy + 2, overlap, overlap,
            2 * qx + 1, 2 * qy + 1, overlap, total_quad_size - 3 * overlap, false, true);
    // src quad 15
    set_quad_area(dir, database, qc, qlm, qs, ql + 1, big_data, big_mask,
            2 * overlap + 2 * quad_size, 2 * overlap + 2 * quad_size, 2 * overlap, 2 * overlap,
            2 * qx + 2, 2 * qy + 2, overlap, overlap,
            2 * qx + 1, 2 * qy + 1, total_quad_size - 3 * overlap, total_quad_size - 3 * overlap, true, true);

    // resample the big quad to the new quad
    bool all_valid = true;
    bool none_valid = true;
    for (int y = 0; y < total_quad_size; y++) {
        for (int x = 0; x < total_quad_size; x++) {
            int i = y * total_quad_size + x;
            int i0 = (2 * y + 0) * 2 * total_quad_size + (2 * x + 0);
            int i1 = (2 * y + 0) * 2 * total_quad_size + (2 * x + 1);
            int i2 = (2 * y + 1) * 2 * total_quad_size + (2 * x + 0);
            int i3 = (2 * y + 1) * 2 * total_quad_size + (2 * x + 1);
            uint8_t m0 = big_mask.ptr<uint8_t>()[i0];
            uint8_t m1 = big_mask.ptr<uint8_t>()[i1];
            uint8_t m2 = big_mask.ptr<uint8_t>()[i2];
            uint8_t m3 = big_mask.ptr<uint8_t>()[i3];
            int valid_values = (m0 > 0 ? 1 : 0) + (m1 > 0 ? 1 : 0) + (m2 > 0 ? 1 : 0) + (m3 > 0 ? 1 : 0);
            if (valid_values >= 2) {
                none_valid = false;
                mask.ptr<uint8_t>()[y * total_quad_size + x] = 0xff;
                if (database.category() == ecmdb::category_texture) {
                    // we have type uint8
                    for (int c = 0; c < database.channels(); c++) {
                        int x = 0;
                        if (m0 > 0)
                            x += srgb::nonlinear_to_linear(big_data.ptr<uint8_t>()[database.channels() * i0 + c]);
                        if (m1 > 0)
                            x += srgb::nonlinear_to_linear(big_data.ptr<uint8_t>()[database.channels() * i1 + c]);
                        if (m2 > 0)
                            x += srgb::nonlinear_to_linear(big_data.ptr<uint8_t>()[database.channels() * i2 + c]);
                        if (m3 > 0)
                            x += srgb::nonlinear_to_linear(big_data.ptr<uint8_t>()[database.channels() * i3 + c]);
                        data.ptr<uint8_t>()[database.channels() * i + c] = srgb::linear_to_nonlinear(x / valid_values
                                + (x % valid_values > valid_values / 3 ? 1 : 0));
                    }
                } else if (database.type() == ecmdb::type_uint8) {
                    for (int c = 0; c < database.channels(); c++) {
                        int x = 0;
                        if (m0 > 0)
                            x += big_data.ptr<uint8_t>()[database.channels() * i0 + c];
                        if (m1 > 0)
                            x += big_data.ptr<uint8_t>()[database.channels() * i1 + c];
                        if (m2 > 0)
                            x += big_data.ptr<uint8_t>()[database.channels() * i2 + c];
                        if (m3 > 0)
                            x += big_data.ptr<uint8_t>()[database.channels() * i3 + c];
                        data.ptr<uint8_t>()[database.channels() * i + c] = x / valid_values
                            + (x % valid_values > valid_values / 3 ? 1 : 0);
                    }
                } else if (database.type() == ecmdb::type_int16) {
                    for (int c = 0; c < database.channels(); c++) {
                        float x = 0.0f;
                        if (m0 > 0)
                            x += big_data.ptr<int16_t>()[database.channels() * i0 + c];
                        if (m1 > 0)
                            x += big_data.ptr<int16_t>()[database.channels() * i1 + c];
                        if (m2 > 0)
                            x += big_data.ptr<int16_t>()[database.channels() * i2 + c];
                        if (m3 > 0)
                            x += big_data.ptr<int16_t>()[database.channels() * i3 + c];
                        data.ptr<int16_t>()[database.channels() * i + c] = std::round(x / valid_values);
                    }
                } else {
                    for (int c = 0; c < database.channels(); c++) {
                        float x = 0.0f;
                        if (m0 > 0)
                            x += big_data.ptr<float>()[database.channels() * i0 + c];
                        if (m1 > 0)
                            x += big_data.ptr<float>()[database.channels() * i1 + c];
                        if (m2 > 0)
                            x += big_data.ptr<float>()[database.channels() * i2 + c];
                        if (m3 > 0)
                            x += big_data.ptr<float>()[database.channels() * i3 + c];
                        data.ptr<float>()[database.channels() * i + c] = x / valid_values;
                    }
                }
            } else {
                all_valid = false;
                mask.ptr<uint8_t>()[y * total_quad_size + x] = 0x00;
                std::memset(data.ptr((y * total_quad_size + x) * database.element_size()),
                        0, database.element_size());
            }
        }
    }

    // save
    if (!none_valid) {
        ecmdb::metadata quad_metadata;
        metadata::set_quad_meta(database, data, mask, &quad_metadata);
        fio::mkdir_p(dir, base_filename.substr(0, base_filename.find_last_of('/')));
        database.save_quad(quad_filename, data.ptr(), mask.ptr<const uint8_t>(), all_valid, &quad_metadata,
                lossy_compression ? 2 : 1, lossy_compression_quality);
    }
}
