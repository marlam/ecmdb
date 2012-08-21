/*
 * Copyright (C) 2008, 2009, 2010, 2011, 2012
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

#include <cerrno>
#include <cmath>
#include <limits>

#include <omp.h>

#include <gdal_priv.h>
#include <ogr_geometry.h>
#include <gdalwarper.h>

#include <ecmdb/ecmdb.h>

#include "dbg.h"
#include "msg.h"
#include "str.h"
#include "opt.h"
#include "exc.h"
#include "fio.h"
#include "blob.h"
#include "sys.h"
#include "thread.h"

#include "uuid.h"
#include "quadlist.h"
#include "metadata.h"
#include "compression-info.h"


extern "C" void ecmdb_add_help(void)
{
    msg::req_txt(
            "add DIR [--srs=SRS] [--min=MIN[,...]] [--max=MAX[,...]] datafile...\n"
            "\n"
            "Add the given datafiles to the database at the directory DIR.\n"
            "Currently, only GDAL-readable datafiles with georeference information can be added.\n"
            "The --srs option allows to set/override the spatial reference system (SRS) of the datafile(s), "
            "in the same way that the '-a_srs' option of gdal_translate works.\n"
            "When the --min or --max options are used, then all input values that are less than the minimum "
            "or greater than the maximum are ignored. You can give a single value for all channels, or a list "
            "of values for each channel.");
}

static void add_file(const ecmdb& database, const std::string& database_dir, bool lossy_compression, int lossy_compression_quality,
        const std::string& filename,
        const std::string& addfilename, FILE* addfile,
        const std::string& srs_wkt, const std::vector<float>& minimum, const std::vector<float>& maximum);

extern "C" int ecmdb_add(int argc, char* argv[])
{
    std::vector<opt::option *> options;
    opt::info help("help", '\0', opt::optional);
    options.push_back(&help);
    opt::tuple<float> minimum("min", '\0', opt::optional, std::vector<float>());
    options.push_back(&minimum);
    opt::tuple<float> maximum("max", '\0', opt::optional, std::vector<float>());
    options.push_back(&maximum);
    opt::string srs("srs", '\0', opt::optional);
    options.push_back(&srs);
    std::vector<std::string> arguments;
    if (!opt::parse(argc, argv, options, 2, -1, arguments)) {
        return 1;
    }
    if (help.value()) {
        ecmdb_add_help();
        return 0;
    }

    try {
        ecmdb database;
        database.open(arguments[0]);
        bool lossy_compression;
        int lossy_compression_quality;
        get_compression_info(arguments[0], &lossy_compression, &lossy_compression_quality);

        GDALAllRegister();
        GDALSetCacheMax64(sys::total_ram() / 8);
        std::string srs_wkt;
        if (srs.value().length() > 0) {
            OGRSpatialReference ogr_srs;
            char *tmpstr = NULL;
            if (ogr_srs.SetFromUserInput(srs.value().c_str()) != OGRERR_NONE
                    || ogr_srs.exportToWkt(&tmpstr) != OGRERR_NONE
                    || !tmpstr) {
                throw exc("GDAL cannot understand spatial reference system (SRS)");
            }
            srs_wkt = std::string(tmpstr);
            OGRFree(tmpstr);
        }

        class uuid uuid;
        uuid.generate();
        std::string addfilename = quadlist::get_addfilename(arguments[0], uuid.to_string());
        FILE* addfile = quadlist::open_addfile(addfilename);

        exc e;
        size_t arguments_index = 1;
        #pragma omp parallel for
        for (size_t i = 1; i < arguments.size(); i++) {
            size_t j = atomic::fetch_and_inc(&arguments_index);
            #pragma omp flush (e)
            if (e.empty()) {
                try {
                    add_file(database, arguments[0], lossy_compression, lossy_compression_quality,
                            arguments[j], addfilename, addfile,
                            srs_wkt, minimum.value(), maximum.value());
                }
                catch (std::exception& _e) {
                    e = _e;
                    #pragma omp flush (e)
                }
            }
        }
        if (!e.empty()) {
            throw e;
        }

        fio::close(addfile, addfilename);
    }
    catch (std::exception& e) {
        msg::err_txt("%s", e.what());
        return 1;
    }

    return 0;
}

static bool get_quads(const ecmdb& database,
        GDALDataset *gdal_dataset,
        OGRCoordinateTransformation *ogr_transf,
        int quads_tl[2], int quads_br[2]);

static void add_quad(const ecmdb& database, const std::string& database_dir, bool lossy_compression, int lossy_compression_quality,
        int qs, int ql, int qx, int qy,
        const std::string& addfilename, FILE* addfile,
        GDALDataset* gdal_dataset, const std::string& wkt, const std::string& id, const std::string& tempdir,
        const std::string& srs_wkt, const std::vector<float>& minimum, const std::vector<float>& maximum);

void add_file(const ecmdb& database, const std::string& database_dir, bool lossy_compression, int lossy_compression_quality,
        const std::string& filename,
        const std::string& addfilename, FILE* addfile,
        const std::string& srs_wkt, const std::vector<float>& minimum, const std::vector<float>& maximum)
{
    GDALDataset *gdal_dataset;
    GDALDataType gdal_datatype;
    int64_t gdal_width, gdal_height;
    int gdal_channels;
    GDALRasterBand *gdal_channel;
    double gdal_geo_transf[6];
    std::string proj4[6];
    OGRSpatialReference ogr_qsc_ref[6];
    OGRSpatialReference ogr_dataset_ref;
    OGRCoordinateTransformation *ogr_transf[6];

    // Generate a unique id for this data file and this process
    class uuid uuid;
    uuid.generate();
    std::string id = uuid.to_string();
    msg::inf(filename + ": UUID: " + id);

    // Open file
    if (!(gdal_dataset = reinterpret_cast<GDALDataset *>(GDALOpen(filename.c_str(), GA_ReadOnly)))) {
        throw exc(filename + ": not a GDAL-supported file format");
    }
    msg::inf(filename + ": GDAL driver: "
            + reinterpret_cast<const char *>(gdal_dataset->GetDriver()->GetDescription()) + " / "
            + reinterpret_cast<const char *>(gdal_dataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME)));

    // Get size, channels, and datatype
    gdal_width = gdal_dataset->GetRasterXSize();
    gdal_height = gdal_dataset->GetRasterYSize();
    gdal_channels = gdal_dataset->GetRasterCount();
    if (gdal_width < 1 || gdal_height < 1) {
        throw exc(filename + ": no data");
    }
    if (gdal_channels < 1 || (gdal_channels == 1
                && gdal_dataset->GetRasterBand(1)->GetColorInterpretation() == GCI_AlphaBand)) {
        throw exc(filename + ": no data channels");
    }
    gdal_channel = gdal_dataset->GetRasterBand(1);
    gdal_datatype = gdal_channel->GetRasterDataType();
    msg::inf(filename + str::asprintf(": GDAL sees %sx%s pixels, %d channels, data format %s",
                str::from(gdal_width).c_str(), str::from(gdal_height).c_str(),
                gdal_channels, GDALGetDataTypeName(gdal_datatype)));
    for (int c = 0; c < gdal_channels; c++) {
        gdal_channel = gdal_dataset->GetRasterBand(c + 1);
        if (gdal_channel->GetXSize() != gdal_width
                || gdal_channel->GetYSize() != gdal_height) {
            throw exc(filename + ": channels differ in size");
        }
        if (gdal_channel->GetRasterDataType() != gdal_datatype) {
            throw exc(filename + ": channels differ in data type");
        }
    }
    bool gdal_dataset_has_alpha = (gdal_dataset->GetRasterBand(gdal_channels)->GetColorInterpretation() == GCI_AlphaBand
            || (database.category() == ecmdb::category_texture && gdal_channels == database.channels() + 1));
    if (gdal_channels - (gdal_dataset_has_alpha ? 1 : 0) != database.channels()) {
        throw exc(filename + ": number of bands does not match database number of channels");
    }
    if (minimum.size() > 1 && minimum.size() != static_cast<size_t>(gdal_channels)) {
        throw exc(filename + ": number of bands does not match number of given minimum values");
    }
    if (maximum.size() > 1 && maximum.size() != static_cast<size_t>(gdal_channels)) {
        throw exc(filename + ": number of bands does not match number of given maximum values");
    }
    gdal_dataset->GetGeoTransform(gdal_geo_transf);

    // Get georeference information from srs_wkt or datafile
    if (srs_wkt.length() > 0) {
        char *tmpstr = new char[srs_wkt.length() + 1];
        strcpy(tmpstr, srs_wkt.c_str());
        char *p = tmpstr;
        if (ogr_dataset_ref.importFromWkt(&p) != OGRERR_NONE) {
            // should never happen because GDAL gave us the WKT string
            delete[] tmpstr;
            throw exc(filename + ": cannot use the given SRS?!");
        }
        delete[] tmpstr;
    } else {
        if (!gdal_dataset->GetProjectionRef() || !gdal_dataset->GetProjectionRef()[0]) {
            throw exc(filename + ": no projection reference information");
        }
        char *tmpstr = new char[strlen(gdal_dataset->GetProjectionRef()) + 1];
        strcpy(tmpstr, gdal_dataset->GetProjectionRef());
        char *p = tmpstr;
        if (ogr_dataset_ref.importFromWkt(&p) != OGRERR_NONE) {
            delete[] tmpstr;
            throw exc(filename + ": cannot understand projection reference string");
        }
        delete[] tmpstr;
    }
    // Sanity-check ellipsoid axes.
    // This is deactivated because in some cases the GDAL driver intentionally uses
    // different values. Example: the PDS driver intentionally uses a sphere instead
    // of an ellipsoid for Mars MOC image data that uses SIMPLE_CYLINDRICAL projection,
    // since apparently the NASA folks use spherical formulas.
#if 0
    {
        OGRErr ogr_err = OGRERR_NONE;
        double semi_major_axis = ogr_dataset_ref.GetSemiMajor(&ogr_err);
        if (ogr_err != OGRERR_NONE)
            throw exc(filename + ": cannot get semi major axis of ellipsoid");
        double semi_minor_axis = ogr_dataset_ref.GetSemiMinor(&ogr_err);
        if (ogr_err != OGRERR_NONE)
            throw exc(filename + ": cannot get semi minor axis of ellipsoid");
        if (!std::isfinite(semi_major_axis)
                || semi_major_axis < database.semi_major_axis()
                || semi_major_axis > database.semi_major_axis()
                || !std::isfinite(semi_minor_axis)
                || semi_minor_axis < database.semi_minor_axis()
                || semi_minor_axis > database.semi_minor_axis())
            throw exc(filename + ": ellipsoid differs from database ellipsoid");
    }
#endif

    // Set up transformations to our six cube sides
    std::string base_proj4 =
        std::string("+wktext +proj=qsc +units=m")
        + " +a=" + str::from(database.semi_major_axis())
        + " +b=" + str::from(database.semi_minor_axis());
    proj4[0] = base_proj4 + " +lat_0=0   +lon_0=0";     // front
    proj4[1] = base_proj4 + " +lat_0=0   +lon_0=90";    // right
    proj4[2] = base_proj4 + " +lat_0=0   +lon_0=180";   // back
    proj4[3] = base_proj4 + " +lat_0=0   +lon_0=-90";   // left
    proj4[4] = base_proj4 + " +lat_0=90  +lon_0=0";     // top
    proj4[5] = base_proj4 + " +lat_0=-90 +lon_0=0";     // bottom
    std::string wkt[6];
    for (int i = 0; i < 6; i++) {
        ogr_qsc_ref[i].importFromProj4(proj4[i].c_str());
        char *tmpstr;
        ogr_qsc_ref[i].exportToWkt(&tmpstr);
        wkt[i] = std::string(tmpstr);
        CPLFree(tmpstr);
        if (!(ogr_transf[i] = OGRCreateCoordinateTransformation(&ogr_dataset_ref, &ogr_qsc_ref[i]))) {
            throw exc(filename + ": cannot create transformation from dataset to ECM side(s)");
        }
    }

    // Loop over all sides of the Ellipsoid Cube Map and add the quads
    const char *side_name[6] = { "front", "right", "back", "left", "top", "bottom" };
    std::string tempdir = fio::mktempdir("ecmdb-add-");
    for (int side = 0; side < 6; side++) {
        int quads_tl[2], quads_br[2];
        if (get_quads(database, gdal_dataset, ogr_transf[side], quads_tl, quads_br)) {
            msg::inf(filename + ": side " + side_name[side]
                    + str::asprintf(": quads from %d,%d to %d,%d",
                        quads_tl[0], quads_tl[1], quads_br[0], quads_br[1]));
            for (int qy = quads_tl[1]; qy <= quads_br[1]; qy++) {
                for (int qx = quads_tl[0]; qx <= quads_br[0]; qx++) {
                    add_quad(database, database_dir, lossy_compression, lossy_compression_quality,
                            side, database.levels() - 1, qx, qy,
                            addfilename, addfile, gdal_dataset, wkt[side], id, tempdir, srs_wkt, minimum, maximum);
                }
            }
        } else {
            msg::inf(filename + ": side " + side_name[side] + ": no quads");
        }
    }
    fio::rm_r(tempdir);

    // Cleanup
    for (int i = 0; i < 6; i++) {
        OCTDestroyCoordinateTransformation(ogr_transf[i]);
    }
    GDALClose(gdal_dataset);
}

bool get_quads(const ecmdb& database,
        GDALDataset *gdal_dataset,
        OGRCoordinateTransformation *ogr_transf,
        int quads_tl[2], int quads_br[2])
{
#if 1
    int quads_in_level = (1 << (database.levels() - 1));
    int w = gdal_dataset->GetRasterXSize();
    int h = gdal_dataset->GetRasterYSize();
    double geo_transform[6];
    gdal_dataset->GetGeoTransform(geo_transform);

    // Test a subset of points in the data set.
    static const int testpoints_per_quadsize = 4;
    int point_dist = std::max(1, database.quad_size() / testpoints_per_quadsize);
    int testpoints = checked_mul(w / point_dist + 2, h / point_dist + 2);
    blob bx(testpoints, sizeof(double));
    blob by(testpoints, sizeof(double));
    blob br(testpoints, sizeof(int));
    int n = 0;
    // Interior
    for (int y = point_dist; y < h; y += point_dist) {
        for (int x = point_dist; x < w; x += point_dist) {
            assert(n < testpoints);
            bx.ptr<double>()[n] = geo_transform[0] + x * geo_transform[1] + y * geo_transform[2];
            by.ptr<double>()[n] = geo_transform[3] + x * geo_transform[4] + y * geo_transform[5];
            n++;
        }
    }
    // Corners
    assert(n + 4 < testpoints);
    bx.ptr<double>()[n] = geo_transform[0] +   0   * geo_transform[1] +   0   * geo_transform[2];
    by.ptr<double>()[n] = geo_transform[3] +   0   * geo_transform[4] +   0   * geo_transform[5];
    n++;
    bx.ptr<double>()[n] = geo_transform[0] + (w-1) * geo_transform[1] +   0   * geo_transform[2];
    by.ptr<double>()[n] = geo_transform[3] + (w-1) * geo_transform[4] +   0   * geo_transform[5];
    n++;
    bx.ptr<double>()[n] = geo_transform[0] +   0   * geo_transform[1] + (h-1) * geo_transform[2];
    by.ptr<double>()[n] = geo_transform[3] +   0   * geo_transform[4] + (h-1) * geo_transform[5];
    n++;
    bx.ptr<double>()[n] = geo_transform[0] + (w-1) * geo_transform[1] + (h-1) * geo_transform[2];
    by.ptr<double>()[n] = geo_transform[3] + (w-1) * geo_transform[4] + (h-1) * geo_transform[5];
    n++;

    // Transform the points to side coordinates and determine the minimum
    // and maximum quad.
    ogr_transf->TransformEx(n, bx.ptr<double>(), by.ptr<double>(), NULL, br.ptr<int>());
    quads_br[0] = -1;
    quads_br[1] = -1;
    quads_tl[0] = (1 << (database.levels() - 1));
    quads_tl[1] = (1 << (database.levels() - 1));
    for (int i = 0; i < n; i++) {
        if (br.ptr<int>()[i]) {
            double bxx = +bx.ptr<double>()[i] / database.semi_major_axis() / 2.0 + 0.5;
            double byy = -by.ptr<double>()[i] / database.semi_major_axis() / 2.0 + 0.5;
            if (bxx >= 0.0 && bxx <= 1.0 && byy >= 0.0 && byy <= 1.0) {
                int qx = bxx * quads_in_level;
                int qy = byy * quads_in_level;
                // Also include all 8 surrounding quads just to be sure.
                // If they contain no data, add_quad() will not write them.
                if (qx - 1 < quads_tl[0])
                    quads_tl[0] = qx - 1;
                if (qx + 1 > quads_br[0])
                    quads_br[0] = qx + 1;
                if (qy - 1 < quads_tl[1])
                    quads_tl[1] = qy - 1;
                if (qy + 1 > quads_br[1])
                    quads_br[1] = qy + 1;
            }
        }
    }
    if (quads_tl[0] < 0)
        quads_tl[0] = 0;
    if (quads_tl[1] < 0)
        quads_tl[1] = 0;
    if (quads_br[0] >= quads_in_level)
        quads_br[0] = quads_in_level - 1;
    if (quads_br[1] >= quads_in_level)
        quads_br[1] = quads_in_level - 1;
    return (quads_br[0] >= quads_tl[0]);
#endif
#if 0
    quads_tl[0] = 0;
    quads_tl[1] = 0;
    quads_br[0] = (1 << (database.levels() - 1)) - 1;
    quads_br[1] = (1 << (database.levels() - 1)) - 1;
    return true;
#endif
}

void add_quad(const ecmdb& database, const std::string& database_dir, bool lossy_compression, int lossy_compression_quality,
        int qs, int ql, int qx, int qy,
        const std::string& addfilename, FILE* addfile,
        GDALDataset* gdal_dataset, const std::string& wkt, const std::string& id, const std::string& tempdir,
        const std::string& srs_wkt, const std::vector<float>& minimum, const std::vector<float>& maximum)
{
    int quads_in_level = (1 << ql);
    int total_quad_size = database.quad_size() + 2 * database.overlap();

    // Set up alpha and no-data information. Alpha is assumed to always be in the last band.
    bool src_has_alpha = (gdal_dataset->GetRasterBand(gdal_dataset->GetRasterCount())->GetColorInterpretation() == GCI_AlphaBand
            || (database.category() == ecmdb::category_texture && gdal_dataset->GetRasterCount() == database.channels() + 1));
    int dst_rastercount = gdal_dataset->GetRasterCount() + (src_has_alpha ? 0 : 1);
    assert(dst_rastercount == database.channels() + 1);
    int *src_bands = static_cast<int *>(CPLMalloc((dst_rastercount - 1) * sizeof(int)));
    int *dst_bands = static_cast<int *>(CPLMalloc((dst_rastercount - 1) * sizeof(int)));
    if (!src_bands || !dst_bands) {
        throw exc(ENOMEM);
    }
    for (int i = 0; i < dst_rastercount - 1; i++) {
        src_bands[i] = i + 1;
        dst_bands[i] = i + 1;
    }
    double *nodata_real = static_cast<double *>(CPLMalloc(gdal_dataset->GetRasterCount() * sizeof(double)));
    double *nodata_imag = static_cast<double *>(CPLMalloc(gdal_dataset->GetRasterCount() * sizeof(double)));
    if (!nodata_real || !nodata_imag) {
        throw exc(ENOMEM);
    }
    for (int i = 0; i < gdal_dataset->GetRasterCount(); i++) {
        nodata_real[i] = gdal_dataset->GetRasterBand(i + 1)->GetNoDataValue(NULL);
        nodata_imag[i] = 0.0;
    }

    // Create output dataset
    std::string dstfilename = tempdir + "/x";
    try { fio::remove(dstfilename); } catch (...) {}
    GDALDataset *gdal_dst_dataset = NULL;
    GDALDriverH gdal_driver = GDALGetDriverByName("MEM");
    if (gdal_driver) {
        gdal_dst_dataset = reinterpret_cast<GDALDataset *>(GDALCreate(gdal_driver, dstfilename.c_str(),
                    total_quad_size, total_quad_size, dst_rastercount,
                    gdal_dataset->GetRasterBand(1)->GetRasterDataType(), NULL));
    }
    if (!gdal_dst_dataset) {
        throw exc("GDAL cannot create data set for quad");
    }
    gdal_dst_dataset->SetProjection(wkt.c_str());

    // Set output area
    double qxx = (static_cast<double>(qx) / quads_in_level - 0.5) * 2.0 * database.semi_major_axis();
    double qyy = (static_cast<double>(qy) / quads_in_level - 0.5) * 2.0 * database.semi_major_axis();
    double ps = 2.0 * database.semi_major_axis() / quads_in_level / database.quad_size();
    double geo_transform[6] = { qxx, ps, 0.0, -qyy, 0.0, -ps };
    gdal_dst_dataset->SetGeoTransform(geo_transform);

    // Set warping options
    GDALWarpOptions *warp_options = GDALCreateWarpOptions();
    warp_options->papszWarpOptions = CSLSetNameValue(warp_options->papszWarpOptions,
            "SOURCE_EXTRA", str::from(std::min(total_quad_size / 2, 64)).c_str());
    warp_options->papszWarpOptions = CSLSetNameValue(warp_options->papszWarpOptions,
            "SAMPLE_GRID", "YES");
    warp_options->dfWarpMemoryLimit = sys::total_ram() / 8 / omp_get_num_procs();
    warp_options->eResampleAlg = GRA_Bilinear;
    warp_options->hSrcDS = reinterpret_cast<GDALDatasetH>(gdal_dataset);
    warp_options->hDstDS = reinterpret_cast<GDALDatasetH>(gdal_dst_dataset);
    warp_options->nBandCount = dst_rastercount - 1;
    warp_options->panSrcBands = src_bands;
    warp_options->panDstBands = dst_bands;
    warp_options->nSrcAlphaBand = src_has_alpha ? gdal_dataset->GetRasterCount() : 0;
    warp_options->nDstAlphaBand = dst_rastercount;
    warp_options->padfSrcNoDataReal = nodata_real;
    warp_options->padfSrcNoDataImag = nodata_imag;
    warp_options->padfDstNoDataReal = NULL;
    warp_options->padfDstNoDataImag = NULL;
    warp_options->pTransformerArg = GDALCreateGenImgProjTransformer(
            warp_options->hSrcDS, srs_wkt.length() > 0 ? srs_wkt.c_str() : NULL,
            warp_options->hDstDS, NULL,
            TRUE, 0.0, 0);
    warp_options->pfnTransformer = GDALGenImgProjTransform;

    // Warp
    GDALWarpOperation warp_op;
    warp_op.Initialize(warp_options);
    warp_op.ChunkAndWarpImage(0, 0, total_quad_size, total_quad_size);
    GDALDestroyGenImgProjTransformer(warp_options->pTransformerArg);
    //GDALDestroyWarpOptions(warp_options);
    CSLDestroy(warp_options->papszWarpOptions);
    CPLFree(warp_options);

    // Get the mask
    blob mask(database.mask_size());
    GDALRasterBand *gdal_dst_mask_band = gdal_dst_dataset->GetRasterBand(dst_rastercount);
    CPLErr e = gdal_dst_mask_band->RasterIO(
            GF_Read, 0, 0, total_quad_size, total_quad_size,
            mask.ptr(), total_quad_size, total_quad_size, GDT_Byte, 0, 0);
    if (e != CE_None) {
        throw exc(dstfilename + ": GDAL read error");
    }
    // Check if we have valid data
    bool all_valid = true;
    bool none_valid = true;
    for (int i = 0; i < database.total_quad_size() * database.total_quad_size(); i++) {
        if (mask.ptr<uint8_t>()[i]) {
            none_valid = false;
        } else {
            all_valid = false;
        }
        if (!all_valid && !none_valid)
            break;
    }

    // Save result
    blob data;
    if (!none_valid) {
        data.resize(database.data_size());
        CPLErr e = gdal_dst_dataset->RasterIO(
                GF_Read, 0, 0, total_quad_size, total_quad_size,
                data.ptr(), total_quad_size, total_quad_size,
                (database.type() == ecmdb::type_uint8 ? GDT_Byte
                 : database.type() == ecmdb::type_int16 ? GDT_Int16 : GDT_Float32),
                database.channels(), NULL,
                /* pixel space */ database.element_size(),
                /* line space  */ database.element_size() * database.total_quad_size(),
                /* band space  */ database.type_size());
        if (e != CE_None) {
            throw exc(dstfilename + ": GDAL read error");
        }
        // Apply the given min/max boundaries, and check for valid values
        if (database.type() == ecmdb::type_float32 || minimum.size() > 0 || maximum.size() > 0) {
            all_valid = true;
            none_valid = true;
            for (int e = 0; e < database.total_quad_size() * database.total_quad_size(); e++) {
                if (mask.ptr<uint8_t>()[e]) {
                    bool e_valid = true;
                    for (int c = 0; c < database.channels() && e_valid; c++) {
                        float mi = (minimum.size() == 0 ? std::numeric_limits<float>::quiet_NaN() : minimum.size() > 1 ? minimum[c] : minimum[0]);
                        float ma = (maximum.size() == 0 ? std::numeric_limits<float>::quiet_NaN() : maximum.size() > 1 ? maximum[c] : maximum[0]);
                        float x;
                        if (database.type() == ecmdb::type_uint8) {
                            x = data.ptr<uint8_t>()[database.channels() * e + c];
                        } else if (database.type() == ecmdb::type_int16) {
                            x = data.ptr<int16_t>()[database.channels() * e + c];
                        } else {
                            assert(database.type() == ecmdb::type_float32);
                            x = data.ptr<float>()[database.channels() * e + c];
                        }
                        x = database.data_offset() + database.data_factor() * x;
                        if (!std::isfinite(x))
                            e_valid = false;
                        if (minimum.size() > 0 && x < mi)
                            e_valid = false;
                        if (maximum.size() > 0 && x > ma)
                            e_valid = false;
                    }
                    if (!e_valid) {
                        mask.ptr<uint8_t>()[e] = 0x00;
                    }
                }
                if (mask.ptr<uint8_t>()[e]) {
                    none_valid = false;
                } else {
                    all_valid = false;
                }
            }
        }
    }
    if (!none_valid) {
        // Set all invalid values to zero in the data, for better compression.
        // Also make sure that mask values are either 0x00 or 0xff.
        for (int e = 0; e < database.total_quad_size() * database.total_quad_size(); e++) {
            if (mask.ptr<uint8_t>()[e]) {
                mask.ptr<uint8_t>()[e] = 0xff;
            } else {
                std::memset(data.ptr(e * database.element_size()), 0, database.element_size());
            }
        }
        // Save metadata
        ecmdb::quad_metadata quad_meta;
        metadata::set_quad_meta(database, data, mask, &quad_meta);
        // Save the quad
        std::string basename = ecmdb::quad_filename(qs, ql, qx, qy);
        basename = basename.substr(0, basename.length() - 4);   // remove ".gta"
        fio::mkdir_p(database_dir, basename.substr(0, basename.find_last_of('/')));
        std::string quadfilename = database_dir + '/' + basename + '.' + id + ".gta";
        database.save_quad(quadfilename, data.ptr(), mask.ptr<const uint8_t>(), all_valid, &quad_meta,
                lossy_compression ? 2 : 1, lossy_compression_quality);
        // Save an entry in the addfile
        quadlist::add_entry(addfilename, addfile, quadlist::entry(id, qs, qx, qy));
    }

    // Cleanup
    GDALClose(gdal_dst_dataset);
}
