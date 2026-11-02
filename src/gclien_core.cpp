// Copyright (c) 2026 Jun'ya Takakura, National Institute for Environmental Studies (NIES)
// SPDX-License-Identifier: BSD-3-Clause
//
// gclien_core.cpp
// Implementation of the processing shared by WFDE5 / ISIMIP3b.

#include <iostream>
#include <filesystem>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <math.h>
#include <vector>
#include <map>
#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <omp.h>

#include <netcdf.h>
#include <xgboost/c_api.h>

#include "wbgtn.h"
#include "gclien_core.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::array;
using namespace std::chrono;
using namespace std::chrono_literals;


// Logging
int g_log_utc_offset = 9;  // Overwritten after the config is loaded

static auto jst_now() {
    return floor<seconds>(system_clock::now() + std::chrono::hours{g_log_utc_offset});
}

void information(string msg) {
    cout << "[" << std::format("{:%Y-%m-%d %H:%M:%S}", jst_now()) << "] [Info] " << msg << endl;
}

// For printing config contents. Many lines are printed together at startup, so a timestamp on each line is redundant.
// Unlike information(), output has no timestamp (content only).
void config_line(string msg) {
    cout << msg << endl;
}

// Error output for NetCDF operations
[[noreturn]] void nc_handle_error(int status, string filepath) {
    cerr << "[" << std::format("{:%Y-%m-%d %H:%M:%S}", jst_now()) << "] [Error] " << nc_strerror(status) << " " << filepath << endl;
    exit(EXIT_FAILURE);
}

// Error output for non-NetCDF operations
[[noreturn]] void handle_error(string msg, string filepath) {
    cerr << "[" << std::format("{:%Y-%m-%d %H:%M:%S}", jst_now()) << "] [Error] " << msg << " " << filepath << endl;
    exit(EXIT_FAILURE);
}

// Non-throwing numeric parsing (avoids the exceptions of std::stoi/std::stod).
// Like std::stoi/std::stod, the leading number is interpreted and trailing extra characters are ignored.
bool parse_int(const std::string& s, int& out) {
    const char* b = s.data();
    const char* e = b + s.size();
    auto res = std::from_chars(b, e, out);
    return res.ec == std::errc{};
}

bool parse_double(const std::string& s, double& out) {
    const char* b = s.c_str();
    char* end = nullptr;
    out = std::strtod(b, &end);
    return end != b;
}

// Get an environment variable
string get_environment_variable(string const& variable_name) {
    size_t size;
    getenv_s(&size, nullptr, 0, variable_name.c_str());
    if (size > 0) {
        char* tmpvar = new char[size];
        errno_t result = getenv_s(&size, tmpvar, size, variable_name.c_str());
        string var = (result == 0 ? string(tmpvar) : "");
        delete[] tmpvar;
        return var;
    }
    return "";
}

// Check whether a string is numeric (an empty string is not numeric)
bool is_number(const string& str)
{
    if (str.empty()) return false;
    for (char const& c : str) {
        if (std::isdigit(c) == 0) return false;
    }
    return true;
}

// Load a config file (INI parser)
// Checks that the file exists and returns false if it cannot be opened (no exceptions).
bool load_ini(const std::string& filepath, IniData& data) {
    data.clear();
    std::error_code ec;
    if (!std::filesystem::exists(filepath, ec) || ec) {
        return false;
    }
    std::ifstream f(filepath);
    if (!f.is_open()) {
        return false;
    }
    std::string line, section;
    auto trim = [](std::string& s) {
        s.erase(0, s.find_first_not_of(" \t\r\n"));
        auto last = s.find_last_not_of(" \t\r\n");
        if (last != std::string::npos) s.erase(last + 1);
        else s.clear();
    };
    while (std::getline(f, line)) {
        trim(line);
        if (line.empty() || line[0] == ';' || line[0] == '#') continue;
        if (line[0] == '[') {
            auto end = line.find(']');
            if (end != std::string::npos) section = line.substr(1, end - 1);
        } else {
            auto eq = line.find('=');
            if (eq != std::string::npos && !section.empty()) {
                std::string key = line.substr(0, eq);
                std::string val = line.substr(eq + 1);
                trim(key);
                trim(val);
                data[section][key] = val;
            }
        }
    }
    return true;
}

std::string get_config_str(const IniData& ini, const char* section, const char* key) {
    auto s = ini.find(section);
    if (s == ini.end()) return "";
    auto v = s->second.find(key);
    if (v == s->second.end()) return "";
    return v->second;
}

// Read a required string key. Exits immediately with a clear message if it is missing or empty.
// Silently continuing with an empty string would cause obscure failures later in nc_open/XGBoost.
std::string get_config_str_required(const IniData& ini, const char* section, const char* key) {
    std::string val = get_config_str(ini, section, key);
    if (val.empty()) {
        handle_error(std::string("Config [") + section + "] " + key + " not found or empty.", "");
    }
    return val;
}

int cfg_stoi(const IniData& ini, const char* section, const char* key) {
    std::string val = get_config_str(ini, section, key);
    if (val.empty())
        handle_error(std::string("Config [") + section + "] " + key + " not found or empty.", "");
    int out = 0;
    if (!parse_int(val, out))
        handle_error(std::string("Config [") + section + "] " + key + " has invalid integer value: '" + val + "'.", "");
    return out;
}

double cfg_stod(const IniData& ini, const char* section, const char* key) {
    std::string val = get_config_str(ini, section, key);
    if (val.empty())
        handle_error(std::string("Config [") + section + "] " + key + " not found or empty.", "");
    double out = 0.0;
    if (!parse_double(val, out))
        handle_error(std::string("Config [") + section + "] " + key + " has invalid double value: '" + val + "'.", "");
    return out;
}

double cfg_stod_default(const IniData& ini, const char* section, const char* key, double def) {
    std::string val = get_config_str(ini, section, key);
    if (val.empty()) return def;
    double out = 0.0;
    if (!parse_double(val, out))
        handle_error(std::string("Config [") + section + "] " + key + " has invalid double value: '" + val + "'.", "");
    return out;
}

int cfg_stoi_default(const IniData& ini, const char* section, const char* key, int def) {
    std::string val = get_config_str(ini, section, key);
    if (val.empty()) return def;
    int out = 0;
    if (!parse_int(val, out))
        handle_error(std::string("Config [") + section + "] " + key + " has invalid integer value: '" + val + "'.", "");
    return out;
}


// Bilinear interpolation
// Weighted mean with bilinear weights over the valid (non-missing) grid points. Returns fillVal if all are missing.
static double getBilinearInterpolation(double x, double y, double x1, double x2, double y1, double y2, double q11, double q12, double q21, double q22, double fillVal) {
    double w11 = (x2 - x) * (y2 - y);
    double w12 = (x2 - x) * (y - y1);
    double w21 = (x - x1) * (y2 - y);
    double w22 = (x - x1) * (y - y1);

    double a = 0, b = 0;
    if (q11 != fillVal) { a += w11 * q11; b += w11; }
    if (q12 != fillVal) { a += w12 * q12; b += w12; }
    if (q21 != fillVal) { a += w21 * q21; b += w21; }
    if (q22 != fillVal) { a += w22 * q22; b += w22; }

    return (b != 0.0) ? a / b : fillVal;
}

/**
* @brief Get the list of files under a folder
* @param[in]    folderPath  folder path
* @param[out]   file_names  list of file names
* return        true: success, false: failure
*/
bool getFileNames(std::string folderPath, std::vector<std::string>& file_names)
{
    using namespace std::filesystem;
    std::error_code err;
    directory_iterator iter(folderPath, err), end;  // ec overload (does not throw)
    if (err) {
        return false;
    }

    for (; iter != end && !err; iter.increment(err)) {
        const directory_entry entry = *iter;

        file_names.push_back(entry.path().string());
    }

    /* Error handling */
    if (err) {
        return false;
    }
    return true;
}


int save_as_ncdf(string netcdf_filepath_out, const std::vector<float>& vars_to_save, const string varname, const std::vector<double>& vec_time, const std::vector<double>& vec_lat, const std::vector<double>& vec_lon, int compression_level, int quantize_nsd = 0)
{
    // Create and define the result file (WBGT NetCDF file)
// -------------------------------------------------------------------------------------------
    // Create a new NetCDF file -> get NCID
    int ncid_out, dimid_time, dimid_lat, dimid_lon, varid, varid_time, varid_lat, varid_lon;
    int dimids[3], dimids_time[1], dimids_lat[1], dimids_lon[1];
    int status;

    status = nc_create(netcdf_filepath_out.c_str(), NC_NETCDF4, &ncid_out);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    // Define dimensions
    status = nc_def_dim(ncid_out, "time", vec_time.size(), &dimid_time);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);
    status = nc_def_dim(ncid_out, "latitude", vec_lat.size(), &dimid_lat);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);
    status = nc_def_dim(ncid_out, "longitude", vec_lon.size(), &dimid_lon);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    dimids[0] = dimid_time;
    dimids[1] = dimid_lat;
    dimids[2] = dimid_lon;

    dimids_time[0] = dimid_time;
    dimids_lat[0] = dimid_lat;
    dimids_lon[0] = dimid_lon;

    // Define the data variable to write
    status = nc_def_var(ncid_out, varname.c_str(), NC_FLOAT, 3, dimids, &varid);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    status = nc_def_var(ncid_out, "time", NC_DOUBLE, 1, dimids_time, &varid_time);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    status = nc_def_var(ncid_out, "latitude", NC_DOUBLE, 1, dimids_lat, &varid_lat);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    status = nc_def_var(ncid_out, "longitude", NC_DOUBLE, 1, dimids_lon, &varid_lon);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);


    // Set the fill value
    float fill_float = NC_FILL_FLOAT;
    status = nc_put_att_float(ncid_out, varid, "_FillValue", NC_FLOAT, 1, &fill_float);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    // Compression settings
    int shuffle       = (compression_level != 0) ? 1 : 0;
    int deflate       = (compression_level != 0) ? 1 : 0;
    int deflate_level = compression_level;

    status = nc_def_var_deflate(ncid_out, varid, shuffle, deflate, deflate_level);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    // Quantization (lossy). Rounding off low-order bits while keeping quantize_nsd significant digits before
    // deflate makes the files much smaller (write time is unchanged). GranularBitRound balances error and compression ratio well.
    // Not applied when quantize_nsd <= 0 (lossless, legacy behavior).
    if (quantize_nsd > 0) {
        status = nc_def_var_quantize(ncid_out, varid, NC_QUANTIZE_GRANULARBR, quantize_nsd);
        if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);
    }

    // End define mode
    status = nc_enddef(ncid_out);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    // Write dimension data
    status = nc_put_var_double(ncid_out, varid_time, &vec_time[0]);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    status = nc_put_var_double(ncid_out, varid_lat, &vec_lat[0]);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    status = nc_put_var_double(ncid_out, varid_lon, &vec_lon[0]);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    // Write the WBGT data (main data; checked because a missed write failure would silently leave an incomplete file)
    status = nc_put_var_float(ncid_out, varid, &vars_to_save[0]);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    // Close the file
    status = nc_close(ncid_out);
    if (status != NC_NOERR) nc_handle_error(status, netcdf_filepath_out);

    return 0;
}

// Read an NC_FLOAT variable from a file. When open_required=false, returns false if the file does not exist.
bool load_nc_2d_float(const string& filepath, const string& varname,
                      std::vector<float>& data, float* fill_value,
                      bool open_required = true,
                      size_t start_lat = 0, size_t start_lon = 0,
                      size_t count_lat = 0, size_t count_lon = 0)
{
    int status, ncid;
    status = nc_open(filepath.c_str(), NC_NOWRITE, &ncid);
    if (status != NC_NOERR) {
        if (open_required) nc_handle_error(status, filepath);
        return false;
    }
    information(filepath);

    int nvars;
    status = nc_inq(ncid, nullptr, &nvars, nullptr, nullptr);
    if (status != NC_NOERR) nc_handle_error(status, filepath);

    char var_name[NC_MAX_NAME];
    nc_type var_xtype;
    int var_ndims, var_dimids[NC_MAX_VAR_DIMS], var_natts;
    int var_id = -1;
    for (int i = 0; i < nvars; i++) {
        status = nc_inq_var(ncid, i, var_name, &var_xtype, &var_ndims, var_dimids, &var_natts);
        if (status != NC_NOERR) nc_handle_error(status, filepath);
        if (strcmp(var_name, varname.c_str()) == 0 && var_xtype == NC_FLOAT)
            var_id = i;
    }
    if (var_id < 0) handle_error("NetCDFfile invalid variable dimid.", filepath);

    if (nc_get_att_float(ncid, var_id, "_FillValue", fill_value)) {
        std::cerr << "Error getting _FillValue attribute" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (count_lat > 0 && count_lon > 0) {
        // Subset read: read all data once and then extract in memory
        // Not using nc_get_vara_float guarantees numerical agreement with nc_get_var_float
        int ndims, dimids[NC_MAX_VAR_DIMS];
        status = nc_inq_varndims(ncid, var_id, &ndims);
        if (status != NC_NOERR || ndims < 2) nc_handle_error(status, filepath);
        status = nc_inq_var(ncid, var_id, nullptr, nullptr, nullptr, dimids, nullptr);
        if (status != NC_NOERR) nc_handle_error(status, filepath);
        size_t full_lat, full_lon;
        nc_inq_dimlen(ncid, dimids[ndims - 2], &full_lat);
        nc_inq_dimlen(ncid, dimids[ndims - 1], &full_lon);
        std::vector<float> full_data(full_lat * full_lon);
        status = nc_get_var_float(ncid, var_id, full_data.data());
        if (status != NC_NOERR) nc_handle_error(status, filepath);
        for (size_t r = 0; r < count_lat; r++)
            for (size_t c = 0; c < count_lon; c++)
                data[r * count_lon + c] = full_data[(r + start_lat) * full_lon + (c + start_lon)];
    } else {
        status = nc_get_var_float(ncid, var_id, data.data());
        if (status != NC_NOERR) nc_handle_error(status, filepath);
    }

    status = nc_close(ncid);
    if (status != NC_NOERR) nc_handle_error(status, filepath);

    return true;
}


int load_txt_map(string filepath_map_input,  std::map<string, string> &map_from_file, std::vector<string> &keyvec)
{
    std::ifstream fin;
    string str_key, str_val;

    map_from_file.clear();
    keyvec.clear();

    // Checks that the file exists and returns non-zero if it cannot be opened (no exceptions).
    std::error_code ec;
    if (!std::filesystem::exists(filepath_map_input, ec) || ec) {
        return -1;
    }
    fin.open(filepath_map_input.c_str());
    if (!fin.is_open()) {
        return -1;
    }
    while (fin >> str_key >> str_val) {
        keyvec.push_back(str_key);
        map_from_file[str_key] = str_val;
    }
    fin.close();

    return 0;
}

// For debugging
void show_map(const std::map<string, string> &map_to_show, const std::vector<string> &keyvec_to_show)
{
    string str_key, str_val;

    for (int i = 0; i < keyvec_to_show.size(); i++) {
        str_key = keyvec_to_show[i];
        auto it = map_to_show.find(str_key);
        str_val = (it != map_to_show.end()) ? it->second : "";
        config_line(str_key + ":" + str_val);
    }
}

// For debugging: print features
void show_features(const std::vector<float> features_all, const std::vector<string> feature_name_vec, const std::vector<double>& vec_lon_fine, const std::vector<double>& vec_lat_fine, const double lon_show, const double lat_show)
{
    double elon, elon_min = 999, elat, elat_min = 999;
    int idx_lat_show = 0, idx_lon_show = 0;

    for (int i = 0; i < vec_lat_fine.size(); i++) {
        elat = abs(vec_lat_fine[i] - lat_show);
        if (i == 0 || elat < elat_min) {
            elat_min = elat;
            idx_lat_show = i;
        }
    }

    for (int i = 0; i < vec_lon_fine.size(); i++) {
        elon = abs(vec_lon_fine[i] - lon_show);
        if (i == 0 || elon < elon_min) {
            elon_min = elon;
            idx_lon_show = i;
        }
    }

    cout << "idx_lon:" << idx_lon_show << endl;
    cout << "lon;" << lon_show << endl;
    cout << "idx_lat:" << idx_lat_show << endl;
    cout << "lat;" << lat_show << endl;

    size_t row_cnt = (size_t)idx_lat_show * (size_t)vec_lon_fine.size() + (size_t)idx_lon_show;
    size_t idx_feature_vec_first = (size_t)feature_name_vec.size() * row_cnt;

    for (int i = 0; i < feature_name_vec.size(); i++) {
        cout << feature_name_vec[i] << ": " << features_all[idx_feature_vec_first + (size_t)i] << endl;
    }
    cout << endl;
}


// Rearrange features into model input order
void arrange_features(const std::vector<std::vector<float>> &features_normal, const std::vector<std::vector<float>> features_cza_original, const std::vector<double> &vec_lon_fine, const std::vector<double> &vec_lat_fine, int hour, int num_threads, std::vector<float> &features_all)
{
    int hour_lst;
    int total_feature_num = features_normal.size() + 24 + 24;
    int lon_len_fine = vec_lon_fine.size();

    information(std::format("Arranging features for hour {}", hour));

    for (int lon_cnt_fine = 0; lon_cnt_fine < lon_len_fine; lon_cnt_fine++) {

        double lon_fine = vec_lon_fine[lon_cnt_fine];

        hour_lst = hour + (int)round(lon_fine / 15.0);
        hour_lst = (hour_lst + 24) % 24;

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
        for (int lat_cnt_fine = 0; lat_cnt_fine < vec_lat_fine.size(); lat_cnt_fine++) {

            double lat_fine = vec_lat_fine[lat_cnt_fine];

            int col_cnt;
            size_t row_cnt = (size_t)lon_len_fine * (size_t)lat_cnt_fine + (size_t)lon_cnt_fine;
            size_t idx_feature_vec_first = (size_t)total_feature_num * (size_t)row_cnt;
            size_t idx_feature_vec;

            // normal features
            for (col_cnt = 0; col_cnt < features_normal.size(); col_cnt++) {
                idx_feature_vec = idx_feature_vec_first + (size_t)col_cnt;

                features_all[idx_feature_vec] = features_normal[col_cnt][row_cnt];
            }

            // local solar hour dummy
            for (int h = 0; h < 24; h++) {
                col_cnt = features_normal.size() + h;
                idx_feature_vec = idx_feature_vec_first + (size_t)col_cnt;

                if (h == hour_lst) {
                    features_all[idx_feature_vec] = 1;
                }
                else {
                    features_all[idx_feature_vec] = 0;
                }

            }

            // czas for past 24 hours
            for (int h = 0; h < 24; h++) {
                col_cnt = features_normal.size() + 24 + h;
                idx_feature_vec = idx_feature_vec_first + (size_t)col_cnt;

                int h_past = ((24 + hour) - h) % 24;

                features_all[idx_feature_vec] = features_cza_original[h_past][row_cnt];
            }
        }
    }
}


// Load the common config sections
CoreConfig load_core_config(const IniData& ini) {
    CoreConfig cfg;
    cfg.output_dir      = get_config_str_required(ini, "System", "OutPutDir");
    cfg.adjust_dir      = get_config_str_required(ini, "System", "AdjustDir");
    // Check that the output folders exist before computing (fail-fast). Output is written after a long
    // computation, so without stopping here a typo in a folder name would waste the entire run.
    // nc_create errors come from HDF5 and are hard to diagnose, so check explicitly.
    std::error_code ec_dir;
    if (!std::filesystem::is_directory(cfg.output_dir, ec_dir)) {
        handle_error("[System] OutPutDir does not exist or is not a directory.", cfg.output_dir);
    }
    if (!std::filesystem::is_directory(cfg.adjust_dir, ec_dir)) {
        handle_error("[System] AdjustDir does not exist or is not a directory.", cfg.adjust_dir);
    }
    cfg.nc_comp_level     = cfg_stoi(ini, "System", "NcCompLevel");
    cfg.utc_offset_hours  = cfg_stoi(ini, "System", "UtcOffsetHours");
    cfg.num_threads     = cfg_stoi(ini, "System", "NumThreads");
    if (cfg.num_threads <= 0)
        handle_error("[System] NumThreads must be > 0, got: " + std::to_string(cfg.num_threads), "");
    cfg.avg             = cfg_stoi(ini, "CalcWBGT", "avg");
    cfg.urban           = cfg_stoi(ini, "CalcWBGT", "urban");
    cfg.zspeed          = cfg_stod(ini, "CalcWBGT", "zspeed");
    cfg.dT              = cfg_stod(ini, "CalcWBGT", "dT");
    cfg.threshold_tair  = cfg_stod(ini, "CalcWBGT", "ThresholdTair");
    cfg.pgw             = cfg_stod(ini, "CalcWBGT", "PGW");
    cfg.lon_len_fine    = cfg_stoi(ini, "Resolution", "longitude");
    cfg.lat_len_fine    = cfg_stoi(ini, "Resolution", "latitude");
    cfg.lon_min = cfg_stod(ini, "Resolution", "lon_min");
    cfg.lon_max = cfg_stod(ini, "Resolution", "lon_max");
    cfg.lat_min = cfg_stod(ini, "Resolution", "lat_min");
    cfg.lat_max = cfg_stod(ini, "Resolution", "lat_max");
    // [Output] (optional): whether to save. Legacy behavior when keys are missing (save tair and sun, not shade).
    cfg.save_tair       = cfg_stoi_default(ini, "Output", "SaveTair",      1) != 0;
    cfg.save_wbgt_sun   = cfg_stoi_default(ini, "Output", "SaveWbgtSun",   1) != 0;
    cfg.save_wbgt_shade = cfg_stoi_default(ini, "Output", "SaveWbgtShade", 0) != 0;
    // Number of significant digits for quantization. 0 = disabled (lossless, legacy behavior). Negative values are treated as disabled.
    cfg.quantize_nsd    = cfg_stoi_default(ini, "Output", "QuantizeNSD", 0);
    if (cfg.quantize_nsd < 0) cfg.quantize_nsd = 0;
    // If no variable is selected for saving, nothing would remain after a long computation, so stop before computing
    if (!cfg.save_tair && !cfg.save_wbgt_sun && !cfg.save_wbgt_shade)
        handle_error("[Output] At least one of SaveTair / SaveWbgtSun / SaveWbgtShade must be 1.", "");

    if (cfg.lon_len_fine <= 0)
        handle_error("[Resolution] longitude must be > 0, got: " + std::to_string(cfg.lon_len_fine), "");
    if (cfg.lat_len_fine <= 0)
        handle_error("[Resolution] latitude must be > 0, got: " + std::to_string(cfg.lat_len_fine), "");
    if (cfg.lon_min < -180.0 || cfg.lon_min > 180.0)
        handle_error("[Resolution] lon_min out of range [-180, 180]: " + std::to_string(cfg.lon_min), "");
    if (cfg.lon_max < -180.0 || cfg.lon_max > 180.0)
        handle_error("[Resolution] lon_max out of range [-180, 180]: " + std::to_string(cfg.lon_max), "");
    if (cfg.lat_min < -90.0 || cfg.lat_min > 90.0)
        handle_error("[Resolution] lat_min out of range [-90, 90]: " + std::to_string(cfg.lat_min), "");
    if (cfg.lat_max < -90.0 || cfg.lat_max > 90.0)
        handle_error("[Resolution] lat_max out of range [-90, 90]: " + std::to_string(cfg.lat_max), "");
    if (cfg.lon_min >= cfg.lon_max)
        handle_error("[Resolution] lon_min must be < lon_max: " + std::to_string(cfg.lon_min) + " >= " + std::to_string(cfg.lon_max), "");
    if (cfg.lat_min >= cfg.lat_max)
        handle_error("[Resolution] lat_min must be < lat_max: " + std::to_string(cfg.lat_min) + " >= " + std::to_string(cfg.lat_max), "");

    cfg.landmask_filepath     = get_config_str_required(ini, "FeatureFile", "landmask");
    cfg.feature_list_filepath = get_config_str_required(ini, "FeatureFile", "FeatureList");
    cfg.model_tair   = get_config_str_required(ini, "XgbModel", "air_temperature");
    cfg.model_relhum = get_config_str_required(ini, "XgbModel", "relative_humidity");
    cfg.model_psurf  = get_config_str_required(ini, "XgbModel", "station_pressure");
    cfg.model_wind   = get_config_str_required(ini, "XgbModel", "wind_speed");
    // Apply the configured value to log timestamps immediately (it is also set in run_wbgt_pipeline,
    // but this also covers earlier logs such as log_core_config)
    g_log_utc_offset = cfg.utc_offset_hours;
    return cfg;
}

// Print the effective config values per section.
// Records at startup the values actually in effect, including keys filled with defaults ([Output], etc.).
void log_core_config(const CoreConfig& cfg, const string& extra_system) {
    config_line(std::format("[System] OutPutDir = {}", cfg.output_dir));
    config_line(std::format("[System] AdjustDir = {}", cfg.adjust_dir));
    config_line(std::format("[System] NcCompLevel = {}, NumThreads = {}, UtcOffsetHours = {}",
        cfg.nc_comp_level, cfg.num_threads, cfg.utc_offset_hours));
    if (!extra_system.empty()) {
        config_line(extra_system);
    }
    config_line(std::format("[CalcWBGT] avg = {}, urban = {}, zspeed = {}, dT = {}, ThresholdTair = {}, PGW = {}",
        cfg.avg, cfg.urban, cfg.zspeed, cfg.dT, cfg.threshold_tair, cfg.pgw));
    config_line(std::format("[Output] SaveTair = {}, SaveWbgtSun = {}, SaveWbgtShade = {}, QuantizeNSD = {}",
        (int)cfg.save_tair, (int)cfg.save_wbgt_sun, (int)cfg.save_wbgt_shade, cfg.quantize_nsd));
    config_line(std::format("[Resolution] longitude = {}, latitude = {}, lon = [{}, {}], lat = [{}, {}]",
        cfg.lon_len_fine, cfg.lat_len_fine, cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max));
    config_line(std::format("[FeatureFile] landmask = {}", cfg.landmask_filepath));
    config_line(std::format("[FeatureFile] FeatureList = {}", cfg.feature_list_filepath));
    config_line(std::format("[XgbModel] air_temperature = {}", cfg.model_tair));
    config_line(std::format("[XgbModel] relative_humidity = {}", cfg.model_relhum));
    config_line(std::format("[XgbModel] station_pressure = {}", cfg.model_psurf));
    config_line(std::format("[XgbModel] wind_speed = {}", cfg.model_wind));
}

struct AdjustData {
    std::vector<float> tair, relhum, psurf, wind;
    float fill_tair   = NC_FILL_FLOAT;
    float fill_relhum = NC_FILL_FLOAT;
    float fill_psurf  = NC_FILL_FLOAT;
    float fill_wind   = NC_FILL_FLOAT;
};

// Run XGBoost prediction and copy into dst[offset..] after checking the return value, null, and length.
// Without these checks, a failed Predict leaves out_result null (null dereference), or an unexpected
// out_len overruns the dst buffer (B-1).
static void xgb_predict_and_copy(BoosterHandle booster, DMatrixHandle dmatrix,
                                 std::vector<float>& dst, size_t offset, size_t expected_len,
                                 const char* name)
{
    bst_ulong out_len = 0;
    const float* out_result = nullptr;
    if (XGBoosterPredict(booster, dmatrix, 0, 0, 0, &out_len, &out_result) != 0) {
        handle_error(std::string("XGBoosterPredict failed for ") + name + ": " + XGBGetLastError(), "");
    }
    if (out_result == nullptr) {
        handle_error(std::string("XGBoosterPredict returned a null result for ") + name + ".", "");
    }
    if (out_len != expected_len) {
        handle_error(std::string("XGBoosterPredict returned unexpected length for ") + name +
                     " (got " + std::to_string(out_len) + ", expected " + std::to_string(expected_len) + ").", "");
    }
    std::copy(out_result, out_result + out_len, dst.begin() + offset);
}

static AdjustData compute_or_load_adjust(const CoreConfig& cfg, int year, int month,
                                   const std::vector<double>& vec_lon_fine,
                                   const std::vector<double>& vec_lat_fine)
{
    const string delimiter = "/";

    AdjustData adj;
    adj.tair.assign(  size_t(24) * cfg.lat_len_fine * cfg.lon_len_fine, NC_FILL_FLOAT);
    adj.relhum.assign(size_t(24) * cfg.lat_len_fine * cfg.lon_len_fine, NC_FILL_FLOAT);
    adj.psurf.assign( size_t(24) * cfg.lat_len_fine * cfg.lon_len_fine, NC_FILL_FLOAT);
    adj.wind.assign(  size_t(24) * cfg.lat_len_fine * cfg.lon_len_fine, NC_FILL_FLOAT);

    // Append domain information to the file name (a different domain is managed as a separate file)
    string domain_tag = "_lat_" + std::to_string(cfg.start_lat_idx) + "_" + std::to_string(cfg.start_lat_idx + cfg.lat_len_fine)
                      + "_lon_" + std::to_string(cfg.start_lon_idx) + "_" + std::to_string(cfg.start_lon_idx + cfg.lon_len_fine);
    string path_adj_tair   = cfg.adjust_dir + delimiter + "adjust_tair_"   + std::to_string(month) + domain_tag + ".nc";
    string path_adj_relhum = cfg.adjust_dir + delimiter + "adjust_relhum_" + std::to_string(month) + domain_tag + ".nc";
    string path_adj_psurf  = cfg.adjust_dir + delimiter + "adjust_psurf_"  + std::to_string(month) + domain_tag + ".nc";
    string path_adj_wind   = cfg.adjust_dir + delimiter + "adjust_wind_"   + std::to_string(month) + domain_tag + ".nc";

    bool adj_exist = true;
    adj_exist &= load_nc_2d_float(path_adj_tair,   "adjust_tair",   adj.tair,   &adj.fill_tair,   false);
    adj_exist &= load_nc_2d_float(path_adj_relhum, "adjust_relhum", adj.relhum, &adj.fill_relhum, false);
    adj_exist &= load_nc_2d_float(path_adj_psurf,  "adjust_psurf",  adj.psurf,  &adj.fill_psurf,  false);
    adj_exist &= load_nc_2d_float(path_adj_wind,   "adjust_wind",   adj.wind,   &adj.fill_wind,   false);

    if (!adj_exist) {
        information("Calculated adjustments do not exist. Calculate them.");

        std::map<string, string> feature_file_map;
        std::vector<string> feature_name_vec;
        if (load_txt_map(cfg.feature_list_filepath, feature_file_map, feature_name_vec) != 0) {
            handle_error("Feature list file not found or cannot be opened.", cfg.feature_list_filepath);
        }

        int normal_feature_num = 0;
        for (int i = 0; i < (int)feature_name_vec.size(); i++) {
            const string& fn = feature_name_vec[i];
            if (fn.starts_with("hour.local.dummy") || fn.starts_with("cza.past")) {
            } else if (fn.starts_with("temperature.monthly.mean") || fn.starts_with("rainfall.monthly.mean")) {
                feature_file_map[fn] = feature_file_map[fn] + "_" + std::to_string(month) + ".nc";
                normal_feature_num += 1;
            } else {
                normal_feature_num += 1;
            }
        }

        int total_feature_num = normal_feature_num + 24 + 24;
        config_line("[FeatureFile] feature -> source file:");
        show_map(feature_file_map, feature_name_vec);

        BoosterHandle booster_tair, booster_relhum, booster_psurf, booster_wind;
        DMatrixHandle dmatrix;
        std::vector<float> data_xgb(total_feature_num);
        XGDMatrixCreateFromMat(data_xgb.data(), 1, total_feature_num, -1, &dmatrix);

        if (XGBoosterCreate(&dmatrix, 1, &booster_tair) != 0 ||
            XGBoosterLoadModel(booster_tair, cfg.model_tair.c_str()) != 0)
            handle_error("failed to load tair model", cfg.model_tair);
        else information("tair model loaded.");

        if (XGBoosterCreate(&dmatrix, 1, &booster_relhum) != 0 ||
            XGBoosterLoadModel(booster_relhum, cfg.model_relhum.c_str()) != 0)
            handle_error("failed to load relhum model", cfg.model_relhum);
        else information("relhum model loaded.");

        if (XGBoosterCreate(&dmatrix, 1, &booster_psurf) != 0 ||
            XGBoosterLoadModel(booster_psurf, cfg.model_psurf.c_str()) != 0)
            handle_error("failed to load psurf model", cfg.model_psurf);
        else information("psurf model loaded.");

        if (XGBoosterCreate(&dmatrix, 1, &booster_wind) != 0 ||
            XGBoosterLoadModel(booster_wind, cfg.model_wind.c_str()) != 0)
            handle_error("failed to load wind model", cfg.model_wind);
        else information("wind model loaded.");

        std::vector<std::vector<float>> features_normal(normal_feature_num, std::vector<float>((size_t)cfg.lon_len_fine * cfg.lat_len_fine, NC_FILL_FLOAT));
        std::vector<std::vector<float>> features_cza_original(24, std::vector<float>((size_t)cfg.lon_len_fine * cfg.lat_len_fine, NC_FILL_FLOAT));
        std::vector<float> fillvalue_feature_vec(normal_feature_num);

        for (int i = 0; i < normal_feature_num; i++) {
            float fillvalue_temp;
            load_nc_2d_float(feature_file_map[feature_name_vec[i]], feature_name_vec[i], features_normal[i], &fillvalue_temp,
                             true, cfg.start_lat_idx, cfg.start_lon_idx, (size_t)cfg.lat_len_fine, (size_t)cfg.lon_len_fine);
            fillvalue_feature_vec[i] = fillvalue_temp;
        }

        information("solar position calculation start.");
        for (int hour = 0; hour < 24; hour++) {
            for (int lat_cnt_fine = 0; lat_cnt_fine < cfg.lat_len_fine; lat_cnt_fine++) {
                double lat_fine = vec_lat_fine[lat_cnt_fine];
#pragma omp parallel for num_threads(cfg.num_threads) schedule(dynamic, 1)
                for (int lon_cnt_fine = 0; lon_cnt_fine < cfg.lon_len_fine; lon_cnt_fine++) {
                    double lon_fine = vec_lon_fine[lon_cnt_fine];
                    double dday = 15 + (double)hour / 24.0;
                    double solar, cza, fdir;
                    calc_solar_parameters(year, month, dday, lat_fine, lon_fine, &solar, &cza, &fdir);
                    size_t idx = size_t(lat_cnt_fine) * size_t(cfg.lon_len_fine) + size_t(lon_cnt_fine);
                    features_cza_original[hour][idx] = (float)cza;
                }
            }
        }
        information("solar position calculated.");

        std::vector<float> features_all((size_t)total_feature_num * cfg.lon_len_fine * cfg.lat_len_fine, NC_FILL_FLOAT);
        for (int hour = 0; hour < 24; hour++) {
            arrange_features(features_normal, features_cza_original, vec_lon_fine, vec_lat_fine, hour, cfg.num_threads, features_all);

            XGDMatrixCreateFromMat(features_all.data(), (size_t)cfg.lon_len_fine * cfg.lat_len_fine, total_feature_num, -9999, &dmatrix);

            size_t offset = (size_t)hour * size_t(cfg.lat_len_fine) * size_t(cfg.lon_len_fine);
            size_t expected_len = (size_t)cfg.lat_len_fine * size_t(cfg.lon_len_fine);

            information("tair prediction started");
            xgb_predict_and_copy(booster_tair, dmatrix, adj.tair, offset, expected_len, "tair");
            information("tair prediction completed");

            information("relhum prediction started");
            xgb_predict_and_copy(booster_relhum, dmatrix, adj.relhum, offset, expected_len, "relhum");
            information("relhum prediction completed");

            information("psurf prediction started");
            xgb_predict_and_copy(booster_psurf, dmatrix, adj.psurf, offset, expected_len, "psurf");
            information("psurf prediction completed");

            information("wind prediction started");
            xgb_predict_and_copy(booster_wind, dmatrix, adj.wind, offset, expected_len, "wind");
            information("wind prediction completed");

            XGDMatrixFree(dmatrix);
        }

        std::vector<double> vec_hour = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
        save_as_ncdf(path_adj_tair,   adj.tair,   "adjust_tair",   vec_hour, vec_lat_fine, vec_lon_fine, 1);
        save_as_ncdf(path_adj_relhum, adj.relhum, "adjust_relhum", vec_hour, vec_lat_fine, vec_lon_fine, 1);
        save_as_ncdf(path_adj_psurf,  adj.psurf,  "adjust_psurf",  vec_hour, vec_lat_fine, vec_lon_fine, 1);
        save_as_ncdf(path_adj_wind,   adj.wind,   "adjust_wind",   vec_hour, vec_lat_fine, vec_lon_fine, 1);
        information("saved adjustment.");

        XGBoosterFree(booster_tair);
        XGBoosterFree(booster_relhum);
        XGBoosterFree(booster_psurf);
        XGBoosterFree(booster_wind);
    } else {
        information("Calculated adjustments already exist. Reuse them.");
    }

    return adj;
}

static void write_results(const CoreConfig& cfg, const string& year_month,
                   const std::vector<float>& tair_fine,
                   const std::vector<float>& wbgt_sun, const std::vector<float>& wbgt_shade,
                   const std::vector<double>& vec_time,
                   const std::vector<double>& vec_lat_fine,
                   const std::vector<double>& vec_lon_fine)
{
    const string delimiter = "/";
    string domain_tag = "_lat_" + std::to_string(cfg.start_lat_idx) + "_" + std::to_string(cfg.start_lat_idx + cfg.lat_len_fine)
                      + "_lon_" + std::to_string(cfg.start_lon_idx) + "_" + std::to_string(cfg.start_lon_idx + cfg.lon_len_fine);
    // ISIMIP3b output names are separated by GCM and scenario (to avoid collisions in the same output directory)
    const string& src = cfg.source_tag;
    // With PGW, insert pgw<value>_ into the file name
    string pgw_tag;
    if (cfg.pgw != 0.0) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << cfg.pgw;
        pgw_tag = "pgw" + oss.str() + "_";
    }
    // Build the WBGT (sun/shade) output path from the variant token
    auto wbgt_path = [&](const string& variant) {
        return cfg.output_dir + delimiter + "wbgt_" + variant + "_downscaled_" + pgw_tag + src + year_month + domain_tag + ".nc";
    };

    if (cfg.save_tair) {
        string path_tair = cfg.output_dir + delimiter + "tair_downscaled_" + src + year_month + domain_tag + ".nc";
        save_as_ncdf(path_tair, tair_fine, "tair", vec_time, vec_lat_fine, vec_lon_fine, cfg.nc_comp_level, cfg.quantize_nsd);
    }
    if (cfg.save_wbgt_sun) {
        save_as_ncdf(wbgt_path("sun"), wbgt_sun, "wbgt_sun", vec_time, vec_lat_fine, vec_lon_fine, cfg.nc_comp_level, cfg.quantize_nsd);
    }
    if (cfg.save_wbgt_shade) {
        save_as_ncdf(wbgt_path("shade"), wbgt_shade, "wbgt_shade", vec_time, vec_lat_fine, vec_lon_fine, cfg.nc_comp_level, cfg.quantize_nsd);
    }
}


// Main WBGT pipeline
int run_wbgt_pipeline(CoreConfig& cfg, int year, int month,
                      ClimateVar& tair, ClimateVar& psurf, ClimateVar& qair,
                      ClimateVar& swdown, ClimateVar& wind)
{
    g_log_utc_offset = cfg.utc_offset_hours;

    // Convert the month to a zero-padded string and join it with the year (used in output file names)
    // -------------------------------------------------------------------------------------------
    std::ostringstream year_month_s;
    year_month_s << year << std::setw(2) << std::setfill('0') << month;
    const string year_month = year_month_s.str();

    // Reference time
    // -------------------------------------------------------------------------------------------
    const sys_time<minutes> systemtime1900 = sys_days{ 1900y / 1 / 1 } + 0h + 0min;
    const sys_time<minutes> systemtime1970 = sys_days{ 1970y / 1 / 1 } + 0h + 0min;

    // Define the interpolated latitude/longitude (supports subsets via lon_min/max, lat_min/max)
    // -------------------------------------------------------------------------------------------
    std::vector<double> vec_lon_fine;
    std::vector<double> vec_lat_fine;
    {
        const double dlon_full = 360.0 / cfg.lon_len_fine;
        const double dlat_full = 180.0 / cfg.lat_len_fine;
        bool lon_started = false;
        for (int i = 0; i < cfg.lon_len_fine; i++) {
            double lon = -180.0 + (dlon_full * i) + (dlon_full / 2);
            if (lon >= cfg.lon_min && lon <= cfg.lon_max) {
                if (!lon_started) { cfg.start_lon_idx = (size_t)i; lon_started = true; }
                vec_lon_fine.push_back(lon);
            }
        }
        bool lat_started = false;
        for (int i = 0; i < cfg.lat_len_fine; i++) {
            double lat = -90.0 + (dlat_full * i) + (dlat_full / 2);
            if (lat >= cfg.lat_min && lat <= cfg.lat_max) {
                if (!lat_started) { cfg.start_lat_idx = (size_t)i; lat_started = true; }
                vec_lat_fine.push_back(lat);
            }
        }
    }
    if (vec_lon_fine.empty() || vec_lat_fine.empty()) {
        handle_error("No grid cells in the specified domain. Check lon_min/max and lat_min/max in config.", "");
    }
    cfg.lon_len_fine = (int)vec_lon_fine.size();
    cfg.lat_len_fine = (int)vec_lat_fine.size();
    information(std::format("Domain: lon [{:.4f}, {:.4f}] ({} cells), lat [{:.4f}, {:.4f}] ({} cells)",
        vec_lon_fine.front(), vec_lon_fine.back(), cfg.lon_len_fine,
        vec_lat_fine.front(), vec_lat_fine.back(), cfg.lat_len_fine));

    // Load bias corrections from the cache, or compute them with XGBoost
    // -------------------------------------------------------------------------------------------
    AdjustData adj = compute_or_load_adjust(cfg, year, month, vec_lon_fine, vec_lat_fine);

    // Load the land mask
    // -------------------------------------------------------------------------------------------
    float fill_value_landmask;
    std::vector<float> vec_landmask((size_t)cfg.lon_len_fine * cfg.lat_len_fine, NC_FILL_FLOAT);
    load_nc_2d_float(cfg.landmask_filepath, "LandUse.sum.Type1.0", vec_landmask, &fill_value_landmask,
                     true, cfg.start_lat_idx, cfg.start_lon_idx, (size_t)cfg.lat_len_fine, (size_t)cfg.lon_len_fine);

    // NetCDF C library check
    // -------------------------------------------------------------------------------------------
    information(std::format("Loading NetCDF files.(target = {}/{}, library version = {})", year, month, nc_inq_libvers()));

    const int time_len = tair.time_len;
    const int lat_len  = tair.lat_len;
    const int lon_len  = tair.lon_len;
    const bool timeunit1900    = tair.timeunit1900;
    const auto& vec_time       = tair.vec_time;
    const auto& vec_lat        = tair.vec_lat;
    const auto& vec_lon        = tair.vec_lon;
    const auto& vec_tair       = tair.data;
    const auto& vec_psurf      = psurf.data;
    const auto& vec_qair       = qair.data;
    const auto& vec_swdown     = swdown.data;
    const auto& vec_wind       = wind.data;
    const float fillvalue_tair   = tair.fillvalue;
    const float fillvalue_psurf  = psurf.fillvalue;
    const float fillvalue_qair   = qair.fillvalue;
    const float fillvalue_swdown = swdown.fillvalue;
    const float fillvalue_wind   = wind.fillvalue;

    // Warn if the output domain extends beyond the coverage of the input grid (up to the cell edges). Cells beyond it
    // remain missing due to the interpolation bounds check (e.g. ISIMIP3b lat coverage of about -56.25 to 83.75).
    // -------------------------------------------------------------------------------------------
    if (vec_lat.size() > 1 && vec_lon.size() > 1) {
        const double half_dlat = std::abs(vec_lat[1] - vec_lat[0]) / 2.0;
        const double half_dlon = std::abs(vec_lon[1] - vec_lon[0]) / 2.0;
        const double in_lat_min = std::min(vec_lat.front(), vec_lat.back()) - half_dlat;
        const double in_lat_max = std::max(vec_lat.front(), vec_lat.back()) + half_dlat;
        const double in_lon_min = std::min(vec_lon.front(), vec_lon.back()) - half_dlon;
        const double in_lon_max = std::max(vec_lon.front(), vec_lon.back()) + half_dlon;
        if (vec_lat_fine.front() < in_lat_min || vec_lat_fine.back() > in_lat_max) {
            information(std::format("Warning: output lat range [{:.4f}, {:.4f}] exceeds input coverage [{:.4f}, {:.4f}]. Cells outside will be left as missing values.",
                vec_lat_fine.front(), vec_lat_fine.back(), in_lat_min, in_lat_max));
        }
        if (vec_lon_fine.front() < in_lon_min || vec_lon_fine.back() > in_lon_max) {
            information(std::format("Warning: output lon range [{:.4f}, {:.4f}] exceeds input coverage [{:.4f}, {:.4f}]. Cells outside will be left as missing values.",
                vec_lon_fine.front(), vec_lon_fine.back(), in_lon_min, in_lon_max));
        }
    }

    // Start message
    // -------------------------------------------------------------------------------------------
    information(std::format("WBGT calculation for all time steps start.(target = {}/{}, time = {}, lat = {}, lon = {})", year, month, time_len, lat_len, lon_len));

    clock_t start = clock();

    // Result arrays (allocated only for variables flagged for saving, to reduce peak memory)
    // With _HAS_EXCEPTIONS=0, an allocation failure terminates immediately (0xC0000409) without a message.
    // Log the required size beforehand so the cause can be identified.
    const size_t fine_size = size_t(time_len) * size_t(cfg.lat_len_fine) * size_t(cfg.lon_len_fine);
    const int n_save_vars = (cfg.save_tair ? 1 : 0) + (cfg.save_wbgt_sun ? 1 : 0) + (cfg.save_wbgt_shade ? 1 : 0);
    information(std::format("Allocating result arrays: {:.1f} GiB ({} output variable(s) x {:.1f} GiB).",
        (double)fine_size * sizeof(float) * n_save_vars / (1024.0 * 1024.0 * 1024.0),
        n_save_vars,
        (double)fine_size * sizeof(float) / (1024.0 * 1024.0 * 1024.0)));
    std::vector<float> vec_tair_fine      (cfg.save_tair       ? fine_size : 0, NC_FILL_FLOAT);
    std::vector<float> vec_wbgt_sun_fine  (cfg.save_wbgt_sun   ? fine_size : 0, NC_FILL_FLOAT);
    std::vector<float> vec_wbgt_shade_fine(cfg.save_wbgt_shade ? fine_size : 0, NC_FILL_FLOAT);
    information("Result arrays allocated.");

    // WBGT calculation
    // -------------------------------------------------------------------------------------------
    const double kelvin = 273.15;

    // For computing lat/lon indices (derived from the input grid; supports descending latitude and partial domains)
    // The same formula also matches WFDE5 (ascending lat from -89.75, ascending lon from -179.75).
    const double lat0_in = vec_lat.empty() ? 0.0 : vec_lat[0];
    const double dlat_in = (vec_lat.size() > 1) ? (vec_lat[1] - vec_lat[0]) : 0.5;
    const double lon0_in = vec_lon.empty() ? 0.0 : vec_lon[0];
    const double dlon_in = (vec_lon.size() > 1) ? (vec_lon[1] - vec_lon[0]) : 0.5;

    // Time loop
    for (int time_cnt = 0; time_cnt < time_len; time_cnt++) {

        // Convert the time value to a convenient form and split it into day, hour, and minute
        sys_time<seconds> utc_time;
        if (timeunit1900) {
            utc_time = systemtime1900 + hours{ int(vec_time[time_cnt]) };
        }
        else {
            // Seconds since 1970 exceed INT_MAX in 2038, so use long long (B-3)
            utc_time = systemtime1970 + seconds{ (long long)vec_time[time_cnt] };
        }

        // Show the time being processed
        information(std::format("Calculation ({:%Y-%m-%d %H:%M})", utc_time));

        // Time for the WBGT calculation (decomposed directly with chrono; no exceptions)
        int arg_gmt = 0;
        auto utc_days = floor<days>(utc_time);
        year_month_day ymd{ utc_days };
        hh_mm_ss hms{ utc_time - utc_days };
        int arg_year = int(ymd.year());
        int arg_month = (int)(unsigned)ymd.month();
        int arg_day = (int)(unsigned)ymd.day();
        int arg_hour = (int)hms.hours().count();
        int arg_minute = (int)hms.minutes().count();

        clock_t start1 = clock();

        // Latitude loop
        for (int lat_cnt_fine = 0; lat_cnt_fine < cfg.lat_len_fine; lat_cnt_fine++) {
            double lat_fine = vec_lat_fine[lat_cnt_fine];

            // Compute the latitude index
            double index_lat = (lat_fine - lat0_in) / dlat_in;
            int lowerIndex_lat = int(std::floor(index_lat));  // Lower index
            int upperIndex_lat = lowerIndex_lat + 1;   // Upper index

#pragma omp parallel for num_threads(cfg.num_threads) schedule(dynamic, 1)
            for (int lon_cnt_fine = 0; lon_cnt_fine < cfg.lon_len_fine; lon_cnt_fine++) {

                double lon_fine = vec_lon_fine[lon_cnt_fine];

                size_t idx_landmask_fine = size_t(lat_cnt_fine) * size_t(cfg.lon_len_fine) + size_t(lon_cnt_fine);
                size_t vec_idx_fine = size_t(time_cnt) * size_t(cfg.lat_len_fine) * size_t(cfg.lon_len_fine) + size_t(lat_cnt_fine) * size_t(cfg.lon_len_fine) + size_t(lon_cnt_fine);

                double landmask = vec_landmask[idx_landmask_fine];

                // Skip non-land cells
                if (landmask > 99.0) {
                    continue;
                }

                // Compute the longitude index
                double index_lon = (lon_fine - lon0_in) / dlon_in;
                int lowerIndex_lon = int(std::floor(index_lon));  // Lower index
                int upperIndex_lon = lowerIndex_lon +1;   // Upper index

                size_t vec_idx_q11 = 0;
                size_t vec_idx_q12 = 0;
                size_t vec_idx_q21 = 0;
                size_t vec_idx_q22 = 0;

                // index for q11
                if (lowerIndex_lon >= 0 && lowerIndex_lat >= 0 && lowerIndex_lon < lon_len && lowerIndex_lat < lat_len &&
                    upperIndex_lon >= 0 && upperIndex_lat >= 0 && upperIndex_lon < lon_len && upperIndex_lat < lat_len)
                {
                    vec_idx_q11 = size_t(time_cnt) * size_t(lat_len * lon_len) + size_t(lowerIndex_lat) * size_t(lon_len) + size_t(lowerIndex_lon);
                    vec_idx_q12 = size_t(time_cnt) * size_t(lat_len * lon_len) + size_t(upperIndex_lat) * size_t(lon_len) + size_t(lowerIndex_lon);
                    vec_idx_q21 = size_t(time_cnt) * size_t(lat_len * lon_len) + size_t(lowerIndex_lat) * size_t(lon_len) + size_t(upperIndex_lon);
                    vec_idx_q22 = size_t(time_cnt) * size_t(lat_len * lon_len) + size_t(upperIndex_lat) * size_t(lon_len) + size_t(upperIndex_lon);
                }else{
                    continue;
                }

                // Parameters for bilinear interpolation
                double lon_x1 = vec_lon[lowerIndex_lon];
                double lon_x2 = vec_lon[upperIndex_lon];
                double lat_y1 = vec_lat[lowerIndex_lat];
                double lat_y2 = vec_lat[upperIndex_lat];
                double q11 = 0.0, q12 = 0.0, q21 = 0.0, q22 = 0.0;
                double tair_fine0, psurf_fine0, qair_fine0, swdown_fine0, wind_fine0, relhum_fine0, es_fine0;
                double tair_fine1, psurf_fine1, wind_fine1, swdown_fine1, relhum_fine1;

                // Tair
                q11 = vec_tair[vec_idx_q11];
                q12 = vec_tair[vec_idx_q12];
                q21 = vec_tair[vec_idx_q21];
                q22 = vec_tair[vec_idx_q22];
                tair_fine0 = getBilinearInterpolation(lon_fine, lat_fine, lon_x1, lon_x2, lat_y1, lat_y2, q11, q12, q21, q22, fillvalue_tair);

                // Skip if there is no value
                if (tair_fine0 == fillvalue_tair) {
                    continue;
                }

                // Psurf
                q11 = vec_psurf[vec_idx_q11];
                q12 = vec_psurf[vec_idx_q12];
                q21 = vec_psurf[vec_idx_q21];
                q22 = vec_psurf[vec_idx_q22];
                psurf_fine0 = getBilinearInterpolation(lon_fine, lat_fine, lon_x1, lon_x2, lat_y1, lat_y2, q11, q12, q21, q22, fillvalue_psurf);

                // Qair (WFDE5: specific humidity / ISIMIP3b: relative humidity %)
                q11 = vec_qair[vec_idx_q11];
                q12 = vec_qair[vec_idx_q12];
                q21 = vec_qair[vec_idx_q21];
                q22 = vec_qair[vec_idx_q22];
                qair_fine0 = getBilinearInterpolation(lon_fine, lat_fine, lon_x1, lon_x2, lat_y1, lat_y2, q11, q12, q21, q22, fillvalue_qair);

                // SWdown
                q11 = vec_swdown[vec_idx_q11];
                q12 = vec_swdown[vec_idx_q12];
                q21 = vec_swdown[vec_idx_q21];
                q22 = vec_swdown[vec_idx_q22];
                swdown_fine0 = getBilinearInterpolation(lon_fine, lat_fine, lon_x1, lon_x2, lat_y1, lat_y2, q11, q12, q21, q22, fillvalue_swdown);

                // Wind Speed
                q11 = vec_wind[vec_idx_q11];
                q12 = vec_wind[vec_idx_q12];
                q21 = vec_wind[vec_idx_q21];
                q22 = vec_wind[vec_idx_q22];
                wind_fine0 = getBilinearInterpolation(lon_fine, lat_fine, lon_x1, lon_x2, lat_y1, lat_y2, q11, q12, q21, q22, fillvalue_wind);

                size_t vec_idx_adj_fine = size_t(arg_hour) * size_t(cfg.lat_len_fine) * size_t(cfg.lon_len_fine) + size_t(lat_cnt_fine) * size_t(cfg.lon_len_fine) + size_t(lon_cnt_fine);

                double tair_adj   = adj.tair[vec_idx_adj_fine];
                double relhum_adj = adj.relhum[vec_idx_adj_fine];
                double psurf_adj  = adj.psurf[vec_idx_adj_fine];
                double wind_adj   = adj.wind[vec_idx_adj_fine];

                tair_fine1   = tair_fine0 - kelvin + tair_adj;
                wind_fine1   = wind_fine0 + wind_adj;
                swdown_fine1 = swdown_fine0;

                if (cfg.humidity_is_rh) {
                    // ISIMIP3b: the humidity slot holds relative humidity (%). Pressure is in hPa.
                    relhum_fine0 = qair_fine0;
                } else {
                    // WFDE5: compute relative humidity from specific humidity.
                    es_fine0     = 6.1094 * exp((17.625 * (tair_fine0 - kelvin)) / ((tair_fine0 - kelvin) + 243.04));
                    relhum_fine0 = 100 * ((psurf_fine0 / 100.0 * qair_fine0) / ((0.622 + 0.378 * qair_fine0) * es_fine0));
                }

                if (cfg.pressure_in_hpa) {
                    psurf_fine1  = psurf_fine0 + psurf_adj;
                } else {
                    psurf_fine1  = psurf_fine0 / 100.0 + psurf_adj;
                }

                relhum_fine1 = relhum_fine0 + relhum_adj;

                if (relhum_fine1 > 100.0) {
                    relhum_fine1 = 100.0;
                }
                else if (relhum_fine1 < 0.0) {
                    relhum_fine1 = 0.0;
                }

                if (wind_fine1 < 0.0) {
                    wind_fine1 = 0.0;
                }

                if (cfg.save_tair) {
                    vec_tair_fine[vec_idx_fine] = (float)tair_fine1;
                }

                // WBGT (sun): computed with actual solar radiation
                if (cfg.save_wbgt_sun) {
                    double est_speed = 0, Tg = 0, Tnwb = 0, Twbg = 0;
                    calc_wbgt(arg_year, arg_month, arg_day, arg_hour, arg_minute, arg_gmt,
                        cfg.avg, lat_fine, lon_fine, swdown_fine1, psurf_fine1, tair_fine1 + cfg.pgw, relhum_fine1, wind_fine1,
                        cfg.zspeed, cfg.dT, cfg.urban, &est_speed, &Tg, &Tnwb, NULL, &Twbg);
                    vec_wbgt_sun_fine[vec_idx_fine] = (float)Twbg;
                }

                // WBGT (shade): set only solar radiation to 0, removing the radiant heat load of direct, diffuse, and ground-reflected radiation
                if (cfg.save_wbgt_shade) {
                    double est_speed = 0, Tg = 0, Tnwb = 0, Twbg = 0;
                    calc_wbgt(arg_year, arg_month, arg_day, arg_hour, arg_minute, arg_gmt,
                        cfg.avg, lat_fine, lon_fine, 0.0, psurf_fine1, tair_fine1 + cfg.pgw, relhum_fine1, wind_fine1,
                        cfg.zspeed, cfg.dT, cfg.urban, &est_speed, &Tg, &Tnwb, NULL, &Twbg);
                    vec_wbgt_shade_fine[vec_idx_fine] = (float)Twbg;
                }
            }
        }

        clock_t calc_end1 = clock();
        information(std::format("WBGT calculation for 1 time step completed. (duration = {:.3f}sec.)", (double)(calc_end1 - start1) / CLOCKS_PER_SEC));
    }

    // Completion message
    // -------------------------------------------------------------------------------------------
    clock_t calc_end = clock();
    information(std::format("WBGT calculation for all time steps completed.(target = {}/{}, duration = {:.3f}sec.)", year, month, (double)(calc_end - start) / CLOCKS_PER_SEC));

    // Write the results
    information("Writing WBGT data start.");
    start = clock();

    write_results(cfg, year_month, vec_tair_fine, vec_wbgt_sun_fine, vec_wbgt_shade_fine, vec_time, vec_lat_fine, vec_lon_fine);

    calc_end = clock();
    information(std::format("Writing WBGT data completed. (duration = {:.3f}sec.)", (double)(calc_end - start) / CLOCKS_PER_SEC));

    return 0;
}
