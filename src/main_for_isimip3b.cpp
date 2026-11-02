// Copyright (c) 2026 Jun'ya Takakura, National Institute for Environmental Studies (NIES)
// SPDX-License-Identifier: BSD-3-Clause
//
// main_for_isimip3b.cpp
// Main for ISIMIP3b input (temporally downscaled using W5E5).
// Performs ISIMIP3b-specific NetCDF reading, config parsing, and argument handling,
// then calls the shared pipeline run_wbgt_pipeline().
//
// ISIMIP3b data characteristics:
//   - Path      : {Isimip3bDir}/{scenario}/{gcm}/{lowercase gcm}_r1i1p1f{N}_w5e5_{scenario}_{var}_global_hourly_{start}_{end}.nc
//   - Variables : tas(K) ps(hPa) hurs(relative humidity %) rsds(W/m2) sfcwind(m/s) (pr/rlds are not used for WBGT)
//   - Storage   : ushort + scale_factor (packed), _FillValue=65535
//   - Dimensions: time (multiple years, hourly) x lat (descending, partial coverage) x lon
//   - Time units: "hours since {chunk start year}-01-01" (proleptic_gregorian)
// Each file contains multiple years, so only the time slice of the target year-month is read.

#include <iostream>
#include <format>
#include <filesystem>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <chrono>

#include <netcdf.h>

#include "gclien_core.h"

using std::cout;
using std::endl;
using std::string;
using namespace std::chrono;
using namespace std::chrono_literals;


// Find the chunk file containing the target year in the GCM directory. Returns its start year in chunk_start_year.
static string find_isimip3b_file(const string& gcm_dir, const string& gcm_lower,
                                 const string& scenario, const string& varname,
                                 int year, int& chunk_start_year)
{
    namespace fs = std::filesystem;
    // The ensemble member identifier (r1i1p1f1 / r1i1p1f2, etc.) differs by GCM, so it is not hard-coded;
    // instead, match head ("{gcm}_r1i1p1f") and mid ("_w5e5_{scenario}_{var}_global_hourly_") with only
    // digits (the number after f) in between (e.g. UKESM1-0-LL uses f2).
    const string head = gcm_lower + "_r1i1p1f";
    const string mid = string("_w5e5_") + scenario + "_" + varname + "_global_hourly_";
    std::error_code ec;
    if (!fs::is_directory(gcm_dir, ec))
        handle_error("ISIMIP3b directory not found.", gcm_dir);

    for (const auto& entry : fs::directory_iterator(gcm_dir, ec)) {
        const string fname = entry.path().filename().string();
        if (fname.rfind(head, 0) != 0) continue;                        // Prefix match ({gcm}_r1i1p1f)
        auto mpos = fname.find(mid, head.size());
        if (mpos == string::npos) continue;                             // Middle part match
        // Between head and mid there must be the ensemble number (digits only, at least one)
        if (mpos <= head.size()) continue;
        bool member_ok = true;
        for (size_t i = head.size(); i < mpos; ++i)
            if (!std::isdigit((unsigned char)fname[i])) { member_ok = false; break; }
        if (!member_ok) continue;
        if (fname.size() < 3 || fname.compare(fname.size() - 3, 3, ".nc") != 0) continue;
        // Parse the trailing "{start}_{end}.nc"
        string years = fname.substr(mpos + mid.size());
        years = years.substr(0, years.size() - 3);                      // Strip ".nc"
        auto us = years.find('_');
        if (us == string::npos) continue;
        int ys = 0, ye = 0;
        if (!parse_int(years.substr(0, us), ys)) continue;
        if (!parse_int(years.substr(us + 1), ye)) continue;
        if (year >= ys && year <= ye) {
            chunk_start_year = ys;
            return entry.path().string();
        }
    }
    handle_error("No ISIMIP3b file covering year " + std::to_string(year) +
                 " for variable '" + varname + "'", gcm_dir);
    return "";  // Not reached
}


// Read a text attribute (empty string if the attribute is missing or unreadable; no exceptions)
static string get_att_text(int ncid, int varid, const char* att_name)
{
    size_t len = 0;
    if (nc_inq_attlen(ncid, varid, att_name, &len) != NC_NOERR || len == 0) return "";
    std::vector<char> buf(len + 1, '\0');
    if (nc_get_att_text(ncid, varid, att_name, buf.data()) != NC_NOERR) return "";
    string s(buf.data(), len);
    while (!s.empty() && (s.back() == '\0' || s.back() == ' ')) s.pop_back();
    return s;
}


// Read the target year-month slice of ISIMIP3b data and return it as a ClimateVar.
// When ref != nullptr, dimensions and coordinates are shared with ref (pass tas as the reference).
static ClimateVar load_isimip3b_var(const string& isimip3b_dir, const string& scenario, const string& gcm,
                                    const string& varname, int year, int month, const ClimateVar* ref = nullptr)
{
    const char time_str[] = "time", lat_str[] = "lat", lon_str[] = "lon";

    // Lowercase the GCM name (for file names)
    string gcm_lower = gcm;
    std::transform(gcm_lower.begin(), gcm_lower.end(), gcm_lower.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });

    const string gcm_dir = isimip3b_dir + "/" + scenario + "/" + gcm;
    int chunk_start_year = 0;
    string filepath = find_isimip3b_file(gcm_dir, gcm_lower, scenario, varname, year, chunk_start_year);

    ClimateVar result;
    int status, ncid;
    status = nc_open(filepath.c_str(), NC_NOWRITE, &ncid);
    if (status != NC_NOERR) nc_handle_error(status, filepath);

    int ndims, nvars, ngatts, unlimdimid;
    status = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
    if (status != NC_NOERR) nc_handle_error(status, filepath);

    int file_time_len = 0, file_lat_len = 0, file_lon_len = 0;
    for (int i = 0; i < ndims; i++) {
        char dim_name[NC_MAX_NAME];
        size_t dim_len;
        status = nc_inq_dim(ncid, i, dim_name, &dim_len);
        if (status != NC_NOERR) nc_handle_error(status, filepath);
        if      (strcmp(dim_name, time_str) == 0) file_time_len = static_cast<int>(dim_len);
        else if (strcmp(dim_name, lat_str)  == 0) file_lat_len  = static_cast<int>(dim_len);
        else if (strcmp(dim_name, lon_str)  == 0) file_lon_len  = static_cast<int>(dim_len);
    }
    if (file_time_len == 0) handle_error("ISIMIP3b invalid time dimlength.",      filepath);
    if (file_lat_len  == 0) handle_error("ISIMIP3b invalid latitude dimlength.",  filepath);
    if (file_lon_len  == 0) handle_error("ISIMIP3b invalid longitude dimlength.", filepath);

    // Compute the time offset and number of hours of the target month (from Jan 1 00:00 of the chunk start year, proleptic_gregorian)
    const sys_days chunk_base  = sys_days{ std::chrono::year{chunk_start_year} / 1 / 1 };
    const sys_days month_start = sys_days{ std::chrono::year{year} / month / 1 };
    const int ny = (month == 12) ? year + 1 : year;
    const int nm = (month == 12) ? 1 : month + 1;
    const sys_days month_end   = sys_days{ std::chrono::year{ny} / nm / 1 };
    const long long hour_offset    = duration_cast<hours>(month_start - chunk_base).count();
    const long long hours_in_month = duration_cast<hours>(month_end - month_start).count();

    if (hour_offset < 0 || hour_offset + hours_in_month > file_time_len)
        handle_error("ISIMIP3b time range out of file bounds for " +
                     std::to_string(year) + "-" + std::to_string(month), filepath);

    if (ref != nullptr) {
        if (file_lat_len != ref->lat_len)       handle_error("ISIMIP3b latitude length mismatch.",  filepath);
        if (file_lon_len != ref->lon_len)       handle_error("ISIMIP3b longitude length mismatch.", filepath);
        if ((int)hours_in_month != ref->time_len) handle_error("ISIMIP3b time length mismatch.",    filepath);
    }
    result.time_len = (int)hours_in_month;
    result.lat_len  = file_lat_len;
    result.lon_len  = file_lon_len;

    // Variable IDs / coordinate IDs
    char var_name[NC_MAX_NAME];
    nc_type var_xtype;
    int var_ndims, var_dimids[NC_MAX_VAR_DIMS], var_natts;
    int data_var_id = -1, lat_var_id = -1, lon_var_id = -1, time_var_id = -1;
    for (int i = 0; i < nvars; i++) {
        status = nc_inq_var(ncid, i, var_name, &var_xtype, &var_ndims, var_dimids, &var_natts);
        if (status != NC_NOERR) nc_handle_error(status, filepath);
        if      (strcmp(var_name, varname.c_str()) == 0) data_var_id = i;
        else if (strcmp(var_name, lat_str) == 0)         lat_var_id  = i;
        else if (strcmp(var_name, lon_str) == 0)         lon_var_id  = i;
        else if (strcmp(var_name, time_str) == 0)        time_var_id = i;
    }
    if (data_var_id < 0) handle_error("ISIMIP3b variable not found: " + varname, filepath);

    // Check that the time calendar and units match the assumptions of the month-offset calculation (B-4).
    // Assumptions: units = "hours since {chunk start year}-01-01 ..." and a Gregorian-type calendar.
    // With noleap / 360_day files, the slice position would shift while still passing the bounds check,
    // silently reading data at a different time, so stop immediately if the assumptions do not hold.
    if (time_var_id >= 0) {
        const string calendar = get_att_text(ncid, time_var_id, "calendar");
        if (!calendar.empty() && calendar != "proleptic_gregorian" &&
            calendar != "gregorian" && calendar != "standard") {
            handle_error("ISIMIP3b unsupported time calendar '" + calendar +
                         "' (expected proleptic_gregorian).", filepath);
        }
        const string units = get_att_text(ncid, time_var_id, "units");
        const string expected_prefix = "hours since " + std::to_string(chunk_start_year) + "-01-01";
        if (!units.empty() && units.rfind(expected_prefix, 0) != 0) {
            handle_error("ISIMIP3b unexpected time units '" + units +
                         "' (expected '" + expected_prefix + " ...').", filepath);
        }
    }

    // Unpacking parameters (scale=1, offset=0 if the attributes are absent)
    double scale_factor = 1.0, add_offset = 0.0;
    nc_get_att_double(ncid, data_var_id, "scale_factor", &scale_factor);
    nc_get_att_double(ncid, data_var_id, "add_offset",   &add_offset);
    unsigned short fill_raw = 65535;
    nc_get_att_ushort(ncid, data_var_id, "_FillValue", &fill_raw);
    result.fillvalue = NC_FILL_FLOAT;

    // Read only the monthly slice (ushort) -> apply scale/offset and expand to float
    size_t start[3] = { (size_t)hour_offset, 0, 0 };
    size_t count[3] = { (size_t)hours_in_month, (size_t)file_lat_len, (size_t)file_lon_len };
    size_t n = (size_t)hours_in_month * (size_t)file_lat_len * (size_t)file_lon_len;
    std::vector<unsigned short> raw(n);
    status = nc_get_vara_ushort(ncid, data_var_id, start, count, raw.data());
    if (status != NC_NOERR) nc_handle_error(status, filepath);

    result.data.assign(n, NC_FILL_FLOAT);
    for (size_t i = 0; i < n; i++)
        result.data[i] = (raw[i] == fill_raw) ? NC_FILL_FLOAT : (float)(raw[i] * scale_factor + add_offset);

    if (ref == nullptr) {
        result.vec_lat.assign(file_lat_len, NC_FILL_DOUBLE);
        status = nc_get_var_double(ncid, lat_var_id, result.vec_lat.data());
        if (status != NC_NOERR) nc_handle_error(status, filepath);

        result.vec_lon.assign(file_lon_len, NC_FILL_DOUBLE);
        status = nc_get_var_double(ncid, lon_var_id, result.vec_lon.data());
        if (status != NC_NOERR) nc_handle_error(status, filepath);

        // Store time converted to "hours since 1900-01-01" (to reuse the downstream time conversion)
        const sys_days epoch1900 = sys_days{ std::chrono::year{1900} / 1 / 1 };
        const long long h1900 = duration_cast<hours>(month_start - epoch1900).count();
        result.vec_time.assign((size_t)hours_in_month, 0.0);
        for (long long i = 0; i < hours_in_month; i++)
            result.vec_time[(size_t)i] = (double)(h1900 + i);
        result.timeunit1900 = true;
    }

    status = nc_close(ncid);
    if (status != NC_NOERR) nc_handle_error(status, filepath);

    return result;
}


// Load the ISIMIP3b config (common part + ISIMIP3b-specific keys)
struct IsimipConfig {
    CoreConfig core;
    string isimip3b_dir;
};

static IsimipConfig load_isimip3b_config(const string& config_filepath) {
    IniData ini;
    if (!load_ini(config_filepath, ini)) {
        handle_error("Config file not found or cannot be opened.", config_filepath);
    }
    IsimipConfig ic;
    ic.core = load_core_config(ini);
    ic.isimip3b_dir = get_config_str(ini, "System", "Isimip3bDir");
    if (ic.isimip3b_dir.empty())
        handle_error("[System] Isimip3bDir is required for ISIMIP3b mode.", config_filepath);
    // ISIMIP3b: pressure in hPa, humidity as relative humidity.
    ic.core.pressure_in_hpa = true;
    ic.core.humidity_is_rh  = true;
    information("Config file: " + config_filepath);
    log_core_config(ic.core, std::format("[System] Isimip3bDir = {}", ic.isimip3b_dir));
    return ic;
}


// Main calculation
static int calculation(int year, int month, string config_filepath0,
                       const string& scenario, const string& gcm) {

    IsimipConfig ic = load_isimip3b_config("./" + config_filepath0);
    // Output file name prefix (separated by GCM and scenario)
    ic.core.source_tag = gcm + "_" + scenario + "_";

    // Read each climate file (dimensions and coordinates are shared, using tas as the reference)
    //   tas->Tair, ps->PSurf, hurs->Qair (humidity slot, relative humidity %), rsds->SWdown, sfcwind->Wind
    ClimateVar tair   = load_isimip3b_var(ic.isimip3b_dir, scenario, gcm, "tas",     year, month);
    ClimateVar psurf  = load_isimip3b_var(ic.isimip3b_dir, scenario, gcm, "ps",      year, month, &tair);
    ClimateVar qair   = load_isimip3b_var(ic.isimip3b_dir, scenario, gcm, "hurs",    year, month, &tair);
    ClimateVar swdown = load_isimip3b_var(ic.isimip3b_dir, scenario, gcm, "rsds",    year, month, &tair);
    ClimateVar wind   = load_isimip3b_var(ic.isimip3b_dir, scenario, gcm, "sfcwind", year, month, &tair);

    return run_wbgt_pipeline(ic.core, year, month, tair, psurf, qair, swdown, wind);
}


int main(int argc, char* argv[])
{
    information("gclien-isa-isimip (ISIMIP3b input) started.");

    // Move to the directory containing the executable
    {
        std::error_code ec;
        auto exe_path = std::filesystem::canonical(argv[0], ec);
        if (ec) handle_error("Failed to get executable path.", "");
        std::filesystem::current_path(exe_path.parent_path(), ec);
        if (ec) handle_error("Failed to change to executable directory.", "");
    }

    // Check the number of arguments: {year} {month} {config} {scenario} {gcm}
    if (argc != 6) {
        handle_error("Usage: gclien-isa-isimip.exe {year} {month} {config} {scenario} {gcm}", "");
    }

    // Argument 1, year: numeric check
    if (!is_number(argv[1])) {
        handle_error("Set a numeric value for argument(1) year.", "");
    }

    // Argument 2, month: numeric check
    if (!is_number(argv[2])) {
        handle_error("Set a numeric value for argument(2) month.", "");
    }

    // Convert year and month to integers
    int arg_year = 0, arg_month = 0;
    if (!parse_int(argv[1], arg_year)) {
        handle_error("Set a valid integer for argument(1) year.", "");
    }
    if (!parse_int(argv[2], arg_month)) {
        handle_error("Set a valid integer for argument(2) month.", "");
    }

    const string arg_scenario = argv[4];
    const string arg_gcm      = argv[5];

    // Validate arguments 4 and 5, scenario / gcm (A-5)
    // A typo would otherwise lead to indirect errors such as "directory not found", so stop early by
    // requiring an exact match with known values (case-sensitive, since they are used directly in paths).
    {
        static const char* const valid_scenarios[] = { "historical", "ssp126", "ssp370", "ssp585" };
        static const char* const valid_gcms[] = { "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL" };
        bool scenario_ok = false;
        for (const char* s : valid_scenarios) if (arg_scenario == s) { scenario_ok = true; break; }
        if (!scenario_ok) {
            handle_error("Unknown scenario for argument(4): '" + arg_scenario +
                         "'. Valid values: historical, ssp126, ssp370, ssp585", "");
        }
        bool gcm_ok = false;
        for (const char* g : valid_gcms) if (arg_gcm == g) { gcm_ok = true; break; }
        if (!gcm_ok) {
            handle_error("Unknown GCM for argument(5): '" + arg_gcm +
                         "'. Valid values (case-sensitive): GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL", "");
        }
    }

    // Range check of year and month (ISIMIP3b: historical 1850-2014 / ssp 2015-2100)
    if (arg_year < 1850 || arg_year > 2100) {
        handle_error("Year out of range (1850-2100): " + std::to_string(arg_year), "");
    }
    if (arg_month < 1 || arg_month > 12) {
        handle_error("Month out of range (1-12): " + std::to_string(arg_month), "");
    }

    // Run the calculation
    int status = calculation(arg_year, arg_month, argv[3], arg_scenario, arg_gcm);
    if (status != 0) handle_error("Failed to process WBGT calculation.", "");

    return 0;
}
