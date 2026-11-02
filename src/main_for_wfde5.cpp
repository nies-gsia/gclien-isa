// Copyright (c) 2026 Jun'ya Takakura, National Institute for Environmental Studies (NIES)
// SPDX-License-Identifier: BSD-3-Clause
//
// main_for_wfde5.cpp
// Main for WFDE5 input. Performs WFDE5-specific NetCDF reading, config parsing, and argument handling,
// then calls the shared pipeline run_wbgt_pipeline().

#include <iostream>
#include <format>
#include <filesystem>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstring>

#include <netcdf.h>

#include "gclien_core.h"

using std::cout;
using std::endl;
using std::string;


// Read a WFDE5 NetCDF variable. When ref != nullptr, dimensions and coordinates are shared with ref.
ClimateVar load_wfde5_var(const string& filepath, const string& varname, const ClimateVar* ref = nullptr)
{
    const char time_str[] = "time", lat_str[] = "lat", lon_str[] = "lon";

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
    if (file_time_len == 0) handle_error("NetCDFfile invalid time dimlength.",      filepath);
    if (file_lat_len  == 0) handle_error("NetCDFfile invalid latitude dimlength.",  filepath);
    if (file_lon_len  == 0) handle_error("NetCDFfile invalid longitude dimlength.", filepath);

    int ref_time_len, ref_lat_len, ref_lon_len, time_ldiff;
    if (ref == nullptr) {
        ref_time_len = file_time_len;
        ref_lat_len  = file_lat_len;
        ref_lon_len  = file_lon_len;
        time_ldiff   = 0;
        result.time_len = ref_time_len;
        result.lat_len  = ref_lat_len;
        result.lon_len  = ref_lon_len;
    } else {
        ref_time_len = ref->time_len;
        ref_lat_len  = ref->lat_len;
        ref_lon_len  = ref->lon_len;
        if (file_time_len > ref_time_len)
            handle_error("Time length does not match Tair NetCDFfile.",      filepath);
        if (file_lat_len != ref_lat_len)
            handle_error("Latitude length does not match Tair NetCDFfile.",  filepath);
        if (file_lon_len != ref_lon_len)
            handle_error("Longitude length does not match Tair NetCDFfile.", filepath);
        time_ldiff = ref_time_len - file_time_len;
    }

    char var_name[NC_MAX_NAME];
    nc_type var_xtype;
    int var_ndims, var_dimids[NC_MAX_VAR_DIMS], var_natts;
    int data_var_id = -1;
    int time_var_id = -1, lat_var_id = -1, lon_var_id = -1;

    for (int i = 0; i < nvars; i++) {
        status = nc_inq_var(ncid, i, var_name, &var_xtype, &var_ndims, var_dimids, &var_natts);
        if (status != NC_NOERR) nc_handle_error(status, filepath);
        if (strcmp(var_name, varname.c_str()) == 0 && var_xtype == NC_FLOAT)
            data_var_id = i;
        if (ref == nullptr) {
            if      (strcmp(var_name, time_str) == 0 && var_xtype == NC_INT64)  time_var_id = i;
            else if (strcmp(var_name, lat_str)  == 0 && var_xtype == NC_DOUBLE) lat_var_id  = i;
            else if (strcmp(var_name, lon_str)  == 0 && var_xtype == NC_DOUBLE) lon_var_id  = i;
        }
    }
    if (data_var_id < 0) handle_error("NetCDFfile invalid variable dimid.", filepath);
    if (ref == nullptr) {
        if (time_var_id < 0) handle_error("NetCDFfile invalid time dimid.",     filepath);
        if (lat_var_id  < 0) handle_error("NetCDFfile invalid latitude dimid.", filepath);
        if (lon_var_id  < 0) handle_error("NetCDFfile invalid longitude dimid.", filepath);
    }

    if (nc_get_att_float(ncid, data_var_id, "_FillValue", &result.fillvalue)) {
        std::cerr << "Error getting _FillValue attribute" << std::endl;
        exit(EXIT_FAILURE);
    }

    result.data.assign((size_t)ref_time_len * ref_lat_len * ref_lon_len, NC_FILL_FLOAT);
    status = nc_get_var_float(ncid, data_var_id,
        result.data.data() + (size_t)time_ldiff * ref_lat_len * ref_lon_len);
    if (status != NC_NOERR) nc_handle_error(status, filepath);

    if (ref == nullptr) {
        result.vec_time.assign(ref_time_len, NC_FILL_DOUBLE);
        status = nc_get_var_double(ncid, time_var_id, result.vec_time.data());
        if (status != NC_NOERR) nc_handle_error(status, filepath);

        result.vec_lat.assign(ref_lat_len, NC_FILL_DOUBLE);
        status = nc_get_var_double(ncid, lat_var_id, result.vec_lat.data());
        if (status != NC_NOERR) nc_handle_error(status, filepath);

        result.vec_lon.assign(ref_lon_len, NC_FILL_DOUBLE);
        status = nc_get_var_double(ncid, lon_var_id, result.vec_lon.data());
        if (status != NC_NOERR) nc_handle_error(status, filepath);

        nc_type tmp_xtypep;
        size_t unit_strlen;
        status = nc_inq_att(ncid, time_var_id, "units", &tmp_xtypep, &unit_strlen);
        if (status != NC_NOERR) nc_handle_error(status, filepath);
        char* time_unit = new char[unit_strlen + 1];
        status = nc_get_att_text(ncid, time_var_id, "units", time_unit);
        if (status != NC_NOERR) nc_handle_error(status, filepath);
        time_unit[unit_strlen] = '\0';
        result.timeunit1900 = (strcmp(time_unit, "hours since 1900-01-01") == 0);
        delete[] time_unit;
    }

    status = nc_close(ncid);
    if (status != NC_NOERR) nc_handle_error(status, filepath);

    return result;
}


// Load the WFDE5 config (common part + WFDE5-specific keys)
struct WfdeConfig {
    CoreConfig core;
    string wfde5_dir, filename_str, file_version;
};

static WfdeConfig load_wfde5_config(const string& config_filepath) {
    IniData ini;
    if (!load_ini(config_filepath, ini)) {
        handle_error("Config file not found or cannot be opened.", config_filepath);
    }
    WfdeConfig wc;
    wc.core = load_core_config(ini);
    wc.wfde5_dir    = get_config_str_required(ini, "System", "Wfde5Dir");
    wc.filename_str = get_config_str(ini, "System", "FileNameStr");
    wc.file_version = get_config_str(ini, "System", "FileVersion");
    // WFDE5: pressure in Pa, humidity as specific humidity. Output name prefix is "wfde5_" (counterpart of gcm_scenario_ for ISIMIP3b).
    wc.core.pressure_in_hpa = false;
    wc.core.humidity_is_rh  = false;
    wc.core.source_tag      = "wfde5_";
    information("Config file: " + config_filepath);
    log_core_config(wc.core, std::format("[System] Wfde5Dir = {}, FileNameStr = {}, FileVersion = {}",
        wc.wfde5_dir, wc.filename_str, wc.file_version));
    return wc;
}


// Main calculation
static int calculation(int year, int month, string config_filepath0) {

    WfdeConfig wc = load_wfde5_config("./" + config_filepath0);

    const string delimiter = "/";

    // Convert the month to a zero-padded string and join it with the year (used in input file names)
    std::ostringstream year_month_s;
    year_month_s << year << std::setw(2) << std::setfill('0') << month;
    const string year_month = year_month_s.str();

    // NetCDF file paths
    string fp_psurf  = wc.wfde5_dir + delimiter + "PSurf"  + delimiter + "PSurf"  + wc.filename_str + year_month + wc.file_version + ".nc";
    string fp_qair   = wc.wfde5_dir + delimiter + "Qair"   + delimiter + "Qair"   + wc.filename_str + year_month + wc.file_version + ".nc";
    string fp_swdown = wc.wfde5_dir + delimiter + "SWdown" + delimiter + "SWdown" + wc.filename_str + year_month + wc.file_version + ".nc";
    string fp_tair   = wc.wfde5_dir + delimiter + "Tair"   + delimiter + "Tair"   + wc.filename_str + year_month + wc.file_version + ".nc";
    string fp_wind   = wc.wfde5_dir + delimiter + "Wind"   + delimiter + "Wind"   + wc.filename_str + year_month + wc.file_version + ".nc";

    // Read each climate file
    ClimateVar tair   = load_wfde5_var(fp_tair,   "Tair");
    ClimateVar psurf  = load_wfde5_var(fp_psurf,  "PSurf",  &tair);
    ClimateVar qair   = load_wfde5_var(fp_qair,   "Qair",   &tair);
    ClimateVar swdown = load_wfde5_var(fp_swdown, "SWdown", &tair);
    ClimateVar wind   = load_wfde5_var(fp_wind,   "Wind",   &tair);

    return run_wbgt_pipeline(wc.core, year, month, tair, psurf, qair, swdown, wind);
}


int main(int argc, char* argv[])
{
    information("gclien-isa-wfde5 (WFDE5 input) started.");

    // Move to the directory containing the executable
    {
        std::error_code ec;
        auto exe_path = std::filesystem::canonical(argv[0], ec);
        if (ec) handle_error("Failed to get executable path.", "");
        std::filesystem::current_path(exe_path.parent_path(), ec);
        if (ec) handle_error("Failed to change to executable directory.", "");
    }

    // Check the number of arguments: {year} {month} {config}
    if (argc != 4) {
        handle_error("Usage: gclien-isa-wfde5.exe {year} {month} {config}", "");
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

    // Range check of year and month
    if (arg_year < 1950 || arg_year > 2049) {
        handle_error("Year out of range (1950-2049): " + std::to_string(arg_year), "");
    }
    if (arg_month < 1 || arg_month > 12) {
        handle_error("Month out of range (1-12): " + std::to_string(arg_month), "");
    }

    // Run the calculation
    int status = calculation(arg_year, arg_month, argv[3]);
    if (status != 0) handle_error("Failed to process WBGT calculation.", "");

    return 0;
}
