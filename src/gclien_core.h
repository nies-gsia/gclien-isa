// Copyright (c) 2026 Jun'ya Takakura, National Institute for Environmental Studies (NIES)
// SPDX-License-Identifier: BSD-3-Clause
//
// gclien_core.h
// Shared module providing the processing common to WFDE5 and ISIMIP3b inputs (logging, INI parsing,
// NetCDF I/O, bias-correction computation, bilinear interpolation, WBGT calculation loop, and output).
//
// Format-specific processing (command-line arguments, format-specific config keys, reading input NetCDF)
// is implemented in main_for_wfde5.cpp / main_for_isimip3b.cpp, which pass the loaded meteorological
// data to run_wbgt_pipeline() as ClimateVar.

#pragma once

#include <string>
#include <vector>
#include <map>
#include <netcdf.h>

using std::string;

// ---- Logging (the UTC offset is overwritten with the config value in run_wbgt_pipeline) ----
extern int g_log_utc_offset;
void information(string msg);
void config_line(string msg);
[[noreturn]] void nc_handle_error(int status, string filepath);
[[noreturn]] void handle_error(string msg, string filepath);

// ---- Utilities ----
bool is_number(const string& str);
string get_environment_variable(string const& variable_name);
// Non-throwing parsers (true on success, false if not convertible). Unlike std::stoi/std::stod, they never throw.
bool parse_int(const std::string& s, int& out);
bool parse_double(const std::string& s, double& out);

// ---- INI parser ----
using IniData = std::map<std::string, std::map<std::string, std::string>>;
// Read a config file. Returns true on success, false if the file is missing or cannot be opened
// (no exceptions; the caller must check the return value and call handle_error).
bool load_ini(const std::string& filepath, IniData& data);
std::string get_config_str(const IniData& ini, const char* section, const char* key);
// Read a required string key. Exits immediately via handle_error if it is missing or empty.
std::string get_config_str_required(const IniData& ini, const char* section, const char* key);
int cfg_stoi(const IniData& ini, const char* section, const char* key);
int cfg_stoi_default(const IniData& ini, const char* section, const char* key, int def);
double cfg_stod(const IniData& ini, const char* section, const char* key);
double cfg_stod_default(const IniData& ini, const char* section, const char* key, double def);

// ---- Input meteorological variable (storage format shared by WFDE5 / ISIMIP3b) ----
struct ClimateVar {
    std::vector<float> data;
    float fillvalue = NC_FILL_FLOAT;
    std::vector<double> vec_time;
    std::vector<double> vec_lat;
    std::vector<double> vec_lon;
    int time_len = 0;
    int lat_len  = 0;
    int lon_len  = 0;
    bool timeunit1900 = false;
};

// ---- Common settings ----
struct CoreConfig {
    string output_dir, adjust_dir;
    string landmask_filepath, feature_list_filepath;
    string model_tair, model_relhum, model_psurf, model_wind;
    int nc_comp_level = 1, num_threads = 1, avg = 0, urban = 0;
    int lon_len_fine = 7200, lat_len_fine = 3600;
    int utc_offset_hours = 9;
    double zspeed = 10.0, dT = 0.0, threshold_tair = 999.0, pgw = 0.0;
    double lon_min = -180.0, lon_max = 180.0, lat_min = -90.0, lat_max = 90.0;
    // Set by run_wbgt_pipeline() at run time (not read from config)
    size_t start_lon_idx = 0, start_lat_idx = 0;
    // Input-format-dependent behavior
    bool pressure_in_hpa = false;  // true: input pressure in hPa (ISIMIP3b), false: Pa (WFDE5)
    bool humidity_is_rh  = false;  // true: humidity slot holds relative humidity % (ISIMIP3b), false: specific humidity (WFDE5)
    string source_tag;             // Output file name prefix (e.g. "GFDL-ESM4_ssp585_")
    // Whether to save each output variable ([Output] section; legacy behavior when keys are missing)
    bool save_tair       = true;   // Air temperature (tair)
    bool save_wbgt_sun   = true;   // WBGT (sun, actual solar radiation)
    bool save_wbgt_shade = false;  // WBGT (shade, approximated with zero solar radiation)
    // Quantization of output NetCDF ([Output] QuantizeNSD): number of significant digits nsd. 0 = disabled (lossless).
    // When > 0, GranularBitRound is applied to the final outputs (tair/wbgt_sun/wbgt_shade),
    // rounding off low-order bits before deflate (lossy, but files become much smaller; write time is unchanged).
    // Not applied to the adjust_*.nc cache (so that recomputed results are not affected).
    int quantize_nsd = 0;
};

// Read the common sections (part of [System] / [CalcWBGT] / [Resolution] / [FeatureFile] / [XgbModel]).
// Format-specific keys, behavior flags, and source_tag are set by the caller.
CoreConfig load_core_config(const IniData& ini);

// Print the effective config values per section (startup log, without timestamps).
// If extra_system is non-empty, it is printed as one line right after the common [System] line
// (each main passes its format-specific keys, e.g. Wfde5Dir / Isimip3bDir).
void log_core_config(const CoreConfig& cfg, const string& extra_system = "");

// Main WBGT pipeline (compute/load bias corrections -> bilinear interpolation -> WBGT calculation -> NetCDF output).
// The input meteorological data (tair/psurf/qair/swdown/wind) are loaded by the caller according to the format.
// qair holds specific humidity for WFDE5 and relative humidity (%) for ISIMIP3b (distinguished by cfg.humidity_is_rh).
int run_wbgt_pipeline(CoreConfig& cfg, int year, int month,
                      ClimateVar& tair, ClimateVar& psurf, ClimateVar& qair,
                      ClimateVar& swdown, ClimateVar& wind);
