# gclien-isa

gclien-isa downscales hourly 0.5° meteorological forcing data to a 0.05° grid
and computes air temperature and the wet-bulb globe temperature (WBGT) on that
grid. The coarse fields are bilinearly interpolated, corrected with XGBoost
models that use high-resolution static features (elevation, land use, etc.),
and passed to the WBGT model of Liljegren et al. (2008). Results are written as
NetCDF files.

Two input datasets are supported, each with its own executable:

| Executable | Input data | Arguments |
|---|---|---|
| `gclien-isa-wfde5.exe` | WFDE5 v2.1 | `{year} {month} {config}` |
| `gclien-isa-isimip.exe` | ISIMIP3b (bias-adjusted with W5E5) | `{year} {month} {config} {scenario} {gcm}` |

Both executables share the downstream processing (bias correction,
interpolation, WBGT calculation, and output) in `src/gclien_core.cpp`.

## Status

As of November 2, 2026, the algorithms and code in this repository have not
been peer-reviewed. They should be regarded as preliminary and may change
without notice. Please keep this in mind when using the code or its outputs.

## Processing overview

1. Read the meteorological variables of the target month from the input NetCDF
   files.
   - WFDE5: `Tair` (K), `PSurf` (Pa), `Qair` (specific humidity), `SWdown`
     (W/m²), `Wind` (m/s); one file per variable and month.
   - ISIMIP3b: `tas` (K), `ps` (hPa), `hurs` (relative humidity, %), `rsds`
     (W/m²), `sfcwind` (m/s); only the target month is read from the
     multi-year files.
2. Compute bias corrections for air temperature, relative humidity, surface
   pressure, and wind speed with the XGBoost models (cached in `AdjustDir`).
3. For each hour and each land cell of the 0.05° grid (parallelized with
   OpenMP):
   - bilinearly interpolate the 0.5° inputs and apply the bias corrections;
   - convert humidity and pressure to relative humidity (%) and hPa;
   - compute WBGT with `calc_wbgt()` (`src/wbgtn.cpp`). WBGT in the sun uses
     the actual solar radiation; WBGT in the shade is approximated by setting
     solar radiation to 0.
4. Write the selected variables to NetCDF.

## Requirements (Windows)

| Component | Version used for development | Notes |
|---|---|---|
| Windows | 11 (x64) | |
| Visual Studio | 2022 or 2026 | "Desktop development with C++" workload; provides the MSVC headers and libraries required by icx, and Ninja |
| Intel oneAPI DPC++/C++ Compiler (`icx`) | 2026.0 | C++20 and OpenMP; oneAPI 2026.0 or later is required for Visual Studio 2026 |
| CMake | 3.20 or later | a version that provides the `Visual Studio 18 2026` generator is needed to build XGBoost with Visual Studio 2026 |
| Ninja | bundled with Visual Studio | |
| XGBoost (C API) | 3.3.0, CPU build | built from source |
| netCDF-C | 4.10.0 | Unidata Windows installer (includes HDF5 and zlib) |

Memory: the result array takes about 72 GiB per output variable for the whole
globe for one month (744 × 3600 × 7200 × 4 bytes). With 128 GB of RAM, compute
global runs one variable at a time (see `[Output]` below). Regional runs (e.g.
Japan) need only a few GiB.

## Setting up the build environment (Windows)

### 1. Visual Studio

Install Visual Studio 2026 or 2022 (Community edition is sufficient) with the
**Desktop development with C++** workload. This also installs CMake and Ninja
under the Visual Studio directory. Note that the installation directory is
`C:\Program Files\Microsoft Visual Studio\18\<edition>` for Visual Studio 2026
and `C:\Program Files\Microsoft Visual Studio\2022\<edition>` for Visual
Studio 2022. The build has been confirmed with both versions.

### 2. Intel oneAPI DPC++/C++ Compiler

Install the Intel oneAPI DPC++/C++ Compiler (stand-alone or as part of the
Intel oneAPI Base Toolkit). The default location is
`C:\Program Files (x86)\Intel\oneAPI\compiler\<version>\`.

### 3. CMake

If you prefer a stand-alone CMake to the one bundled with Visual Studio,
install it from <https://cmake.org/download/> (default location
`C:\Program Files\CMake`).

### 4. netCDF-C

Download the 64-bit Windows installer of netCDF-C 4.10.0 from Unidata
(<https://downloads.unidata.ucar.edu/netcdf/>) and install it with the default
settings (`C:\Program Files\netCDF 4.10.0`). The installation contains
`include\netcdf.h`, `lib\netcdf.lib`, and the runtime DLLs (`netcdf.dll`,
`hdf5.dll`, `hdf5_hl.dll`, `zlib1.dll`) in `bin\`.

### 5. XGBoost

Build the XGBoost shared library from source. In a
"Developer Command Prompt" of your Visual Studio version:

```bat
git clone --recursive https://github.com/dmlc/xgboost.git
cd xgboost
git checkout <release tag, e.g. the one for 3.3.0>
git submodule update --init --recursive
rem Visual Studio 2026:
cmake -S . -B build -G "Visual Studio 18 2026" -A x64
rem Visual Studio 2022:
rem cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

After the build, the XGBoost root directory must contain
`include\xgboost\c_api.h` and `lib\xgboost.lib` / `lib\xgboost.dll`.

### 6. Edit `build.bat`

Open `build.bat` and set the paths at the top to match your installation:

| Variable | Meaning |
|---|---|
| `XGBOOST_DIR` | XGBoost root directory (step 5) |
| `NETCDF_DIR` | netCDF-C installation directory (step 4) |
| `ONEAPI_BIN`, `ONEAPI_LIB` | `bin` and `lib` directories of the oneAPI compiler |
| `CMAKE_BIN` | directory containing `cmake.exe` |
| `NINJA_BIN` | directory containing `ninja.exe` |

Also update the path of the compiler's `env\vars.bat` in the `call` line if
your compiler version or location differs.

`build.bat` is set up for Visual Studio 2022. For Visual Studio 2026, change
the following two lines (`18` is the installation directory of Visual Studio
2026):

```bat
set "NINJA_BIN=C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja"
call "C:\Program Files (x86)\Intel\oneAPI\compiler\2026.0\env\vars.bat" intel64 vs2026
```

### 7. Build

Run the following from a Command Prompt in the repository root:

```bat
build.bat
```

`build.bat` initializes the compiler environment, configures the project with
CMake (Ninja generator, Release), builds both executables into `bin\`, and
copies the required runtime DLLs (oneAPI, XGBoost, netCDF/HDF5/zlib) into
`bin\`.

To call CMake directly instead:

```bat
cmake -S . -B build -G "Ninja" -DCMAKE_CXX_COMPILER=icx -DCMAKE_BUILD_TYPE=Release ^
    -DXGBOOST_DIR="<path>" -DNETCDF_DIR="<path>"
cmake --build build
```

In that case, make sure that `xgboost.dll`, `netcdf.dll` and its dependencies,
and the oneAPI runtime DLLs can be found at run time (in `bin\` or on `PATH`).

## Input data

### Meteorological forcing

**WFDE5 v2.1** (Cucchi et al., 2020), available from the Copernicus Climate
Data Store (<https://doi.org/10.24381/cds.20d54e34>).

> **Note:** WFDE5 v3.1 is now available from the same page, but the software has
> been tested only with v2.1.

The executable reads

```
{Wfde5Dir}/{Var}/{Var}{FileNameStr}{YYYYMM}{FileVersion}.nc
```

for `Var` = `Tair`, `PSurf`, `Qair`, `SWdown`, and `Wind`, e.g.
`E:/WFDE5_v2.1/Tair/Tair_WFDE5_CRU_201008_v2.1.nc`.

**ISIMIP3b hourly data** (Bechtold et al., 2026), available from the ISIMIP
repository (<https://doi.org/10.48364/ISIMIP.170328>). Note that the standard
ISIMIP3b climate forcing is daily; this software requires the hourly dataset,
which was produced by temporally disaggregating the daily ISIMIP3b data.
The executable reads

```
{Isimip3bDir}/{scenario}/{gcm}/{gcm}_r1i1p1f{N}_w5e5_{scenario}_{var}_global_hourly_{start}_{end}.nc
```

(`{gcm}` in lower case in the file name) for `var` = `tas`, `ps`, `hurs`,
`rsds`, and `sfcwind`.

- `scenario`: `historical` (1850–2014), `ssp126`, `ssp370`, `ssp585`
  (2015–2100)
- `gcm`: `GFDL-ESM4`, `IPSL-CM6A-LR`, `MPI-ESM1-2-HR`, `MRI-ESM2-0`,
  `UKESM1-0-LL`

Both arguments are case-sensitive.

### Static features and XGBoost models

The bias correction requires the static feature files (0.05° NetCDF) and the
trained XGBoost models:

- `input/input_features.txt`: tab-separated list of feature names and the paths
  of their NetCDF files (paths relative to `bin\`).
- `input/gridX/`: static feature files (elevation, land use, anthropogenic heat,
  rainfall climatology, etc.) and the land mask.
- `input/xgb_models/`: XGBoost models (JSON) for air temperature, relative
  humidity, station pressure, and wind speed.

The static feature files and the XGBoost models are not included in this
repository because of their size (about 6.8 GB in total). They are available
from Zenodo:
https://doi.org/10.5281/zenodo.23877707

Download the archive and extract it into `input/` so that the directory
layout matches the list above. The paths in `input/input_features.txt` are
relative to `bin\`, so no editing is needed if this layout is kept.

The static feature files were derived from the following public datasets
(see [References](#references)); please also cite them when you use the
features:

| Feature files | Source dataset |
|---|---|
| `Elevation.*` | ETOPO1 1 Arc-Minute Global Relief Model (Amante and Eakins, 2009) |
| `LandUse.*` | MODIS land cover, MCD12C1 (Land_Cover_Type_1_Percent) |
| `AHE.*` | Global anthropogenic heat flux database (Dong et al., 2017) |
| `temperature_*`, `rainfall_*` | WFDE5 v2.1 (Cucchi et al., 2020) |

## Configuration

Copy an example configuration into `bin\` and edit it:

```bat
copy config\config_wfde5.example.ini  bin\config.ini
copy config\config_isimip.example.ini bin\config_isimip.ini
```

The configuration is an INI file (`;` starts a comment). The main sections are:

| Section | Keys | Description |
|---|---|---|
| `[System]` | `OutPutDir`, `AdjustDir` | output and bias-correction cache directories (must exist) |
| | `NumThreads`, `NcCompLevel`, `UtcOffsetHours` | OpenMP threads, deflate level, offset of log timestamps |
| | `Wfde5Dir`, `FileNameStr`, `FileVersion` | WFDE5 only: input location and file name parts |
| | `Isimip3bDir` | ISIMIP3b only: input root directory |
| `[CalcWBGT]` | `avg`, `urban`, `zspeed`, `dT` | parameters of `calc_wbgt()` |
| | `PGW` | pseudo-global-warming offset (°C) added to air temperature in the WBGT calculation; 0 = disabled |
| | `ThresholdTair` | currently unused (the key must be present) |
| `[Output]` | `SaveTair`, `SaveWbgtSun`, `SaveWbgtShade` | variables to compute and save (1/0; defaults 1, 1, 0) |
| | `QuantizeNSD` | significant digits kept by lossy quantization before compression; 0 = lossless (default) |
| `[Resolution]` | `longitude`, `latitude` | size of the global output grid (7200 × 3600 for 0.05°) |
| | `lon_min`, `lon_max`, `lat_min`, `lat_max` | computation domain (required) |
| `[FeatureFile]` | `FeatureList`, `landmask` | feature list file and land mask |
| `[XgbModel]` | `air_temperature`, `relative_humidity`, `station_pressure`, `wind_speed` | XGBoost model files |

See the comments in `config/*.example.ini` for details. The effective
configuration is printed at startup.

`QuantizeNSD = 3` reduces the output size by a factor of about 2.6 with a
maximum error of about 0.03 °C, but does not shorten the write time.

## Running

Pass the configuration as a file name relative to `bin\` (not as an absolute
path): the executable changes the working directory to its own location at
startup and opens `./{config}`.

```bat
cd bin
gclien-isa-wfde5.exe 2010 8 config.ini
gclien-isa-isimip.exe 2100 8 config_isimip.ini ssp585 GFDL-ESM4
```

The first run for a given month and domain computes the bias corrections and
stores them in `AdjustDir` as `adjust_{var}_{month}{domain_tag}.nc`. Later runs
for the same month and domain reuse them. Delete these files if you change the
models or the feature files.

## Output

One NetCDF file per variable is written to `OutPutDir`:

| File | Variable | Description |
|---|---|---|
| `tair_downscaled_{source}{YYYYMM}{domain_tag}.nc` | `tair` | air temperature (°C) |
| `wbgt_sun_downscaled_{source}{YYYYMM}{domain_tag}.nc` | `wbgt_sun` | WBGT in the sun (°C) |
| `wbgt_shade_downscaled_{source}{YYYYMM}{domain_tag}.nc` | `wbgt_shade` | WBGT in the shade, with solar radiation set to 0 (°C) |

- `{source}` is `wfde5_` for WFDE5 and `{gcm}_{scenario}_` for ISIMIP3b. With a
  non-zero `PGW`, `pgw{value}_` is inserted after `downscaled_`.
- `{domain_tag}` is `_lat_{start}_{end}_lon_{start}_{end}`, the index range of
  the domain on the global grid (end exclusive), e.g.
  `_lat_0_3600_lon_0_7200` for the whole globe.
- The variables have dimensions `time × latitude × longitude`. Non-land cells
  and cells without input data are filled with missing values.

The shaded WBGT excludes the radiant heat load of direct, diffuse, and
ground-reflected solar radiation, and therefore represents a lower bound
(complete shade).

## Repository structure

```
build.bat                  Windows build script
CMakeLists.txt             CMake project (two executables)
config/                    example configuration files
input/input_features.txt   list of static feature files
src/
  gclien_core.h/.cpp       shared processing (config, I/O, bias correction,
                           interpolation, WBGT loop, output)
  main_for_wfde5.cpp       main for WFDE5 input
  main_for_isimip3b.cpp    main for ISIMIP3b input
  wbgtn.h/.cpp             WBGT model (Liljegren et al., 2008; modified by NIES)
```

## License

This software is released under the BSD 3-Clause License. See
[LICENSE](LICENSE).

Copyright (c) 2026 Jun'ya Takakura, National Institute for Environmental
Studies (NIES).

`src/wbgtn.cpp` and `src/wbgtn.h` are derived from third-party code and are
also subject to the original licenses reproduced in [LICENSE](LICENSE): the
WBGT code by James C. Liljegren (Copyright (c) 2008, UChicago Argonne, LLC)
and its C version by Max Lieblich (MIT License).

## Acknowledgment

This product includes software produced by UChicago Argonne, LLC under
Contract No. DE-AC02-06CH11357 with the Department of Energy.

## References

- Amante, C., and Eakins, B. W. (2009). ETOPO1 1 Arc-Minute Global Relief
  Model: Procedures, Data Sources and Analysis. NOAA Technical Memorandum
  NESDIS NGDC-24. https://doi.org/10.7289/V5C8276M
- Bechtold, M., Poschlod, B., Otto, C., Volkholz, J., Büchner, M., and Zabel,
  F. (2026). A global hourly ISIMIP3 climate forcing dataset for impact
  modeling. *Earth System Science Data*, 18, 6637–6648.
  https://doi.org/10.5194/essd-18-6637-2026
- Copernicus Climate Change Service (C3S). Near surface meteorological
  variables from 1979 to 2019 derived from bias-corrected reanalysis
  (WFDE5), version 2.1. C3S Climate Data Store (CDS).
  https://doi.org/10.24381/cds.20d54e34
- Cucchi, M., Weedon, G. P., Amici, A., Bellouin, N., Lange, S., Müller
  Schmied, H., Hersbach, H., and Buontempo, C. (2020). WFDE5: bias-adjusted
  ERA5 reanalysis data for impact studies. *Earth System Science Data*, 12,
  2097–2120. https://doi.org/10.5194/essd-12-2097-2020
- Dong, Y., Varquez, A. C. G., and Kanda, M. (2017). Global anthropogenic
  heat flux database with high spatial resolution. *Atmospheric
  Environment*, 150, 276–294. https://doi.org/10.1016/j.atmosenv.2016.11.040
- Friedl, M., and Sulla-Menashe, D. (2022). MODIS/Terra+Aqua Land Cover Type
  Yearly L3 Global 0.05Deg CMG V061. NASA EOSDIS Land Processes DAAC.
  https://doi.org/10.5067/MODIS/MCD12C1.061
- Liljegren, J. C., Carhart, R. A., Lawday, P., Tschopp, S., and Sharp, R.
  (2008). Modeling the wet bulb globe temperature using standard
  meteorological measurements. *Journal of Occupational and Environmental
  Hygiene*, 5(10), 645–655. https://doi.org/10.1080/15459620802310770
