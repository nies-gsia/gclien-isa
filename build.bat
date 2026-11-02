@echo off
rem Copyright (c) 2026 Jun'ya Takakura, National Institute for Environmental Studies (NIES)
rem SPDX-License-Identifier: BSD-3-Clause
setlocal

set "XGBOOST_DIR=D:\workspace\xgboost_cpu\xgboost"
set "NETCDF_DIR=C:\Program Files\netCDF 4.10.0"
set "ONEAPI_BIN=C:\Program Files (x86)\Intel\oneAPI\compiler\2026.0\bin"
set "ONEAPI_LIB=C:\Program Files (x86)\Intel\oneAPI\compiler\2026.0\lib"
set "CMAKE_BIN=C:\Program Files\CMake\bin"
set "NINJA_BIN=C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja"
set BUILD_DIR=build

rem Setup Intel compiler environment.
rem Use compiler vars.bat directly (more reliable than setvars.bat in subprocess contexts).
rem It also calls vcvarsall.bat to set up MSVC headers/libs required by icx.
call "C:\Program Files (x86)\Intel\oneAPI\compiler\2026.0\env\vars.bat" intel64 vs2022

rem Ensure cmake, ninja, and icx are in PATH regardless of whether vars.bat succeeded.
set "PATH=%CMAKE_BIN%;%NINJA_BIN%;%ONEAPI_BIN%;%PATH%"
set "LIB=%ONEAPI_LIB%;%LIB%"

if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"

echo.
echo === CMake configuration ===
echo.

cmake -S . -B "%BUILD_DIR%" -G "Ninja" ^
    "-DCMAKE_CXX_COMPILER=%ONEAPI_BIN%\icx.exe" ^
    -DCMAKE_BUILD_TYPE=Release ^
    "-DXGBOOST_DIR=%XGBOOST_DIR%" ^
    "-DNETCDF_DIR=%NETCDF_DIR%"

if %errorlevel% neq 0 (
    echo.
    echo === CMake configuration FAILED ===
    exit /b %errorlevel%
)

echo.
echo === Build ===
echo.

cmake --build "%BUILD_DIR%"

if %errorlevel% neq 0 (
    echo.
    echo === Build FAILED ===
    exit /b %errorlevel%
)

echo.
echo === Copying runtime DLLs to bin\ ===
copy /Y "%ONEAPI_BIN%\libiomp5md.dll"       bin\
copy /Y "%ONEAPI_BIN%\svml_dispmd.dll"      bin\
copy /Y "%ONEAPI_BIN%\libmmd.dll"           bin\
copy /Y "%XGBOOST_DIR%\lib\xgboost.dll"     bin\
copy /Y "%NETCDF_DIR%\bin\netcdf.dll"       bin\
copy /Y "%NETCDF_DIR%\bin\hdf5.dll"         bin\
copy /Y "%NETCDF_DIR%\bin\hdf5_hl.dll"      bin\
copy /Y "%NETCDF_DIR%\bin\zlib1.dll"        bin\

echo.
echo === Build succeeded: bin\gclien-isa-wfde5.exe / bin\gclien-isa-isimip.exe ===
