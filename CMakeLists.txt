# Copyright (c) 2026 Jun'ya Takakura, National Institute for Environmental Studies (NIES)
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required(VERSION 3.20)
project(gclien-isa CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Output to bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY         ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_SOURCE_DIR}/bin)

# ------------------------------------------------------------
# OpenMP
# ------------------------------------------------------------
find_package(OpenMP REQUIRED)

# ------------------------------------------------------------
# XGBoost  (no CMake config is installed, so locate it directly with find_*)
# ------------------------------------------------------------
set(XGBOOST_DIR "" CACHE PATH "XGBoost root directory")

find_path(XGBOOST_INCLUDE_DIR
    NAMES xgboost/c_api.h
    HINTS "${XGBOOST_DIR}/include"
    NO_DEFAULT_PATH
)
find_library(XGBOOST_LIBRARY
    NAMES xgboost
    HINTS "${XGBOOST_DIR}/lib"
    NO_DEFAULT_PATH
)

if(NOT XGBOOST_INCLUDE_DIR)
    message(FATAL_ERROR "xgboost/c_api.h not found. Set XGBOOST_DIR to the XGBoost root directory.")
endif()
if(NOT XGBOOST_LIBRARY)
    message(FATAL_ERROR "xgboost library not found. Set XGBOOST_DIR to the XGBoost root directory.")
endif()

message(STATUS "XGBoost include : ${XGBOOST_INCLUDE_DIR}")
message(STATUS "XGBoost library : ${XGBOOST_LIBRARY}")

# ------------------------------------------------------------
# NetCDF  (installed Config has hardcoded paths, use find_* directly)
# ------------------------------------------------------------
set(NETCDF_DIR "" CACHE PATH "NetCDF root directory")

find_path(NETCDF_INCLUDE_DIR
    NAMES netcdf.h
    HINTS "${NETCDF_DIR}/include"
    NO_DEFAULT_PATH
)
find_library(NETCDF_LIBRARY
    NAMES netcdf
    HINTS "${NETCDF_DIR}/lib"
    NO_DEFAULT_PATH
)

if(NOT NETCDF_INCLUDE_DIR)
    message(FATAL_ERROR "netcdf.h not found. Set NETCDF_DIR to the NetCDF installation directory.")
endif()
if(NOT NETCDF_LIBRARY)
    message(FATAL_ERROR "netcdf library not found. Set NETCDF_DIR to the NetCDF installation directory.")
endif()

message(STATUS "NetCDF include : ${NETCDF_INCLUDE_DIR}")
message(STATUS "NetCDF library : ${NETCDF_LIBRARY}")

# ------------------------------------------------------------
# Targets
#   gclien-isa-wfde5  : for WFDE5 input
#   gclien-isa-isimip : for ISIMIP3b input
#   Shared processing is consolidated in src/gclien_core.cpp
# ------------------------------------------------------------
add_executable(gclien-isa-wfde5
    src/wbgtn.cpp
    src/gclien_core.cpp
    src/main_for_wfde5.cpp
)

add_executable(gclien-isa-isimip
    src/wbgtn.cpp
    src/gclien_core.cpp
    src/main_for_isimip3b.cpp
)

foreach(tgt gclien-isa-wfde5 gclien-isa-isimip)
    # Exceptions are not used at all (errors are handled by immediate exit via handle_error() or by return values).
    # icx (MSVC-compatible driver) does not enable exceptions by default, so /EHsc is not needed.
    # Define _HAS_EXCEPTIONS=0 to reliably disable exceptions and to notice at build time if throw/try
    # is introduced by mistake (exceptions inside the standard library terminate the process).
    target_compile_definitions(${tgt} PRIVATE _HAS_EXCEPTIONS=0)
    target_include_directories(${tgt} PRIVATE
        ${XGBOOST_INCLUDE_DIR}
        ${NETCDF_INCLUDE_DIR}
    )
    target_link_libraries(${tgt} PRIVATE
        OpenMP::OpenMP_CXX
        ${XGBOOST_LIBRARY}
        ${NETCDF_LIBRARY}
    )
endforeach()

# ------------------------------------------------------------
# Windows-specific: add icons etc. below if needed
# ------------------------------------------------------------
if(WIN32)
    # Add Windows resources here in the future if needed
endif()
