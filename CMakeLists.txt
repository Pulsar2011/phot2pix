# ──────────────────────────────────────────────────────────────────────────────
#  phot2pix — Photon-to-Pixel simulator library
#  Centre de Physique des Particules de Marseille (CPPM)
#
#  Typical usage:
#    cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<prefix>
#    cmake --build build -j$(nproc)
#    cmake --install build
#
#  Options:
#    -DCMAKE_BUILD_TYPE=Debug|Release|RelWithDebInfo
#    -DBUILD_TESTING=ON|OFF          (default ON)
#    -DDeepSkyROOT=<path>            (overrides search path for DST libs)
# ──────────────────────────────────────────────────────────────────────────────
cmake_minimum_required(VERSION 3.20)

project(phot2pix
    VERSION     0.1.0
    DESCRIPTION "Photon-to-Pixel CCD/CMOS sensor simulator library"
    LANGUAGES   CXX
)

# ─── Sanity: out-of-source build enforced ────────────────────────────────────
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds are not allowed.\n"
        "Please run:  cmake -B build  (or any directory other than the source root)."
    )
endif()

# ─── Global C++ standard ─────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD          23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

# Shared libs by default (produces .dylib / .so)
set(BUILD_SHARED_LIBS ON)

# ─── Helpers ─────────────────────────────────────────────────────────────────
include(GNUInstallDirs)         # CMAKE_INSTALL_{LIBDIR,INCLUDEDIR,BINDIR}
include(CheckIPOSupported)      # link-time optimisation probe
include(CMakePackageConfigHelpers)

# ─── Terminal colours (non-Windows only) ─────────────────────────────────────
if(NOT WIN32)
    string(ASCII 27 _Esc)
    set(_R  "${_Esc}[31m")
    set(_G  "${_Esc}[32m")
    set(_Y  "${_Esc}[33m")
    set(_B  "${_Esc}[34m")
    set(_Z  "${_Esc}[m")
endif()

# ─── Default build type ───────────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
        STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# ─── Default install prefix ───────────────────────────────────────────────────
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX
        "${CMAKE_SOURCE_DIR}/../install"
        CACHE PATH "Default install prefix" FORCE)
endif()

# ─── DeepSkyROOT: common search root for sibling DST libraries ───────────────
if(DEFINED ENV{DeepSkyROOT})
    set(DeepSkyROOT "$ENV{DeepSkyROOT}" CACHE PATH "Root of the DeepSky install tree")
endif()
if(NOT DeepSkyROOT)
    set(DeepSkyROOT "${CMAKE_INSTALL_PREFIX}")
endif()

message(STATUS "${_G}──── phot2pix configuration ────────────────────────${_Z}")
message(STATUS "${_G}  Build type    :${_Z} ${CMAKE_BUILD_TYPE}")
message(STATUS "${_G}  Install prefix:${_Z} ${CMAKE_INSTALL_PREFIX}")
message(STATUS "${_G}  DeepSkyROOT  :${_Z} ${DeepSkyROOT}")
message(STATUS "${_G}  Compiler      :${_Z} ${CMAKE_CXX_COMPILER}")
message(STATUS "${_G}  C++ standard  :${_Z} C++${CMAKE_CXX_STANDARD}")

# ─── Interface target: compile options shared by lib + tests ──────────────────
#  Using an INTERFACE target keeps all flags in one place and lets every
#  consumer inherit them via target_link_libraries(... PRIVATE _p2p_opts).
add_library(_p2p_opts INTERFACE)

target_compile_options(_p2p_opts INTERFACE
    -Wall
    -Wno-unknown-pragmas
    -fPIC
    -ffp-contract=on          # honour IEEE contracted FMA where available
)

# ── Build-type flags ──
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(_p2p_opts INTERFACE
        -O0 -g
        -fsanitize=address
        -fno-omit-frame-pointer
    )
    target_link_options(_p2p_opts INTERFACE -fsanitize=address)
else()
    target_compile_options(_p2p_opts INTERFACE -O3 -fomit-frame-pointer)
endif()

# ── Platform / architecture flags ──
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    target_compile_options(_p2p_opts INTERFACE -mcpu=native -mtune=native)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
        set(CMAKE_OSX_ARCHITECTURES "arm64")
        target_compile_definitions(_p2p_opts INTERFACE Darwin DarwinArm64)
    else()
        set(CMAKE_OSX_ARCHITECTURES "x86_64")
        target_compile_options(_p2p_opts INTERFACE -m64)
        target_compile_definitions(_p2p_opts INTERFACE Darwin DarwinX86_64)
    endif()
    # LTO on macOS: thin (debug builds skip it — LLDB can't follow bitcode OSOs)
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(_p2p_opts INTERFACE -flto=thin)
        target_link_options(  _p2p_opts INTERFACE -flto=thin)
    endif()

elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    cmake_policy(SET CMP0069 NEW)
    target_compile_options(_p2p_opts INTERFACE
        -pipe -pthread -fexceptions -Wmisleading-indentation
    )
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
        target_compile_options(_p2p_opts INTERFACE -m64 -mtune=native)
        target_compile_definitions(_p2p_opts INTERFACE Linuxx86_64)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
        target_compile_options(_p2p_opts INTERFACE -mtune=native)
        target_compile_definitions(_p2p_opts INTERFACE LinuxAarch64)
    endif()
    # Full LTO on Linux via CMake's IPO mechanism (respects linker capability)
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_err)
        if(_ipo_ok)
            set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
        else()
            message(STATUS "${_Y}IPO/LTO not supported: ${_ipo_err}${_Z}")
        endif()
    endif()
endif()

if(CI_CIWORKER)
    target_compile_definitions(_p2p_opts INTERFACE CI_CIWORKER)
endif()

# ─── Dependency: nlohmann/json (PixelConfig.h requires it) ───────────────────
find_package(nlohmann_json 3.10 QUIET)
if(nlohmann_json_FOUND)
    message(STATUS "${_G}  nlohmann/json :${_Z} found via find_package")
else()
    message(STATUS "${_Y}  nlohmann/json : not found — fetching via FetchContent${_Z}")
    include(FetchContent)
    FetchContent_Declare(nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG        v3.11.3
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(nlohmann_json)
endif()

# ─── Dependency: DSTmath (dynamic, required) ─────────────────────────────────
find_library(DSTmath_LIB
    NAMES DSTmath libDSTmath
    HINTS "${DeepSkyROOT}/${CMAKE_INSTALL_LIBDIR}"
          "${DeepSkyROOT}/lib"
    REQUIRED
)
# Derive include root from the library location: .../lib/libDSTmath.dylib → .../include
cmake_path(GET DSTmath_LIB PARENT_PATH _dstmath_libdir)
cmake_path(GET _dstmath_libdir PARENT_PATH _dstmath_root)
set(DSTmath_INCLUDE_DIR "${_dstmath_root}/include")

add_library(DST::math SHARED IMPORTED GLOBAL)
set_target_properties(DST::math PROPERTIES
    IMPORTED_LOCATION             "${DSTmath_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${DSTmath_INCLUDE_DIR}"
)
message(STATUS "${_G}  DSTmath lib   :${_Z} ${DSTmath_LIB}")
message(STATUS "${_G}  DSTmath inc   :${_Z} ${DSTmath_INCLUDE_DIR}")

# ─── Dependency: DSTfits (dynamic, required — used for full-matrix generation) ─
find_library(DSTfits_LIB
    NAMES DSTfits libDSTfits
    HINTS "${DeepSkyROOT}/${CMAKE_INSTALL_LIBDIR}"
          "${DeepSkyROOT}/lib"
    REQUIRED
)
cmake_path(GET DSTfits_LIB PARENT_PATH _dstfits_libdir)
cmake_path(GET _dstfits_libdir PARENT_PATH _dstfits_root)
set(DSTfits_INCLUDE_DIR "${_dstfits_root}/include")

add_library(DST::fits SHARED IMPORTED GLOBAL)
set_target_properties(DST::fits PROPERTIES
    IMPORTED_LOCATION             "${DSTfits_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${DSTfits_INCLUDE_DIR}"
)
message(STATUS "${_G}  DSTfits lib   :${_Z} ${DSTfits_LIB}")
message(STATUS "${_G}  DSTfits inc   :${_Z} ${DSTfits_INCLUDE_DIR}")

# ─── Library: phot2pix ────────────────────────────────────────────────────────
file(GLOB_RECURSE PHOT2PIX_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/phot2pix/*.cxx"
)
file(GLOB_RECURSE PHOT2PIX_HEADERS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/phot2pix/*.h"
)

add_library(phot2pix SHARED ${PHOT2PIX_SOURCES})

# PUBLIC include path: during build expose src/, at install time expose include/
target_include_directories(phot2pix
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(phot2pix
    PRIVATE _p2p_opts
    PUBLIC  nlohmann_json::nlohmann_json
    PUBLIC  DST::math
    PUBLIC  DST::fits
)

set_target_properties(phot2pix PROPERTIES
    OUTPUT_NAME  "phot2pix"
    VERSION      ${PROJECT_VERSION}
    SOVERSION    ${PROJECT_VERSION_MAJOR}
    # RPATH baked into the installed .dylib/.so so it finds the DST siblings
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
    BUILD_RPATH   "${CMAKE_CURRENT_BINARY_DIR}"
)

# ─── Tests ────────────────────────────────────────────────────────────────────
option(BUILD_TESTING "Build the GoogleTest suite" ON)

if(BUILD_TESTING)
    # googletest is a git submodule — build it as part of this project,
    # excluded from the default 'all' target and from installation.
    add_subdirectory(googletest EXCLUDE_FROM_ALL)
    enable_testing()
    include(GoogleTest)

    file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/test/*.cxx"
    )

    foreach(_test_src ${TEST_SOURCES})
        get_filename_component(_test_name ${_test_src} NAME_WE)

        add_executable(${_test_name} ${_test_src})

        target_include_directories(${_test_name}
            PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src"
        )
        target_link_libraries(${_test_name}
            PRIVATE phot2pix
            PRIVATE _p2p_opts
            PRIVATE GTest::gtest
            PRIVATE GTest::gtest_main
        )
        # Register tests with CTest (auto-discovers TEST() macros at build time)
        gtest_discover_tests(${_test_name})
    endforeach()
endif()

# ─── Install: library + headers ───────────────────────────────────────────────
install(TARGETS phot2pix
    EXPORT      phot2pixTargets
    LIBRARY     DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE     DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME     DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

install(FILES ${PHOT2PIX_HEADERS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/phot2pix"
)

# ─── Install: CMake package config (enables find_package(phot2pix)) ───────────
install(EXPORT phot2pixTargets
    FILE        phot2pixTargets.cmake
    NAMESPACE   phot2pix::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/phot2pix"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/phot2pixConfigVersion.cmake"
    VERSION      ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/phot2pixConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/phot2pixConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/phot2pix"
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/phot2pixConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/phot2pixConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/phot2pix"
)

# ─── distclean helper ─────────────────────────────────────────────────────────
add_custom_target(distclean
    COMMAND ${CMAKE_COMMAND} -E echo "Removing build artefacts..."
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}"
    COMMAND /usr/bin/env bash -c
        "[ -f '${CMAKE_BINARY_DIR}/install_manifest.txt' ] && xargs rm -f < '${CMAKE_BINARY_DIR}/install_manifest.txt' || true"
    COMMENT "distclean: removing build directory and installed files"
    VERBATIM
)

# ─── Summary ─────────────────────────────────────────────────────────────────
message(STATUS "${_G}────────────────────────────────────────────────────${_Z}")
message(STATUS "${_G}  Tests enabled :${_Z} ${BUILD_TESTING}")
message(STATUS "${_G}  lib output    :${_Z} ${CMAKE_CURRENT_BINARY_DIR}/libphot2pix${CMAKE_SHARED_LIBRARY_SUFFIX}")
message(STATUS "${_G}────────────────────────────────────────────────────${_Z}")
