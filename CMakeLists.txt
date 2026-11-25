# ──────────────────────────────────────────────────────────────────────────────
#  phot2pix — Photon-to-Pixel simulator library
#  Centre de Physique des Particules de Marseille (CPPM)
#
#  Typical usage:
#    cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<prefix>
#    cmake --build build -j$(nproc)
#    cmake --install build
#    cmake --build build --target doc   (Doxygen HTML → build/doc/html)
#
#  Options:
#    -DCMAKE_BUILD_TYPE=Debug|Release|RelWithDebInfo
#    -DBUILD_TESTING=ON|OFF          (default ON)
#    -DDeepSkyROOT=<path>            (overrides search path for DST libs)
# ──────────────────────────────────────────────────────────────────────────────
cmake_minimum_required(VERSION 3.20)

# ─── Version: from exact git tag vX.Y.Z, or fall back to 0.1.0 ──────────────
find_package(Git QUIET)
set(_p2p_version "0.1.0")
if(Git_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --match "v*" --exact-match
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE _git_tag
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _git_rc
    )
    if(_git_rc EQUAL 0 AND _git_tag MATCHES "^v([0-9]+)\\.([0-9]+)\\.([0-9]+)")
        set(_p2p_version "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
    endif()
    # Full describe string for display (includes commit offset + hash when not on a tag).
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --match "v*" --always --dirty
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE PHOT2PIX_VERSION_STRING
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
endif()
if(NOT PHOT2PIX_VERSION_STRING)
    set(PHOT2PIX_VERSION_STRING "${_p2p_version}")
endif()

project(phot2pix
    VERSION     ${_p2p_version}
    DESCRIPTION "Photon-to-Pixel CCD/CMOS sensor simulator library"
    LANGUAGES   CXX
)

# ─── Generated header: version.h ─────────────────────────────────────────────
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/phot2pix/phot2pix.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/phot2pix/phot2pix.h"
    @ONLY
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
message(STATUS "${_G}  Version       :${_Z} ${PROJECT_VERSION}")
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
    message(STATUS "${_Y}  nlohmann/json : not found — fetching via FetchContent, installed with phot2pix${_Z}")
    include(FetchContent)
    # phot2pix links nlohmann_json PUBLIC (PixelConfig.h includes <nlohmann/json.hpp>), so
    # install(EXPORT phot2pixTargets) needs nlohmann_json in an export set. Its own install rules
    # (export set nlohmann_jsonTargets) are off when it is not the top-level project: switch them
    # on, the headers and CMake package are installed with phot2pix and phot2pixConfig.cmake
    # finds them again with find_dependency().
    set(JSON_Install ON CACHE INTERNAL "Install the fetched nlohmann/json with phot2pix")
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

# ─── Dependency: Minuit2 (shared, required by DSTfits) ───────────────────────
# Headers live at ${DeepSkyROOT}/include/Minuit2/Minuit2/*.h — so the include
# path must be ${DeepSkyROOT}/include/Minuit2 for <Minuit2/FCNBase.h> to resolve.
find_library(Minuit2_LIB
    NAMES Minuit2 libMinuit2
    HINTS "${DeepSkyROOT}/${CMAKE_INSTALL_LIBDIR}"
          "${DeepSkyROOT}/lib"
    REQUIRED
)
cmake_path(GET Minuit2_LIB PARENT_PATH _minuit2_libdir)
cmake_path(GET _minuit2_libdir PARENT_PATH _minuit2_root)
set(Minuit2_INCLUDE_DIR "${_minuit2_root}/include/Minuit2")

add_library(ROOT::Minuit2 SHARED IMPORTED GLOBAL)
set_target_properties(ROOT::Minuit2 PROPERTIES
    IMPORTED_LOCATION             "${Minuit2_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${Minuit2_INCLUDE_DIR}"
)
message(STATUS "${_G}  Minuit2 lib   :${_Z} ${Minuit2_LIB}")
message(STATUS "${_G}  Minuit2 inc   :${_Z} ${Minuit2_INCLUDE_DIR}")

# ─── Dependency: cfitsio (static, required by DSTfits) ───────────────────────
find_library(cfitsio_LIB
    NAMES cfitsio libcfitsio
    HINTS "${DeepSkyROOT}/${CMAKE_INSTALL_LIBDIR}"
          "${DeepSkyROOT}/lib"
    REQUIRED
)
add_library(DST::cfitsio STATIC IMPORTED GLOBAL)
set_target_properties(DST::cfitsio PROPERTIES
    IMPORTED_LOCATION             "${cfitsio_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${DeepSkyROOT}/include/cfitsio"
)
message(STATUS "${_G}  cfitsio lib   :${_Z} ${cfitsio_LIB}")

# ─── Dependency: wcslib (static, required by DSTfits) ────────────────────────
find_library(wcslib_LIB
    NAMES wcs libwcs
    HINTS "${DeepSkyROOT}/${CMAKE_INSTALL_LIBDIR}"
          "${DeepSkyROOT}/lib"
    REQUIRED
)
add_library(DST::wcs STATIC IMPORTED GLOBAL)
set_target_properties(DST::wcs PROPERTIES
    IMPORTED_LOCATION             "${wcslib_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${DeepSkyROOT}/include/wcslib"
)
message(STATUS "${_G}  wcslib lib    :${_Z} ${wcslib_LIB}")

# ─── Dependency: DSTfits (dynamic, required — used for full-matrix generation) ─
# Declares cfitsio, wcslib and Minuit2 as interface dependencies so their
# include paths and libraries propagate to all consumers of DST::fits.
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
    INTERFACE_LINK_LIBRARIES      "ROOT::Minuit2;DST::cfitsio;DST::wcs"
)
message(STATUS "${_G}  DSTfits lib   :${_Z} ${DSTfits_LIB}")
message(STATUS "${_G}  DSTfits inc   :${_Z} ${DSTfits_INCLUDE_DIR}")

# ─── Dependency: DSTools headers (required, headers only: libDSTools is not linked) ─
# DetectorEffect.h needs DSTmorpho.h, and the two headers it includes. When they are not
# installed, only these three files are downloaded from the deepskyutilities repository.
# The GitLab repository files API is used rather than the "raw" web link: for a private
# project, the latter answers with the HTML of the sign-in page instead of an error.
# Token, if the project is private: environment variable DEEPSKYFITS_TOKEN (read_repository).
set(DSTools_REF "main" CACHE STRING "Branch, tag or commit of deepskyutilities to download the DSTools headers from")

find_path(DSTools_INCLUDE_DIR
    NAMES DSTools/DSTmorpho.h
    HINTS "${DeepSkyROOT}/${CMAKE_INSTALL_INCLUDEDIR}"
          "${DeepSkyROOT}/include"
)
if(NOT DSTools_INCLUDE_DIR)
    set(_dstools_inc  "${CMAKE_CURRENT_BINARY_DIR}/_deps/DSTools/include")
    set(_dstools_url  "https://gitlab.in2p3.fr/api/v4/projects/deepskyproject%2Fdeepskyutilities/repository/files")
    set(_dstools_auth)
    if(NOT "$ENV{DEEPSKYFITS_TOKEN}" STREQUAL "")
        set(_dstools_auth HTTPHEADER "PRIVATE-TOKEN: $ENV{DEEPSKYFITS_TOKEN}")
    endif()

    message(STATUS "${_Y}  DSTools       : not found — downloading DSTmorpho.h (deepskyutilities, ${DSTools_REF})${_Z}")
    foreach(_h DSTmorpho.h DSTmorph_ops.h DST_version.h)
        file(DOWNLOAD
            "${_dstools_url}/include%2FDSTools%2F${_h}/raw?ref=${DSTools_REF}"
            "${_dstools_inc}/DSTools/${_h}"
            ${_dstools_auth}
            TLS_VERIFY ON
            TIMEOUT    60
            STATUS     _dstools_status
        )
        list(GET _dstools_status 0 _dstools_rc)
        if(NOT _dstools_rc EQUAL 0)
            list(GET _dstools_status 1 _dstools_msg)
            file(REMOVE_RECURSE "${_dstools_inc}")   # no partial file left for the next configure
            message(FATAL_ERROR "DSTools/${_h} is neither installed in ${DeepSkyROOT} nor downloadable: ${_dstools_msg}. "
                                "Install deepskyutilities, or set DEEPSKYFITS_TOKEN to a token that can read it.")
        endif()
    endforeach()

    set(DSTools_INCLUDE_DIR "${_dstools_inc}" CACHE PATH "Directory of DSTools/DSTmorpho.h" FORCE)
    set(DSTools_DOWNLOADED  TRUE CACHE INTERNAL "DSTools headers downloaded by phot2pix")
endif()
message(STATUS "${_G}  DSTools inc   :${_Z} ${DSTools_INCLUDE_DIR}")

# ─── Dependency: DSTCore (dynamic, optional — only needed by src/utilities) ──
# Provides DOption (OptionManager.h) and DSF::gGIT. If absent, the libraries
# still build and the utilities are skipped.
find_library(DSTCore_LIB
    NAMES DSTCore libDSTCore
    HINTS "${DeepSkyROOT}/${CMAKE_INSTALL_LIBDIR}"
          "${DeepSkyROOT}/lib"
)
if(DSTCore_LIB)
    cmake_path(GET DSTCore_LIB PARENT_PATH _dstcore_libdir)
    cmake_path(GET _dstcore_libdir PARENT_PATH _dstcore_root)
    set(DSTCore_INCLUDE_DIR "${_dstcore_root}/include")

    add_library(DST::core SHARED IMPORTED GLOBAL)
    set_target_properties(DST::core PROPERTIES
        IMPORTED_LOCATION             "${DSTCore_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${DSTCore_INCLUDE_DIR}"
    )
    message(STATUS "${_G}  DSTCore lib   :${_Z} ${DSTCore_LIB}")
else()
    message(STATUS "${_Y}  DSTCore lib   : not found — utilities will be skipped${_Z}")
endif()

# ─── Library: phot2pix ────────────────────────────────────────────────────────
file(GLOB_RECURSE PHOT2PIX_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/phot2pix/*.cxx"
)
file(GLOB_RECURSE PHOT2PIX_HEADERS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/phot2pix/*.h"
)

add_library(phot2pix SHARED ${PHOT2PIX_SOURCES})

# PUBLIC include path: during build expose src/ and generated/, at install time expose include/
target_include_directories(phot2pix
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
        $<BUILD_INTERFACE:${DSTools_INCLUDE_DIR}>
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

# ─── Library: nisppixelsimulator ─────────────────────────────────────────────
file(GLOB_RECURSE NISP_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/nisp/*.cxx"
)
file(GLOB_RECURSE NISP_HEADERS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/nisp/*.h"
)

if(NISP_SOURCES)
    add_library(nisppixelsimulator SHARED ${NISP_SOURCES})

    # Expose both src/ (for <nisp/...>) and src/phot2pix/ (for <PixelSimulator.h>
    # and <DetectorSimulator.h> used by the nisp headers directly).
    target_include_directories(nisppixelsimulator
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/phot2pix>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

    target_link_libraries(nisppixelsimulator
        PRIVATE _p2p_opts
        PUBLIC  phot2pix
    )

    set_target_properties(nisppixelsimulator PROPERTIES
        OUTPUT_NAME "nisppixelsimulator"
        VERSION     ${PROJECT_VERSION}
        SOVERSION   ${PROJECT_VERSION_MAJOR}
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
        BUILD_RPATH   "${CMAKE_CURRENT_BINARY_DIR}"
    )

    install(TARGETS nisppixelsimulator
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )
    install(FILES ${NISP_HEADERS}
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nisp"
    )

    message(STATUS "${_G}  nisp lib      :${_Z} nisppixelsimulator")
endif()

# ─── Utilities: nispSimulator ────────────────────────────────────────────────
if(TARGET nisppixelsimulator AND TARGET DST::core)
    add_executable(nispSimulator
        "${CMAKE_CURRENT_SOURCE_DIR}/src/utilities/nispSimulator.cxx"
    )

    target_link_libraries(nispSimulator
        PRIVATE nisppixelsimulator
        PRIVATE DST::core
        PRIVATE _p2p_opts
    )

    set_target_properties(nispSimulator PROPERTIES
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
    )

    install(TARGETS nispSimulator
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )

    message(STATUS "${_G}  utility       :${_Z} nispSimulator")
endif()

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
        # nisp tests additionally link against nisppixelsimulator
        if(TARGET nisppixelsimulator)
            target_link_libraries(${_test_name} PRIVATE nisppixelsimulator)
        endif()
        # Register tests with CTest (auto-discovers TEST() macros at build time)
        gtest_discover_tests(${_test_name})
    endforeach()
endif()

# ─── Documentation: Doxygen (target 'doc', not part of 'all') ────────────────
#  PH2PIX.Doxyfile remains the reference configuration. A small Doxyfile is
#  generated in the build tree: it @INCLUDEs the reference one and overrides
#  only what depends on this build (output dir, generated headers, version,
#  dot location). Doxygen runs from the source root, so the relative paths of
#  PH2PIX.Doxyfile (EXCLUDE, HTML_EXTRA_STYLESHEET, ...) keep working.
find_package(Doxygen OPTIONAL_COMPONENTS dot)

if(DOXYGEN_FOUND)
    if(DOXYGEN_DOT_FOUND)
        set(_p2p_doc_have_dot YES)
        cmake_path(GET DOXYGEN_DOT_EXECUTABLE PARENT_PATH _p2p_doc_dot_path)
    else()
        set(_p2p_doc_have_dot NO)
        set(_p2p_doc_dot_path "")
    endif()

    set(_p2p_doc_dir      "${CMAKE_CURRENT_BINARY_DIR}/doc")
    set(_p2p_doc_doxyfile "${CMAKE_CURRENT_BINARY_DIR}/PH2PIX.Doxyfile")

    # file(CONFIGURE) only rewrites the file when its content changes.
    file(CONFIGURE OUTPUT "${_p2p_doc_doxyfile}" @ONLY CONTENT [=[
@INCLUDE               = "@CMAKE_CURRENT_SOURCE_DIR@/PH2PIX.Doxyfile"
PROJECT_NUMBER         = "@PHOT2PIX_VERSION_STRING@"
INPUT                  = "@CMAKE_CURRENT_SOURCE_DIR@/src" \
                         "@CMAKE_CURRENT_BINARY_DIR@/generated"
OUTPUT_DIRECTORY       = "@_p2p_doc_dir@"
HAVE_DOT               = @_p2p_doc_have_dot@
DOT_PATH               = "@_p2p_doc_dot_path@"
]=])

    add_custom_target(doc
        COMMAND Doxygen::doxygen "${_p2p_doc_doxyfile}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMENT "Generating Doxygen documentation in ${_p2p_doc_dir}/html"
        VERBATIM
    )

    message(STATUS "${_G}  Doxygen       :${_Z} ${DOXYGEN_EXECUTABLE} (dot: ${_p2p_doc_have_dot}) — target 'doc'")
else()
    message(STATUS "${_Y}  Doxygen       : not found — 'doc' target unavailable${_Z}")
endif()

# ─── Install: library + headers ───────────────────────────────────────────────
# Remove stale libphot2pix* files left by any previously installed version.
install(CODE "
    file(GLOB _old_libs \"\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libphot2pix*\")
    foreach(_lib \${_old_libs})
        message(STATUS \"Removing stale: \${_lib}\")
        file(REMOVE \"\${_lib}\")
    endforeach()
")

install(TARGETS phot2pix
    EXPORT      phot2pixTargets
    LIBRARY     DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE     DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME     DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

install(FILES ${PHOT2PIX_HEADERS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/phot2pix"
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/generated/phot2pix/phot2pix.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/phot2pix"
)
# The installed DetectorEffect.h includes <DSTools/DSTmorpho.h>: provide it when it was downloaded.
if(DSTools_DOWNLOADED)
    install(DIRECTORY "${DSTools_INCLUDE_DIR}/DSTools" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endif()

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
