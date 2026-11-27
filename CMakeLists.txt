#-------------------- CMAKE Configuration
# 3.24: the minimum of the other DeepSky projects (the presets need 3.21; the SYSTEM property set on the
# imported DeepSky targets needs 3.25, older versions ignore it).
cmake_minimum_required(VERSION 3.24)

# -Wno-dev by default (developer and deprecation warnings, mostly from the submodules): no need to
# pass it, nor to alias cmake with it, which breaks 'cmake --build' and 'cmake --install'.
# 'cmake -Wdev ...' shows them again.
if(NOT DEFINED CACHE{CMAKE_SUPPRESS_DEVELOPER_WARNINGS})
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON CACHE INTERNAL "Suppress developer warnings (-Wno-dev)")
endif()

# Everything is generated in a build directory, never in the sources: from the project directory,
# 'cmake --preset default' (or 'cmake -B build'); from a build directory, 'cmake <path to project>'.
get_filename_component(_src "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(_bin "${CMAKE_BINARY_DIR}" REALPATH)
if(_src STREQUAL _bin)
    message(FATAL_ERROR
        "In-source builds are not allowed. Remove what this run (and any former in-source build) created:\n"
        "    cmake -P cmake/CleanInSourceBuild.cmake\n"
        "then, from the project directory:\n"
        "    cmake --preset default        (or: cmake -B build)\n"
        "    cmake --build build -j\n"
        "    cmake --build build --target install\n"
        "Other presets: debug (build-debug/), asan (build-asan/).\n"
        "Or from a build directory: mkdir build && cd build && cmake .. && make")
endif()

project(phot2pix
    DESCRIPTION "Photon-to-Pixel CCD/CMOS sensor simulator library"
    LANGUAGES   CXX
)

if(NOT WIN32)
    string(ASCII 27 Esc)
    set(ColourReset "${Esc}[m")
    set(ColourBold  "${Esc}[1m")
    set(Red         "${Esc}[31m")
    set(Green       "${Esc}[32m")
    set(Yellow      "${Esc}[33m")
    set(Blue        "${Esc}[34m")
    set(Magenta     "${Esc}[35m")
    set(Cyan        "${Esc}[36m")
    set(White       "${Esc}[37m")
    set(BoldRed     "${Esc}[1;31m")
    set(BoldGreen   "${Esc}[1;32m")
    set(BoldYellow  "${Esc}[1;33m")
    set(BoldBlue    "${Esc}[1;34m")
    set(BoldMagenta "${Esc}[1;35m")
    set(BoldCyan    "${Esc}[1;36m")
    set(BoldWhite   "${Esc}[1;37m")
endif()

message("${BoldCyan}----------- phot2pix -----------${ColourReset}")
message("${Green}-- Project source directory :${ColourReset} ${CMAKE_SOURCE_DIR}")
message("${Green}-- Project build directory  :${ColourReset} ${CMAKE_BINARY_DIR}")

if(DEFINED ENV{DeepSkyROOT})
    set(DeepSkyROOT $ENV{DeepSkyROOT})
else()
    set(DeepSkyROOT ${CMAKE_SOURCE_DIR}/..)
endif()
get_filename_component(DeepSkyROOT "${DeepSkyROOT}" ABSOLUTE)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${DeepSkyROOT} CACHE PATH "default install path" FORCE)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

#-------------------- Options
option(WITH_ASAN   "Build phot2pix and its tests with AddressSanitizer" OFF)
option(BUILD_TESTS "Build the GoogleTest suite" ON)
option(CI_CIWORKER "Build for the GitLab CI runner" OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CI_CIWORKER)
    add_compile_definitions(CI_CIWORKER)
endif()

# Products of this build: executables in bin/, libraries in lib/ (objects stay in CMakeFiles/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# phot2pix is a C++23 library: no fallback to an older standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# compile_commands.json in the build directory, for clangd / IntelliSense
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#-------------------- COMPILER
# These flags apply to every C++ target of this directory, the googletest submodule included.
set(_release "$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>")
set(_debug   "$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>")

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(OperatingSystem "Linux")
    add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-Wall;-pipe;-pthread;-fexceptions;-Wno-unknown-pragmas;-Wmisleading-indentation;-mtune=native;-ffp-contract=on>")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
        add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-m64>")
        add_compile_definitions(Linuxx86_64)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
        add_compile_definitions(Linuxaarch64)
    endif()
    add_compile_options("$<${_release}:-flto>")
    add_link_options(-pthread "$<$<CONFIG:Release>:-flto>")
    # Interprocedural optimization (LTO) of the release build, when supported
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
    if(ipo_supported)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    endif()
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(OperatingSystem "Mac OS X")
    add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wno-unknown-pragmas;-mcpu=native;-mtune=native;-ffp-contract=on>")
    add_compile_definitions(Darwin)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
        add_compile_definitions(Darwinarm64)
        set(CMAKE_OSX_ARCHITECTURES "arm64")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
        add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-m64>")
        add_compile_definitions(Darwinx86_64)
        set(CMAKE_OSX_ARCHITECTURES "x86_64")
    endif()
    add_compile_options("$<${_release}:-flto=thin>")
    add_link_options("$<$<CONFIG:Release>:-flto=thin>")
endif()

add_compile_options("$<${_release}:-O3;-fomit-frame-pointer>")
# No LTO in debug: on macOS it replaces the .o files with LLVM bitcode, leaving the debug map with
# empty OSO entries that LLDB cannot follow.
add_compile_options("$<${_debug}:-O0;-g>")
add_compile_definitions("$<$<CONFIG:Debug>:__DEBUG__>")

if(WITH_ASAN)
    add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fsanitize=address;-fno-omit-frame-pointer;-g>")
    add_link_options(-fsanitize=address)
endif()

message("${BoldYellow}-- CMAKE INSTALL PREFIX : ${ColourReset}${CMAKE_INSTALL_PREFIX}")
message("${Green}-- DeepSkyROOT :${ColourReset} ${DeepSkyROOT}")
message("${Green}-- BUILD TYPE :${ColourReset} ${CMAKE_BUILD_TYPE}")
message("${Green}-- CXX COMPILER :${ColourReset} ${CMAKE_CXX_COMPILER}")
message("${Green}-- USE CXX STANDARD :${ColourReset} ${CMAKE_CXX_STANDARD}")
if(WITH_ASAN)
    message("${Blue}-- AddressSanitizer ACTIVATED${ColourReset}")
endif()

#-------------------- VERSION FROM THE GIT TAG
find_package(Git)
include(GetGitRevisionDescription)
git_describe(VERSION --always --tags --long)

string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")

# A clone without a reachable tag (e.g. a shallow CI clone) has no numeric version: the libraries
# stay unversioned and no CMake package is installed.
if(VERSION MATCHES "^v[0-9]+\\.[0-9]+\\.[0-9]+")
    set(phot2pix_NUMERIC_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
    set(PROJECT_VERSION "${phot2pix_NUMERIC_VERSION}")
else()
    message(WARNING "No version tag found (git describe gave '${VERSION}'): phot2pix will not be versioned")
    # phot2pix.h still needs numbers for PHOT2PIX_VERSION_MAJOR/MINOR/PATCH
    set(VERSION_MAJOR 0)
    set(VERSION_MINOR 0)
    set(VERSION_PATCH 0)
endif()
set(VERSION_SHORT "v${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
message("-- PKG VERSION : ${VERSION} (${VERSION_SHORT})")

# Read by phot2pix.h.in
set(PHOT2PIX_VERSION_STRING "${VERSION}")
set(PROJECT_VERSION_MAJOR   "${VERSION_MAJOR}")
set(PROJECT_VERSION_MINOR   "${VERSION_MINOR}")
set(PROJECT_VERSION_PATCH   "${VERSION_PATCH}")

# Generated sources stay in the build directory
set(GENERATED_DIR ${CMAKE_BINARY_DIR}/generated)
configure_file(${CMAKE_SOURCE_DIR}/src/phot2pix/phot2pix.h.in ${GENERATED_DIR}/phot2pix/phot2pix.h @ONLY)

#-------------------- DeepSky libraries (installed in DeepSkyROOT)
message("${Red}----------- SEARCH DeepSky libraries -----------${ColourReset}")
# Minimum versions of the DeepSky libraries: the latest tag of each package.
set(DSTmath_MIN_VERSION 0.2.0)
set(DSTfits_MIN_VERSION 1.5.0)
set(DSTCore_MIN_VERSION 2.0.1)   # optional: only nispSimulator needs it
set(DSTools_MIN_VERSION 2.2.0)   # headers only; downloaded when not installed (see below)

# Only look in DeepSkyROOT: the RPATH points there, and CMake would otherwise fall back on another
# install (e.g. found through PATH) when the one in DeepSkyROOT is too old.
foreach(_dst_pkg DSTmath DSTfits)
    # A <pkg>_DIR cached by an earlier configure is used as is, whatever PATHS and NO_DEFAULT_PATH say:
    # forget it, so that a copy installed elsewhere (e.g. in /usr) cannot come back from the cache.
    unset(${_dst_pkg}_DIR CACHE)
    find_package(${_dst_pkg} ${${_dst_pkg}_MIN_VERSION} REQUIRED PATHS ${DeepSkyROOT} NO_DEFAULT_PATH)
    # Imported targets are system targets by default, which turns ${DeepSkyROOT}/include into -isystem:
    # it then loses its -I priority over other installed copies and hides warnings from our own headers.
    set_target_properties(${_dst_pkg}::${_dst_pkg} PROPERTIES SYSTEM OFF)
    message("-- ${Green}Found${ColourReset} ${_dst_pkg} ${${_dst_pkg}_VERSION}: ${${_dst_pkg}_DIR}")
endforeach()

# Optional: without DSTCore (DOption, DSF::gGIT), the libraries are built but not nispSimulator
unset(DSTCore_DIR CACHE)
find_package(DSTCore ${DSTCore_MIN_VERSION} QUIET PATHS ${DeepSkyROOT} NO_DEFAULT_PATH)
if(DSTCore_FOUND)
    set_target_properties(DSTCore::DSTCore PROPERTIES SYSTEM OFF)
    message("-- ${Green}Found${ColourReset} DSTCore ${DSTCore_VERSION}: ${DSTCore_DIR}")
else()
    message("${Yellow}-- DSTCore >= ${DSTCore_MIN_VERSION} not found (versions considered: ${DSTCore_CONSIDERED_VERSIONS}): nispSimulator will not be built${ColourReset}")
endif()

# The DSTfits headers include cfitsio, wcslib and Minuit2 headers, installed beside them by deepskyfits
# with their libraries. Imported target, as the DeepSky ones: phot2pixConfig.cmake defines it again for
# the projects that use the installed phot2pix.
find_library(MINUIT2_LIBRARY Minuit2 PATHS ${DeepSkyROOT}/lib NO_DEFAULT_PATH REQUIRED)
find_library(CFITSIO_LIBRARY cfitsio PATHS ${DeepSkyROOT}/lib NO_DEFAULT_PATH REQUIRED)
find_library(WCSLIB_LIBRARY  wcs     PATHS ${DeepSkyROOT}/lib NO_DEFAULT_PATH REQUIRED)
foreach(_dir cfitsio wcslib Minuit2)
    if(NOT EXISTS ${DeepSkyROOT}/include/${_dir})
        message(FATAL_ERROR "${DeepSkyROOT}/include/${_dir} not found: install deepskyfits first")
    endif()
endforeach()
add_library(DeepSky::externals INTERFACE IMPORTED)
set_target_properties(DeepSky::externals PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${DeepSkyROOT}/include/cfitsio;${DeepSkyROOT}/include/wcslib;${DeepSkyROOT}/include/Minuit2"
    INTERFACE_LINK_LIBRARIES      "${MINUIT2_LIBRARY};${CFITSIO_LIBRARY};${WCSLIB_LIBRARY}"
)
message("-- ${Green}Found${ColourReset} Minuit2: ${MINUIT2_LIBRARY}")
message("-- ${Green}Found${ColourReset} cfitsio: ${CFITSIO_LIBRARY}")
message("-- ${Green}Found${ColourReset} wcslib : ${WCSLIB_LIBRARY}")

#-------------------- DSTools headers (required; headers only: libDSTools is not linked)
# DetectorEffect.h includes DSTools/DSTmorpho.h, which includes DSTmorph_ops.h and DST_version.h.
# - DSTools >= DSTools_MIN_VERSION installed in DeepSkyROOT: its headers are there, in the include
#   directory of DSTmath and DSTfits.
# - No DSTools package: only these three files are downloaded from deepskyutilities (DSTools_REF), rather
#   than building the whole package and its DSTdisplay and Qt6 dependencies, and installed with phot2pix.
# - An older DSTools package: error, since installing newer headers would replace those of its library.
# The GitLab repository files API is used rather than the "raw" web link: for a private project, the
# latter answers with the HTML of the sign-in page instead of an error.
# Token, if the project is private: environment variable DEEPSKYFITS_TOKEN (read_repository).
set(DSTools_REF "main" CACHE STRING "Branch, tag or commit of deepskyutilities to download the DSTools headers from")
set(DSTools_DOWNLOADED FALSE)

unset(DSTools_DIR CACHE)
find_package(DSTools ${DSTools_MIN_VERSION} QUIET PATHS ${DeepSkyROOT} NO_DEFAULT_PATH)
if(DSTools_FOUND)
    message("-- ${Green}Found${ColourReset} DSTools ${DSTools_VERSION} (headers only): ${DSTools_DIR}")
elseif(EXISTS ${DeepSkyROOT}/lib/cmake/DSTools/DSToolsConfig.cmake)
    message(FATAL_ERROR "${DeepSkyROOT} holds a DSTools older than ${DSTools_MIN_VERSION}, or whose package could not "
                        "be loaded (versions considered: ${DSTools_CONSIDERED_VERSIONS}): update deepskyutilities.")
else()
    # Downloaded once per reference: delete this directory to download again
    string(MAKE_C_IDENTIFIER "${DSTools_REF}" _dstools_ref)
    set(DSTools_INCLUDE_DIR "${CMAKE_BINARY_DIR}/_deps/DSTools-${_dstools_ref}/include")
    set(_dstools_url "https://gitlab.in2p3.fr/api/v4/projects/deepskyproject%2Fdeepskyutilities/repository/files")
    set(_dstools_auth)
    if(NOT "$ENV{DEEPSKYFITS_TOKEN}" STREQUAL "")
        set(_dstools_auth HTTPHEADER "PRIVATE-TOKEN: $ENV{DEEPSKYFITS_TOKEN}")
    endif()

    foreach(_h DSTmorpho.h DSTmorph_ops.h DST_version.h)
        if(EXISTS "${DSTools_INCLUDE_DIR}/DSTools/${_h}")
            continue()
        endif()
        message("${Yellow}-- DSTools not found: downloading DSTools/${_h} (deepskyutilities, ${DSTools_REF})${ColourReset}")
        file(DOWNLOAD
            "${_dstools_url}/include%2FDSTools%2F${_h}/raw?ref=${DSTools_REF}"
            "${DSTools_INCLUDE_DIR}/DSTools/${_h}"
            ${_dstools_auth}
            TLS_VERIFY ON
            TIMEOUT    60
            STATUS     _dstools_status
        )
        list(GET _dstools_status 0 _dstools_rc)
        if(NOT _dstools_rc EQUAL 0)
            list(GET _dstools_status 1 _dstools_msg)
            file(REMOVE_RECURSE "${DSTools_INCLUDE_DIR}")   # no partial file left for the next configure
            message(FATAL_ERROR "DSTools/${_h} is neither installed in ${DeepSkyROOT} nor downloadable: ${_dstools_msg}. "
                                "Install deepskyutilities, or set DEEPSKYFITS_TOKEN to a token that can read it.")
        endif()
    endforeach()
    set(DSTools_DOWNLOADED TRUE)
    message("-- ${Green}DSTools headers${ColourReset} downloaded (${DSTools_REF}): ${DSTools_INCLUDE_DIR}")
endif()

#-------------------- nlohmann/json (PixelConfig.h includes it)
find_package(nlohmann_json 3.10 QUIET)
if(nlohmann_json_FOUND)
    message("-- ${Green}Found${ColourReset} nlohmann/json ${nlohmann_json_VERSION}: ${nlohmann_json_DIR}")
else()
    message("${Yellow}-- nlohmann/json not found: fetched with FetchContent, installed with phot2pix${ColourReset}")
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

#-------------------- RPATH
# This build comes first, then DeepSkyROOT: an installed copy must not replace what the tests are testing.
set(CMAKE_BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY};${DeepSkyROOT}/lib")
# Installed: <prefix>/exe and <prefix>/lib find <prefix>/lib, then DeepSkyROOT
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path/../lib;${DeepSkyROOT}/lib")
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib;${DeepSkyROOT}/lib")
endif()
# Also the directories of the external libraries
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    # RPATH rather than RUNPATH: also used to find the dependencies of the DeepSky libraries
    add_link_options(-Wl,--disable-new-dtags)
endif()

#-------------------- libphot2pix
message("${Red}----------- libphot2pix -----------${ColourReset}")
file(GLOB_RECURSE PHOT2PIX_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/phot2pix/*.cxx)
file(GLOB_RECURSE PHOT2PIX_HEADERS CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/phot2pix/*.h)

add_library(phot2pix SHARED ${PHOT2PIX_SOURCES} ${PHOT2PIX_HEADERS})
if(phot2pix_NUMERIC_VERSION)
    set_target_properties(phot2pix PROPERTIES
        VERSION   "${phot2pix_NUMERIC_VERSION}"
        SOVERSION "${VERSION_MAJOR}"
    )
endif()
# Built: src/ and the generated phot2pix.h; installed: include/
target_include_directories(phot2pix PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${GENERATED_DIR}>
    $<INSTALL_INTERFACE:include>
)
if(DSTools_DOWNLOADED)
    target_include_directories(phot2pix PUBLIC $<BUILD_INTERFACE:${DSTools_INCLUDE_DIR}>)
endif()
target_link_libraries(phot2pix PUBLIC
    nlohmann_json::nlohmann_json DSTfits::DSTfits DSTmath::DSTmath DeepSky::externals
)

#-------------------- libnisppixelsimulator (Euclid NISP detectors, built on phot2pix)
file(GLOB_RECURSE NISP_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/nisp/*.cxx)
file(GLOB_RECURSE NISP_HEADERS CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/nisp/*.h)

if(NISP_SOURCES)
    message("${Red}----------- libnisppixelsimulator -----------${ColourReset}")
    add_library(nisppixelsimulator SHARED ${NISP_SOURCES} ${NISP_HEADERS})
    if(phot2pix_NUMERIC_VERSION)
        set_target_properties(nisppixelsimulator PROPERTIES
            VERSION   "${phot2pix_NUMERIC_VERSION}"
            SOVERSION "${VERSION_MAJOR}"
        )
    endif()
    # src/ for <nisp/...>, src/phot2pix/ for <PixelSimulator.h> and <DetectorSimulator.h>, which the
    # nisp headers include directly
    target_include_directories(nisppixelsimulator PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/phot2pix>
        $<INSTALL_INTERFACE:include>
    )
    target_link_libraries(nisppixelsimulator PUBLIC phot2pix)
endif()

#-------------------- nispSimulator (needs DSTCore)
if(TARGET nisppixelsimulator AND DSTCore_FOUND)
    add_executable(nispSimulator ${CMAKE_SOURCE_DIR}/src/utilities/nispSimulator.cxx)
    target_link_libraries(nispSimulator PRIVATE nisppixelsimulator DSTCore::DSTCore)
endif()

#-------------------- Tests
if(BUILD_TESTS)
    message("${Red}----------- GoogleTest -----------${ColourReset}")
    if(NOT EXISTS ${CMAKE_SOURCE_DIR}/googletest/CMakeLists.txt)
        message(FATAL_ERROR "${BoldRed}googletest is missing: run 'git submodule update --init' (or configure with -DBUILD_TESTS=OFF)${ColourReset}")
    endif()
    set(PHOT2PIX_TESTS ON)
    enable_testing()
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(BUILD_GMOCK   OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF)
    add_subdirectory(googletest EXCLUDE_FROM_ALL)
    include(GoogleTest)

    # The tests run from the build directory: whatever they write (test_data/) stays there
    file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/test/*.cxx)
    set(TEST_NAMES "")
    foreach(test_source IN LISTS TEST_SOURCES)
        get_filename_component(tname ${test_source} NAME_WE)
        list(APPEND TEST_NAMES ${tname})
        add_executable(${tname} ${test_source})
        target_link_libraries(${tname} PRIVATE phot2pix gtest gtest_main)
        # The nisp tests also need libnisppixelsimulator
        if(TARGET nisppixelsimulator)
            target_link_libraries(${tname} PRIVATE nisppixelsimulator)
        endif()
        gtest_discover_tests(${tname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    endforeach()
    message("${Blue}-- Test executables ${ColourReset} ${TEST_NAMES} (run with ctest)")
else()
    set(PHOT2PIX_TESTS OFF)
endif()

#-------------------- Doxygen documentation (target: doc, written in build/doc)
# PH2PIX.Doxyfile remains the reference configuration. The Doxyfile generated in the build directory
# @INCLUDEs it and overrides only what depends on this build (version, generated headers, output
# directory, dot). Doxygen runs from the sources, so the relative paths of PH2PIX.Doxyfile
# (EXCLUDE, HTML_EXTRA_STYLESHEET, ...) keep working.
find_package(Doxygen OPTIONAL_COMPONENTS dot)
if(DOXYGEN_FOUND)
    if(DOXYGEN_DOT_FOUND)
        set(_p2p_doc_have_dot YES)
        cmake_path(GET DOXYGEN_DOT_EXECUTABLE PARENT_PATH _p2p_doc_dot_path)
    else()
        set(_p2p_doc_have_dot NO)
        set(_p2p_doc_dot_path "")
    endif()
    set(_p2p_doc_dir      "${CMAKE_BINARY_DIR}/doc")
    set(_p2p_doc_doxyfile "${GENERATED_DIR}/PH2PIX.Doxyfile")

    # file(CONFIGURE) only rewrites the file when its content changes
    file(CONFIGURE OUTPUT "${_p2p_doc_doxyfile}" @ONLY CONTENT [=[
@INCLUDE               = "@CMAKE_SOURCE_DIR@/PH2PIX.Doxyfile"
PROJECT_NUMBER         = "@PHOT2PIX_VERSION_STRING@"
INPUT                  = "@CMAKE_SOURCE_DIR@/src" \
                         "@GENERATED_DIR@"
OUTPUT_DIRECTORY       = "@_p2p_doc_dir@"
HAVE_DOT               = @_p2p_doc_have_dot@
DOT_PATH               = "@_p2p_doc_dot_path@"
]=])
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/PH2PIX.Doxyfile)

    add_custom_target(doc
        COMMAND Doxygen::doxygen "${_p2p_doc_doxyfile}"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen in ${_p2p_doc_dir}/html"
        VERBATIM
    )
endif()

#-------------------- Installation
install(FILES ${PHOT2PIX_HEADERS} ${GENERATED_DIR}/phot2pix/phot2pix.h DESTINATION include/phot2pix)
# The installed DetectorEffect.h includes <DSTools/DSTmorpho.h>: provide it when it was downloaded
if(DSTools_DOWNLOADED)
    install(DIRECTORY ${DSTools_INCLUDE_DIR}/DSTools DESTINATION include)
endif()
install(TARGETS phot2pix
    EXPORT  phot2pixTargets
    RUNTIME DESTINATION exe
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
if(TARGET nisppixelsimulator)
    install(FILES ${NISP_HEADERS} DESTINATION include/nisp)
    install(TARGETS nisppixelsimulator
        RUNTIME DESTINATION exe
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
endif()
if(TARGET nispSimulator)
    install(TARGETS nispSimulator RUNTIME DESTINATION exe)
endif()

if(phot2pix_NUMERIC_VERSION)
    # CMake package, so that dependent projects can call find_package(phot2pix <version>)
    include(CMakePackageConfigHelpers)
    set(phot2pix_CMAKE_DIR lib/cmake/phot2pix)
    install(EXPORT phot2pixTargets
        FILE        phot2pixTargets.cmake
        NAMESPACE   phot2pix::
        DESTINATION ${phot2pix_CMAKE_DIR}
    )
    configure_package_config_file(
        ${CMAKE_SOURCE_DIR}/cmake/phot2pixConfig.cmake.in
        ${GENERATED_DIR}/phot2pixConfig.cmake
        INSTALL_DESTINATION ${phot2pix_CMAKE_DIR}
    )
    write_basic_package_version_file(
        ${GENERATED_DIR}/phot2pixConfigVersion.cmake
        VERSION ${phot2pix_NUMERIC_VERSION}
        COMPATIBILITY AnyNewerVersion
    )
    install(FILES
        ${GENERATED_DIR}/phot2pixConfig.cmake
        ${GENERATED_DIR}/phot2pixConfigVersion.cmake
        DESTINATION ${phot2pix_CMAKE_DIR}
    )

    # CMake never removes files from the install directory: drop the versioned libraries
    # (libphot2pix.0.1.0.dylib, libnisppixelsimulator.so.0.1.0, ...) left behind by previous installs.
    install(CODE "set(phot2pix_KEEP_VERSIONS \"${phot2pix_NUMERIC_VERSION};${VERSION_MAJOR}\")")
    install(CODE [[
        foreach(_p2p_name phot2pix nisppixelsimulator)
            file(GLOB _p2p_libs "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/lib${_p2p_name}.*")
            foreach(_p2p_lib ${_p2p_libs})
                get_filename_component(_p2p_file "${_p2p_lib}" NAME)
                if(_p2p_file MATCHES "^lib${_p2p_name}\\.(so\\.)?([0-9.]*[0-9])(\\.dylib)?$")
                    list(FIND phot2pix_KEEP_VERSIONS "${CMAKE_MATCH_2}" _p2p_keep)
                    if(_p2p_keep EQUAL -1)
                        message(STATUS "Removing old library: ${_p2p_lib}")
                        file(REMOVE "${_p2p_lib}")
                    endif()
                endif()
            endforeach()
        endforeach()
    ]])
endif()

#-------------------- Summary
message("${BoldGreen}----------- phot2pix ${VERSION} configured -----------${ColourReset}")
message("${Green}-- Build directory  :${ColourReset} ${CMAKE_BINARY_DIR}")
message("${Green}-- Build type       :${ColourReset} ${CMAKE_BUILD_TYPE}")
message("${Green}-- AddressSanitizer :${ColourReset} ${WITH_ASAN}")
message("${Green}-- Tests            :${ColourReset} ${PHOT2PIX_TESTS}")
message("${Green}-- phot2pix         :${ColourReset} ${phot2pix_NUMERIC_VERSION}")
message("${Green}-- DSTmath          :${ColourReset} ${DSTmath_VERSION}")
message("${Green}-- DSTfits          :${ColourReset} ${DSTfits_VERSION}")
message("${Green}-- DSTCore          :${ColourReset} ${DSTCore_VERSION}")
if(DSTools_DOWNLOADED)
    message("${Green}-- DSTools headers  :${ColourReset} downloaded (${DSTools_REF})")
else()
    message("${Green}-- DSTools headers  :${ColourReset} ${DSTools_VERSION}")
endif()
if(TARGET nispSimulator)
    message("${Green}-- nispSimulator    :${ColourReset} ON")
else()
    message("${Green}-- nispSimulator    :${ColourReset} OFF (needs DSTCore >= ${DSTCore_MIN_VERSION})")
endif()
message("${Green}-- Install prefix   :${ColourReset} ${CMAKE_INSTALL_PREFIX}")
if(WITH_ASAN)
    message("${Yellow}--                    AddressSanitizer build: installing it replaces the release libphot2pix${ColourReset}")
endif()
message("${BoldGreen}-------------------------------------------------${ColourReset}")
