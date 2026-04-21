
# ==============================================================================
# OASIS - Offshore Advanced SImulation Software
# Cross-platform CMake build system (Linux / Windows)
# ==============================================================================
CMAKE_MINIMUM_REQUIRED(VERSION 3.14)

# ==============================================================================
# 1. COMPILER SETUP (must be before project() call)
# ==============================================================================
if(UNIX)
  # Linux / cluster: use environment-provided compilers if available
  if(DEFINED ENV{EBROOTGCCCORE})
    set(CMAKE_C_COMPILER   "$ENV{EBROOTGCCCORE}/bin/gcc"  CACHE STRING "")
    set(CMAKE_CXX_COMPILER "$ENV{EBROOTGCCCORE}/bin/g++"  CACHE STRING "")
    set(CMAKE_LINKER        "/usr/bin/ld"                  CACHE STRING "")
  endif()
endif()

# Load vcpkg toolchain from CMakeUserConfig.cmake BEFORE project() so that
# vcpkg's find_package integration is available for all dependencies.
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(_user_config "${CMAKE_CURRENT_SOURCE_DIR}/CMakeUserConfig.cmake")
  if(EXISTS "${_user_config}")
    file(STRINGS "${_user_config}" _tc_lines REGEX "CMAKE_TOOLCHAIN_FILE")
    foreach(_l IN LISTS _tc_lines)
      if(NOT _l MATCHES "^[ \t]*#")  # skip comments
        string(REGEX MATCH "\"([^\"]+)\"" _ "${_l}")
        if(CMAKE_MATCH_1 AND EXISTS "${CMAKE_MATCH_1}")
          set(CMAKE_TOOLCHAIN_FILE "${CMAKE_MATCH_1}" CACHE STRING "vcpkg toolchain")
          message(STATUS "vcpkg toolchain (from CMakeUserConfig): ${CMAKE_MATCH_1}")
          break()
        endif()
      endif()
    endforeach()
  endif()
endif()

# Enable both C and CXX (C needed for FindHDF5 test compilation)
PROJECT(OASIS C CXX)

# ==============================================================================
# 2. USER CONFIGURATION
# ==============================================================================
# Load user-local paths if present (not tracked by git)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeUserConfig.cmake")
  message(STATUS "Loading user configuration from CMakeUserConfig.cmake")
  include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeUserConfig.cmake")
endif()

# Build options
option(OASIS_USE_OPENFAST   "Enable OpenFAST / FASTurbine wind turbine coupling" OFF)
option(OASIS_USE_SUPERLU    "Enable SuperLU sparse solver (via Armadillo)"       ON)
option(OASIS_USE_STL_READER "Enable STL mesh reader for BodyMesh"               ON)

# Auto-detect: enable OpenFAST on Linux cluster if env vars are set
if(DEFINED ENV{EBROOTOpenFAST} AND NOT OASIS_USE_OPENFAST)
  set(OASIS_USE_OPENFAST ON CACHE BOOL "Enable OpenFAST" FORCE)
endif()

# Auto-detect: enable OpenFAST on Windows if install directories exist
if(WIN32 AND NOT OASIS_USE_OPENFAST)
  if(DEFINED OPENFAST_ROOT AND DEFINED FASTURBINE_ROOT
     AND EXISTS "${OPENFAST_ROOT}" AND EXISTS "${FASTURBINE_ROOT}")
    set(OASIS_USE_OPENFAST ON CACHE BOOL "Enable OpenFAST" FORCE)
    message(STATUS "Auto-detected OpenFAST and FASTurbine_wrapper on Windows")
  endif()
endif()

# C++17 on Windows for std::filesystem
if(WIN32)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ==============================================================================
# 3. VERSION FROM GIT
# ==============================================================================
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

include(GetGitRevisionDescription)
git_describe(VERSION --tags --dirty=-d)

string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+(.*)" "\\1" VERSION_SHA1 "${VERSION}")
set(VERSION_SHORT "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/version.cpp.in
               ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

# ==============================================================================
# 4. FIND DEPENDENCIES
# ==============================================================================

# --- Armadillo ---
if(DEFINED ENV{EBROOTArmadillo})
  # Linux cluster: use EasyBuild path
  set(ARMADILLO_INCLUDE_DIR "$ENV{EBROOTArmadillo}/include")
  message(STATUS "Armadillo (EasyBuild): ${ARMADILLO_INCLUDE_DIR}")
elseif(DEFINED ARMADILLO_INCLUDE_DIR)
  # Manual path from CMakeUserConfig.cmake
  message(STATUS "Armadillo (manual): ${ARMADILLO_INCLUDE_DIR}")
else()
  # Try find_package
  find_package(Armadillo QUIET)
  if(ARMADILLO_FOUND)
    set(ARMADILLO_INCLUDE_DIR ${ARMADILLO_INCLUDE_DIRS})
    message(STATUS "Armadillo (found): ${ARMADILLO_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR "Armadillo not found. Set ARMADILLO_INCLUDE_DIR in CMakeUserConfig.cmake")
  endif()
endif()

# --- HDF5 ---
if(WIN32)
  find_package(HDF5 COMPONENTS C CXX QUIET)
  if(HDF5_FOUND AND HDF5_C_LIBRARIES AND HDF5_CXX_LIBRARIES)
    message(STATUS "HDF5 (found): ${HDF5_INCLUDE_DIRS}")
  elseif(DEFINED HDF5_INCLUDE_DIR AND DEFINED HDF5_LIBRARIES)
    # Manual path from CMakeUserConfig.cmake
    set(HDF5_FOUND TRUE)
    set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
    message(STATUS "HDF5 (manual): ${HDF5_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR
      "HDF5 with C++ bindings not found.\n"
      "  If using vcpkg, install with:  vcpkg install hdf5[cpp]:x64-windows\n"
      "  Or set HDF5_INCLUDE_DIR and HDF5_LIBRARIES in CMakeUserConfig.cmake")
  endif()
endif()

# --- BLAS/LAPACK (OpenBLAS on Linux, find automatically on Windows) ---
if(DEFINED ENV{EBROOTOPENBLAS})
  set(BLAS_INCLUDE_DIR "$ENV{EBROOTOPENBLAS}/include")
  message(STATUS "OpenBLAS (EasyBuild): ${BLAS_INCLUDE_DIR}")
elseif(WIN32)
  # On Windows with vcpkg or Intel MKL, find_package should work
  find_package(BLAS QUIET)
  find_package(LAPACK QUIET)
  if(BLAS_FOUND)
    message(STATUS "BLAS found: ${BLAS_LIBRARIES}")
  elseif(DEFINED BLAS_INCLUDE_DIR)
    message(STATUS "BLAS (manual): ${BLAS_INCLUDE_DIR}")
  else()
    message(STATUS "BLAS not explicitly found - Armadillo may handle it internally")
  endif()
endif()

# --- SuperLU ---
if(OASIS_USE_SUPERLU)
  if(WIN32)
    # Only look for Release config to avoid vcpkg debug lib errors
    set(_superlu_found FALSE)
    if(DEFINED SUPERLU_LIBRARIES)
      set(_superlu_found TRUE)
    else()
      # Manual discovery instead of find_package(superlu CONFIG) to avoid
      # vcpkg targets file failing when debug libs are not installed.
      find_library(SUPERLU_LIBRARIES NAMES superlu)
      find_path(SUPERLU_INCLUDE_DIR NAMES slu_ddefs.h PATH_SUFFIXES superlu)
      if(SUPERLU_LIBRARIES)
        set(_superlu_found TRUE)
        message(STATUS "SuperLU library: ${SUPERLU_LIBRARIES}")
        if(SUPERLU_INCLUDE_DIR)
          message(STATUS "SuperLU headers: ${SUPERLU_INCLUDE_DIR}")
        endif()
      endif()
    endif()
    if(NOT _superlu_found)
      message(WARNING "SuperLU not found on Windows. Disabling OASIS_USE_SUPERLU. "
                      "Set SUPERLU_INCLUDE_DIR and SUPERLU_LIBRARIES in CMakeUserConfig.cmake to enable.")
      set(OASIS_USE_SUPERLU OFF)
    endif()
  endif()
endif()

# --- OpenMP ---
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found")
else()
  message(STATUS "OpenMP not found - building without parallelization")
endif()

# --- OpenFAST / FASTurbine wrapper ---
if(OASIS_USE_OPENFAST)
  if(WIN32)
    # Windows: search for OpenFAST and FASTurbine wrapper (compiled with Intel oneAPI Fortran)
    # Paths are set in CMakeUserConfig.cmake (OPENFAST_ROOT, FASTURBINE_ROOT, INTEL_ONEAPI_ROOT)
    if(EXISTS "${OPENFAST_ROOT}/include" AND EXISTS "${FASTURBINE_ROOT}/include")
      set(OPENFAST_INCLUDE_DIR "${OPENFAST_ROOT}/include")
      set(FASTURBINE_INCLUDE_DIR "${FASTURBINE_ROOT}/include")
      
      # Find required OpenFAST libraries
      find_library(AERODYN_LIB aerodynlib PATHS "${OPENFAST_ROOT}/lib" NO_DEFAULT_PATH)
      find_library(IFW_LIB ifwlib PATHS "${OPENFAST_ROOT}/lib" NO_DEFAULT_PATH)
      find_library(SERVODYN_LIB servodynlib PATHS "${OPENFAST_ROOT}/lib" NO_DEFAULT_PATH)
      find_library(NWTC_LIB nwtclibs PATHS "${OPENFAST_ROOT}/lib" NO_DEFAULT_PATH)
      find_library(FASTURBINE_LIB fasturbwrapper PATHS "${FASTURBINE_ROOT}/lib" NO_DEFAULT_PATH)
      
      # Find Intel oneAPI Fortran runtime (required for linking OpenFAST Fortran libs with MSVC)
      find_library(IFCORE_LIB libifcoremd PATHS "${INTEL_ONEAPI_ROOT}/lib" NO_DEFAULT_PATH)
      find_library(LIBMMD_LIB libmmd PATHS "${INTEL_ONEAPI_ROOT}/lib" NO_DEFAULT_PATH)
      
      if(AERODYN_LIB AND IFW_LIB AND SERVODYN_LIB AND NWTC_LIB AND FASTURBINE_LIB)
        set(OPENFAST_LIBRARIES ${FASTURBINE_LIB} ${AERODYN_LIB} ${IFW_LIB} ${SERVODYN_LIB} ${NWTC_LIB})
        
        # Add Intel Fortran runtime if found
        if(IFCORE_LIB AND LIBMMD_LIB)
          list(APPEND OPENFAST_LIBRARIES ${IFCORE_LIB} ${LIBMMD_LIB})
          message(STATUS "Intel Fortran runtime: ${INTEL_ONEAPI_ROOT}/lib")
        else()
          message(WARNING "Intel oneAPI Fortran runtime not found at ${INTEL_ONEAPI_ROOT}/lib. "
                          "You may need to run from Intel oneAPI command prompt or set INTEL_ONEAPI_ROOT in CMakeUserConfig.cmake.")
        endif()
        
        message(STATUS "OpenFAST (Windows/Intel): ${OPENFAST_INCLUDE_DIR}")
        message(STATUS "FASTurbine wrapper (Windows): ${FASTURBINE_INCLUDE_DIR}")
      else()
        message(FATAL_ERROR "OpenFAST libraries not found in ${OPENFAST_ROOT}/lib")
      endif()
    else()
      message(FATAL_ERROR "OASIS_USE_OPENFAST=ON but OpenFAST/FASTurbine not found. "
                          "Set OPENFAST_ROOT and FASTURBINE_ROOT in CMakeUserConfig.cmake")
    endif()
  elseif(DEFINED ENV{EBROOTOpenFAST})
    set(OPENFAST_INCLUDE_DIR "$ENV{EBROOTOpenFAST}/install/include")
    message(STATUS "OpenFAST (EasyBuild): ${OPENFAST_INCLUDE_DIR}")
  elseif(DEFINED OPENFAST_INCLUDE_DIR)
    message(STATUS "OpenFAST (manual): ${OPENFAST_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR "OASIS_USE_OPENFAST=ON but OpenFAST not found. "
                        "Set OPENFAST_INCLUDE_DIR in CMakeUserConfig.cmake")
  endif()
endif()

# --- stl_reader ---
if(OASIS_USE_STL_READER)
  # EasyBuild: try both naming conventions (lowercase and uppercase)
  if(DEFINED ENV{EBROOTstl_reader})
    set(_STL_READER_ROOT "$ENV{EBROOTstl_reader}")
  elseif(DEFINED ENV{EBROOTSTL_READER})
    set(_STL_READER_ROOT "$ENV{EBROOTSTL_READER}")
  elseif(DEFINED ENV{EBROOTSTL_MINUS_READER})
    set(_STL_READER_ROOT "$ENV{EBROOTSTL_MINUS_READER}")
  endif()

# --- yaml-cpp ---
if(WIN32)
  find_package(yaml-cpp CONFIG QUIET)
  if(yaml-cpp_FOUND)
    message(STATUS "yaml-cpp (found via vcpkg/config)")
  elseif(DEFINED YAMLCPP_INCLUDE_DIR AND DEFINED YAMLCPP_LIBRARIES)
    set(yaml-cpp_FOUND TRUE)
    message(STATUS "yaml-cpp (manual): ${YAMLCPP_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR
      "yaml-cpp not found.\n"
      "  If using vcpkg, install with:  vcpkg install yaml-cpp:x64-windows\n"
      "  Or set YAMLCPP_INCLUDE_DIR and YAMLCPP_LIBRARIES in CMakeUserConfig.cmake")
  endif()
endif()

  if(DEFINED _STL_READER_ROOT)
    # Search for header in root and include subdirectory
    find_path(STL_READER_INCLUDE_DIR stl_reader.h
      PATHS "${_STL_READER_ROOT}" "${_STL_READER_ROOT}/include"
      NO_DEFAULT_PATH
    )
    if(STL_READER_INCLUDE_DIR)
      message(STATUS "stl_reader (EasyBuild): ${STL_READER_INCLUDE_DIR}")
    else()
      message(WARNING "stl_reader module loaded (EBROOT=${_STL_READER_ROOT}) but stl_reader.h not found there. Disabling.")
      set(OASIS_USE_STL_READER OFF)
    endif()
  elseif(STL_READER_INCLUDE_DIR AND EXISTS "${STL_READER_INCLUDE_DIR}/stl_reader.h")
    message(STATUS "stl_reader (manual): ${STL_READER_INCLUDE_DIR}")
  else()
    # stl_reader is header-only, try paths from CMakeUserConfig.cmake then common locations
    find_path(STL_READER_INCLUDE_DIR stl_reader.h
      PATHS
        ${STL_READER_SEARCH_PATHS}
    )
    if(STL_READER_INCLUDE_DIR)
      message(STATUS "stl_reader found: ${STL_READER_INCLUDE_DIR}")
    else()
      message(WARNING "stl_reader not found. Disabling OASIS_USE_STL_READER. "
                      "Set STL_READER_INCLUDE_DIR in CMakeUserConfig.cmake or place stl_reader.h in third_party/stl_reader/")
      set(OASIS_USE_STL_READER OFF)
    endif()
  endif()
endif()

# ==============================================================================
# 5. SOURCE FILES
# ==============================================================================

# Core source files (always compiled)
set(OASIS_SOURCE_FILES 
    "src/CommonTools.cpp" 
    "src/Logger.cpp"
    "src/main.cpp"
    "src/MathTools.cpp"
    "src/os_tools.cpp"
    "src/BCPs/BCPs.cpp"
    "src/BCPs/Winchies.cpp"
    "src/BCPs/WinchiesController.cpp"
    "src/BCPs/WinchiesControllerConstantTension.cpp"
    "src/BCPs/WinchiesControllerHorizontal.cpp"
    "src/Bodies/Bodies.cpp"
    "src/Bodies/BodyMesh.cpp"
    "src/Exceptions/Exception.cpp"
    "src/Waves/Wave.cpp"
    "src/Simulations/Simulation.cpp"
    "src/Hydro/Morison.cpp"
    "src/Hydro/HydroDatabase.cpp"
    "src/Lines/Lines_Dyn.cpp"
    "src/Lines/Lines_QS.cpp"
    "src/Lines/LineType.cpp"
    "src/ODE_solvers/ODE_solvers.cpp"
    "src/SEM_math/quadrule.cpp"
    "src/Spring/Spring.cpp"
    "src/Spring/SpringType.cpp"
    "src/Sinking/Sinking.cpp"
    "src/SeaFloor/SeaFloor.cpp"
    "src/OWC/OWC.cpp"
)

# Core header files
set(OASIS_HEADER_FILES 
    "src/CommonTools.hpp" 
    "src/Logger.hpp"
    "src/main.cpp"
    "src/MathTools.hpp"
    "src/os_tools.hpp"
    "src/version.hpp"
    "src/BCPs/BCPs.hpp"
    "src/BCPs/Winchies.hpp"
    "src/BCPs/WinchiesController.hpp"
    "src/BCPs/WinchiesControllerConstantTension.hpp"
    "src/BCPs/WinchiesControllerHorizontal.hpp"
    "src/Bodies/Bodies.hpp"
    "src/Bodies/BodyMesh.hpp"
    "src/Exceptions/Exception.hpp"
    "src/Hydro/Morison.hpp"
    "src/Hydro/HydroDatabase.hpp"
    "src/Lines/Lines.hpp"
    "src/ODE_solvers/ODE_solvers.hpp"
    "src/SEM_math/quadrule.hpp"
    "src/Waves/Wave.hpp"
    "src/Simulations/ISimulation.hpp"
    "src/Simulations/Simulation.hpp"
    "src/Spring/Spring.hpp"
    "src/Sinking/Sinking.hpp"
    "src/SeaFloor/SeaFloor.hpp"
    "src/OWC/OWC.hpp"
)

# WindTurbine header is always included (has stub class when OpenFAST is off)
list(APPEND OASIS_HEADER_FILES "src/WindTurbine/WindTurbine.hpp")

# WindTurbine source only compiled with OpenFAST (entire file is #ifdef guarded)
if(OASIS_USE_OPENFAST)
  list(APPEND OASIS_SOURCE_FILES "src/WindTurbine/WindTurbine.cpp")
endif()

set(OASIS_FILES ${OASIS_SOURCE_FILES} ${OASIS_HEADER_FILES})

# ==============================================================================
# 6. COMPILE DEFINITIONS
# ==============================================================================
set(OASIS_DEFINITIONS -DARMA_DONT_USE_WRAPPER -DARMA_USE_HDF5 -DARMA_DONT_PRINT_FAST_MATH_WARNING)

if(OASIS_USE_SUPERLU)
  list(APPEND OASIS_DEFINITIONS -DARMA_USE_SUPERLU)
  # Armadillo includes SuperLU headers as: #include "ARMA_SUPERLU_INCLUDE_DIR/slu_ddefs.h"
  # Point it to the directory containing the slu_*.h files (with trailing /).
  if(DEFINED SUPERLU_INCLUDE_DIR)
    list(APPEND OASIS_DEFINITIONS "-DARMA_SUPERLU_INCLUDE_DIR=${SUPERLU_INCLUDE_DIR}/")
  endif()
endif()

if(OASIS_USE_OPENFAST)
  list(APPEND OASIS_DEFINITIONS -DOASIS_USE_OPENFAST)
endif()

if(OASIS_USE_STL_READER)
  list(APPEND OASIS_DEFINITIONS -DOASIS_USE_STL_READER)
endif()

if(WIN32)
  list(APPEND OASIS_DEFINITIONS -DH5_BUILT_AS_DYNAMIC_LIB)
endif()

ADD_DEFINITIONS(${OASIS_DEFINITIONS})

# ==============================================================================
# 7. BUILD TARGET
# ==============================================================================
ADD_EXECUTABLE(${PROJECT_NAME} ${OASIS_FILES} ${version_file})

# --- Optimization flags ---
target_compile_options(${PROJECT_NAME} PRIVATE
  $<$<AND:$<CONFIG:Release>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>>:-O3 -march=native -ffast-math>
  $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:/O2 /fp:fast>
)

# --- OpenMP ---
if(OpenMP_CXX_FOUND)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

# --- LTO (Release only) ---
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_output)
if(ipo_supported)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
else()
  message(STATUS "LTO not supported: ${ipo_output}")
endif()

# --- Post-build: copy non-vcpkg runtime DLLs on Windows ---
# vcpkg handles its own DLLs via VCPKG_APPLOCAL_DEPS (on by default).
# Here we copy FASTurbine and Intel Fortran runtime DLLs so that
# build/Release/ contains everything needed to run.
if(WIN32 AND OASIS_USE_OPENFAST)
  if(EXISTS "${FASTURBINE_ROOT}/lib/fasturbwrapper.dll")
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${FASTURBINE_ROOT}/lib/fasturbwrapper.dll"
        $<TARGET_FILE_DIR:${PROJECT_NAME}>
      COMMENT "Copying fasturbwrapper.dll")
  endif()
  if(DEFINED INTEL_ONEAPI_ROOT)
    foreach(_dll libifcoremd.dll libifportmd.dll libmmd.dll svml_dispmd.dll)
      if(EXISTS "${INTEL_ONEAPI_ROOT}/bin/${_dll}")
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${INTEL_ONEAPI_ROOT}/bin/${_dll}"
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
          COMMENT "Copying ${_dll}")
      endif()
    endforeach()
  endif()
endif()

# ==============================================================================
# 8. INCLUDE DIRECTORIES
# ==============================================================================
set(OASIS_INCLUDE_DIRS "")
list(APPEND OASIS_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIR})

if(DEFINED BLAS_INCLUDE_DIR)
  list(APPEND OASIS_INCLUDE_DIRS ${BLAS_INCLUDE_DIR})
endif()

if(WIN32 AND HDF5_FOUND)
  list(APPEND OASIS_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
endif()

if(OASIS_USE_OPENFAST AND DEFINED OPENFAST_INCLUDE_DIR)
  list(APPEND OASIS_INCLUDE_DIRS ${OPENFAST_INCLUDE_DIR})
  if(DEFINED FASTURBINE_INCLUDE_DIR)
    list(APPEND OASIS_INCLUDE_DIRS ${FASTURBINE_INCLUDE_DIR})
  endif()
endif()

if(OASIS_USE_STL_READER AND STL_READER_INCLUDE_DIR)
  list(APPEND OASIS_INCLUDE_DIRS ${STL_READER_INCLUDE_DIR})
endif()

if(OASIS_USE_SUPERLU AND DEFINED SUPERLU_INCLUDE_DIR)
  list(APPEND OASIS_INCLUDE_DIRS ${SUPERLU_INCLUDE_DIR})
endif()

if(WIN32 AND yaml-cpp_FOUND AND DEFINED YAMLCPP_INCLUDE_DIR)
  list(APPEND OASIS_INCLUDE_DIRS ${YAMLCPP_INCLUDE_DIR})
endif()

TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC ${OASIS_INCLUDE_DIRS})

# ==============================================================================
# 9. LINK LIBRARIES
# ==============================================================================
if(WIN32)
  # ---- Windows linking ----
  if(HDF5_FOUND)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${HDF5_LIBRARIES})
    if(DEFINED HDF5_CXX_LIBRARIES)
      TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${HDF5_CXX_LIBRARIES})
    endif()
  endif()
  
  if(BLAS_FOUND)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${BLAS_LIBRARIES})
  endif()
  if(LAPACK_FOUND)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES})
  endif()
  
  if(OASIS_USE_SUPERLU AND DEFINED SUPERLU_LIBRARIES)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${SUPERLU_LIBRARIES})
    # SuperLU (vcpkg) may be compiled with METIS support
    find_library(METIS_LIBRARY NAMES metis)
    if(METIS_LIBRARY)
      TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${METIS_LIBRARY})
      message(STATUS "METIS (SuperLU dependency): ${METIS_LIBRARY}")
    endif()
  endif()
  
  if(OASIS_USE_OPENFAST AND DEFINED OPENFAST_LIBRARIES)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${OPENFAST_LIBRARIES})
  endif()

  if(yaml-cpp_FOUND)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE yaml-cpp::yaml-cpp)
  elseif(DEFINED YAMLCPP_LIBRARIES)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${YAMLCPP_LIBRARIES})
  endif()

  # Additional manual Windows libraries from CMakeUserConfig.cmake
  if(DEFINED OASIS_EXTRA_LIBRARIES)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${OASIS_EXTRA_LIBRARIES})
  endif()

else()
  # ---- Linux linking ----
  TARGET_LINK_LIBRARIES(${PROJECT_NAME}
    PRIVATE
    -lpthread
    -lstdc++fs
    -lhdf5
    -lhdf5_cpp
    -lopenblas
  )
  
  if(OASIS_USE_SUPERLU)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE -lsuperlu)
  endif()

  TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE -lyaml-cpp)

  if(OASIS_USE_OPENFAST)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME}
      PRIVATE
      -lgfortran
      -lquadmath
      -laerodynlib
      -lifwlib
      -lservodynlib
      -lnwtclibs
      -lfasturbwrapper
    )
    # OpenMPI link flags (cluster-specific, paths from CMakeUserConfig.cmake)
    if(DEFINED OPENMPI_ROOT AND EXISTS "${OPENMPI_ROOT}")
      set_target_properties(${PROJECT_NAME} PROPERTIES
        COMPILE_FLAGS -I${OPENMPI_ROOT}/include)
      set_target_properties(${PROJECT_NAME} PROPERTIES
        LINK_FLAGS "-L${HWLOC_LIB_DIR} \
        -L${LIBEVENT_LIB_DIR} \
        -Wl,-rpath \
        -Wl,${HWLOC_LIB_DIR} \
        -Wl,-rpath \
        -Wl,${LIBEVENT_LIB_DIR} \
        -Wl,-rpath \
        -Wl,${OPENMPI_ROOT}/lib \
        -Wl,--enable-new-dtags \
        -L${OPENMPI_ROOT}/lib \
        -lmpi")
    endif()
  endif()
endif()

set_property(TARGET ${PROJECT_NAME} PROPERTY LINKER_LANGUAGE CXX)

# ==============================================================================
# 10. STATUS SUMMARY
# ==============================================================================
message(STATUS "")
message(STATUS "=== OASIS Build Configuration ===")
message(STATUS "  Platform:       ${CMAKE_SYSTEM_NAME}")
message(STATUS "  Compiler:       ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Build type:     ${CMAKE_BUILD_TYPE}")
message(STATUS "  OpenFAST:       ${OASIS_USE_OPENFAST}")
message(STATUS "  SuperLU:        ${OASIS_USE_SUPERLU}")
message(STATUS "  STL reader:     ${OASIS_USE_STL_READER}")
message(STATUS "  OpenMP:         ${OpenMP_CXX_FOUND}")
message(STATUS "  Armadillo:      ${ARMADILLO_INCLUDE_DIR}")
if(WIN32 AND HDF5_FOUND)
  message(STATUS "  HDF5:           ${HDF5_INCLUDE_DIRS}")
endif()
message(STATUS "=================================")
