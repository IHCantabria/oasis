
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
if(DEFINED ENV{EBROOTOpenFAST} AND NOT DEFINED OASIS_USE_OPENFAST)
  set(OASIS_USE_OPENFAST ON CACHE BOOL "Enable OpenFAST" FORCE)
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
  if(HDF5_FOUND)
    message(STATUS "HDF5 (found): ${HDF5_INCLUDE_DIRS}")
  elseif(DEFINED HDF5_INCLUDE_DIR AND DEFINED HDF5_LIBRARIES)
    # Manual path from CMakeUserConfig.cmake
    set(HDF5_FOUND TRUE)
    set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
    message(STATUS "HDF5 (manual): ${HDF5_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR "HDF5 not found. Set HDF5_INCLUDE_DIR and HDF5_LIBRARIES in CMakeUserConfig.cmake")
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
      # Try find_package with error handling for missing debug libs
      set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)
      find_package(superlu CONFIG QUIET)
      if(superlu_FOUND)
        set(_superlu_found TRUE)
        # Get the imported target's release library
        if(TARGET superlu::superlu)
          get_target_property(_superlu_lib superlu::superlu IMPORTED_LOCATION_RELEASE)
          if(_superlu_lib)
            set(SUPERLU_LIBRARIES ${_superlu_lib})
          endif()
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
  if(DEFINED ENV{EBROOTOpenFAST})
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
  if(DEFINED ENV{EBROOTstl_reader})
    set(STL_READER_INCLUDE_DIR "$ENV{EBROOTstl_reader}/include")
    message(STATUS "stl_reader (EasyBuild): ${STL_READER_INCLUDE_DIR}")
  elseif(DEFINED STL_READER_INCLUDE_DIR)
    message(STATUS "stl_reader (manual): ${STL_READER_INCLUDE_DIR}")
  else()
    # stl_reader is header-only, try common locations
    find_path(STL_READER_INCLUDE_DIR stl_reader.h
      PATHS
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/stl_reader
        ${CMAKE_CURRENT_SOURCE_DIR}/external/stl_reader
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
    "src/main.cpp"
    "src/MathTools.cpp"
    "src/os_tools.cpp"
    "src/BCPs/BCPs.cpp"
    "src/BCPs/Winchies.cpp"
    "src/BCPs/WinchiesController.cpp"
    "src/Bodies/Bodies.cpp"
    "src/Bodies/BodyMesh.cpp"
    "src/Exceptions/Exception.cpp"
    "src/Waves/Wave.cpp"
    "src/Simulations/Simulation.cpp"
    "src/Hydro/Morison.cpp"
    "src/Hydro/HydroDatabase.cpp"
    "src/Lines/Lines_Dyn.cpp"
    "src/Lines/Lines_QS.cpp"
    "src/ODE_solvers/ODE_solvers.cpp"
    "src/SEM_math/quadrule.cpp"
    "src/Spring/Spring.cpp"
    "src/Sinking/Sinking.cpp"
    "src/SeaFloor/SeaFloor.cpp"
    "src/OWC/OWC.cpp"
)

# Core header files
set(OASIS_HEADER_FILES 
    "src/CommonTools.hpp" 
    "src/main.cpp"
    "src/MathTools.hpp"
    "src/os_tools.hpp"
    "src/version.hpp"
    "src/BCPs/BCPs.hpp"
    "src/BCPs/Winchies.hpp"
    "src/BCPs/WinchiesController.hpp"
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
set(OASIS_DEFINITIONS -DARMA_DONT_USE_WRAPPER -DARMA_USE_HDF5)

if(OASIS_USE_SUPERLU)
  list(APPEND OASIS_DEFINITIONS -DARMA_USE_SUPERLU)
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
  $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang>>:-O3 -march=native -ffast-math>
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
endif()

if(OASIS_USE_STL_READER AND DEFINED STL_READER_INCLUDE_DIR)
  list(APPEND OASIS_INCLUDE_DIRS ${STL_READER_INCLUDE_DIR})
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
    # OpenMPI link flags (cluster-specific)
    if(EXISTS "/sNow/easybuild")
      set_target_properties(${PROJECT_NAME} PROPERTIES
        COMPILE_FLAGS -I/sNow/easybuild/centos/7.3.1611/SandyBridge/software/OpenMPI/4.0.5-GCC-10.2.0/include)
      set_target_properties(${PROJECT_NAME} PROPERTIES
        LINK_FLAGS "-L/sNow/easybuild/centos/7.3.1611/SandyBridge/software/hwloc/2.2.0-GCCcore-10.2.0/lib \
        -L/sNow/easybuild/centos/7.3.1611/SandyBridge/software/libevent/2.1.12-GCCcore-10.2.0/lib64 \
        -Wl,-rpath \
        -Wl,/sNow/easybuild/centos/7.3.1611/SandyBridge/software/hwloc/2.2.0-GCCcore-10.2.0/lib \
        -Wl,-rpath \
        -Wl,/sNow/easybuild/centos/7.3.1611/SandyBridge/software/libevent/2.1.12-GCCcore-10.2.0/lib64 \
        -Wl,-rpath \
        -Wl,/sNow/easybuild/centos/7.3.1611/SandyBridge/software/OpenMPI/4.0.5-GCC-10.2.0/lib \
        -Wl,--enable-new-dtags \
        -L/sNow/easybuild/centos/7.3.1611/SandyBridge/software/OpenMPI/4.0.5-GCC-10.2.0/lib \
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
