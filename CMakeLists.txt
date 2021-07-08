
# Define project and minimum cmake version required
PROJECT(OASIS CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 3.7)

# Set the compiler
set(CMAKE_C_COMPILER $ENV{EBROOTGCCCORE}/bin/gcc)
set(CMAKE_CXX_COMPILER $ENV{EBROOTGCCCORE}/bin/g++)
set(CMAKE_LINKER /usr/bin/ld)

# Appends the cmake/modules path to MAKE_MODULE_PATH variable.
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

# Get Git Version
include(GetGitRevisionDescription)
git_describe(VERSION --tags --dirty=-d)

# Parse the version information into pieces
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+(.*)" "\\1" VERSION_SHA1 "${VERSION}")
set(VERSION_SHORT "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# Configure version file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/version.cpp.in
               ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

# Add target link libraries
set(EXTERNAL_LIBS_INCLUDE_PATH "$<BUILD_INTERFACE:$ENV{EBROOTArmadillo}/include>"
"$<BUILD_INTERFACE:$ENV{EBROOTOPENBLAS}/include>" "$<BUILD_INTERFACE:$ENV{EBROOTOpenFAST}/install/include>")

# Add source files
set(OASIS_SOURCE_FILES 
    "src/CommonTools.cpp" 
    "src/main.cpp"
    "src/MathTools.cpp"
    "src/os_tools.cpp"
    "src/BCPs/BCPs.cpp"
    "src/BCPs/Winchies.cpp"
    "src/BCPs/WinchiesController.cpp"
    "src/Bodies/Bodies.cpp"
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
    "src/WindTurbine/WindTurbine.cpp"
    )

# Add include files
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
    "src/Exceptions/Exception.hpp"
    "src/Hydro/Morison.hpp"
    "src/Hydro/HydroDatabase.hpp"
    "src/Lines/Lines.hpp"
    "src/ODE_solvers/ODE_solvers.hpp"
    "src/SEM_math/quadrule.hpp"
    "src/Waves/Wave.hpp"
    "src/Simulations/Simulation.hpp"
    "src/Spring/Spring.hpp"
    "src/Sinking/Sinking.hpp"
    "src/WindTurbine/WindTurbine.hpp"
    )

# Join source and header files
set(OASIS_FILES ${OASIS_SOURCE_FILES} ${OASIS_HEADER_FILES})
# Add the compilation flags: HDF5 is a dynamic library!
ADD_DEFINITIONS(-DARMA_DONT_USE_WRAPPER -DARMA_USE_HDF5 -DARMA_USE_SUPERLU)
# Add the version_file to the project build or it won't compile.
ADD_EXECUTABLE(${PROJECT_NAME} ${OASIS_FILES} ${version_file})
# Target to the include directories
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC ${EXTERNAL_LIBS_INCLUDE_PATH})
# Target to the link libraries
TARGET_LINK_LIBRARIES(${PROJECT_NAME} -lgfortran -lquadmath -lpthread -lstdc++fs -lhdf5 -lsuperlu -lopenblas -lopenfastcpplib)

