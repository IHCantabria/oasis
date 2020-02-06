
PROJECT(GIT_CONTROL C CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 3.9)

# Appends the cmake/modules path to MAKE_MODULE_PATH variable.
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

# Get Git Version
include(GetGitRevisionDescription)
git_describe(VERSION --tags --dirty=-d)

#parse the version information into pieces.
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+(.*)" "\\1" VERSION_SHA1 "${VERSION}")
set(VERSION_SHORT "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

MESSAGE("GIT_EXECUTABLE: ${GIT_EXECUTABLE}")
MESSAGE("${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE("VERSION: ${VERSION}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/version.cpp.in
                ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

# Add source files
set(SRC_FILES "src/hello.cpp" "src/version.hpp")

MESSAGE(STATUS "DDDDDD")
MESSAGE(STATUS $ENV{EASYBUILD_SOFTWARE_PREFIX})
MESSAGE(STATUS "$ENV{EASYBUILD_SOFTWARE_PREFIX}/Armadillo")

#Add the version_file to the project build or it won't compile.
ADD_DEFINITIONS(-DARMA_DONT_USE_WRAPPER)
ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_FILES} ${version_file})
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:$ENV{EASYBUILD_SOFTWARE_PREFIX}/Armadillo/9.800.4/include>" "$<BUILD_INTERFACE:$ENV{EASYBUILD_SOFTWARE_PREFIX}/OpenBLAS/0.2.19-GCC-6.3.0-2.27-LAPACK-3.7.0/include>")
TARGET_LINK_LIBRARIES(${PROJECT_NAME} "$ENV{EASYBUILD_SOFTWARE_PREFIX}/OpenBLAS/0.2.19-GCC-6.3.0-2.27-LAPACK-3.7.0/lib/libopenblas.a")