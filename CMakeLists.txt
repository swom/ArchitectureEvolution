
cmake_minimum_required (VERSION 3.8)

project("arc_evo_CMake" C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED)
set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
if (WIN32)
  set(CMAKE_CXX_FLAGS_RELEASE "-DNOMINMAX -DNDEBUG ${RELEASE_WARNING_FLAGS} ${WARNING_FLAGS} ${OTHER_FLAGS}")
else()
	set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG ${RELEASE_WARNING_FLAGS} ${WARNING_FLAGS} ${OTHER_FLAGS} -O2 -ffast-math -funroll-loops -Wno-ignored-attributes")
endif()

file(GLOB src
     "${PROJECT_SOURCE_DIR}/*.cpp"
     "${PROJECT_SOURCE_DIR}/*.hpp"
     "${PROJECT_SOURCE_DIR}/*.h"
)

add_compile_options(-Wa,-mbig)
add_executable("arc_evo_cmake" ${src})
include_directories(${CMAKE_SOURCE_DIR}/libs/ $(CMAKE_SYSTEM_PREFIX_PATH}/include/))
install(TARGETS arc_evo_cmake DESTINATION bin)
install(TARGETS arc_evo_cmake CONFIGURATIONS Release DESTINATION bin)
