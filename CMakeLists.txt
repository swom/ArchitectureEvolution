cmake_minimum_required (VERSION 3.8)

project("arc_evo" C CXX)
set(CMAKE_CXX_STANDARD 20)

if(UNIX)
set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(MINGW)
  add_compile_options(-Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Werror -g -Wa,-mbig)
  add_compile_options("$<$<OR:$<CONFIG:Release>,$<CONFIG:MinSizeRel>>:-Ofast>")
 # add_compile_options(-Weffc++)
endif()

file(GLOB src
     "${PROJECT_SOURCE_DIR}/*.cpp"
     "${PROJECT_SOURCE_DIR}/*.h"
     "${PROJECT_SOURCE_DIR}/*.hpp"
)

add_executable("arc_evo" ${src})
