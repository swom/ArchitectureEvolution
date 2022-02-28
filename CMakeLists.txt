cmake_minimum_required (VERSION 3.8)

project("arc_evo" C CXX)
set(CMAKE_CXX_STANDARD 20)

if(UNIX)
set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# the rather controversial -Weffc++ compiler flag is supported
# by GNU g++ and clang++ but not by MSVC
# since we can't pass unsupported flags to the compiler,
# we need to identify the compiler:
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-march=native)
    # add_compile_options(-Weffc++)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # MSVC
    add_compile_options(/GS-)     # /GS-: no security cookie
endif()


if(MINGW)
  add_compile_options(-Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Werror -g -Wa,-mbig)
  add_compile_options("$<$<OR:$<CONFIG:Release>,$<CONFIG:MinSizeRel>>:-Ofast>")
endif()

file(GLOB src
     "${PROJECT_SOURCE_DIR}/*.cpp"
     "${PROJECT_SOURCE_DIR}/*.h"
     "${PROJECT_SOURCE_DIR}/*.hpp"
)

add_executable("arc_evo" ${src})



