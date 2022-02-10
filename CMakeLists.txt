cmake_minimum_required (VERSION 3.8)

project("arc_evo" C CXX)
set(CMAKE_CXX_STANDARD 20)

if(UNIX)
set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

file(GLOB src
     "${PROJECT_SOURCE_DIR}/*.cpp"
     "${PROJECT_SOURCE_DIR}/*.hpp"
     "${PROJECT_SOURCE_DIR}/*.h"
)
if(MINGW)
  add_compile_options(-Wa,-mbig)
endif()


add_executable("arc_evo" ${src})
