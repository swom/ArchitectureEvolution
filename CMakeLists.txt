cmake_minimum_required (VERSION 3.8)

project("test" C CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

file(GLOB src
     "${PROJECT_SOURCE_DIR}/*.cpp"
     "${PROJECT_SOURCE_DIR}/*.hpp"
     "${PROJECT_SOURCE_DIR}/*.h"
)
if(MINGW)
  add_compile_options(-Wa,-mbig)
endif()


add_executable("test" ${src})
