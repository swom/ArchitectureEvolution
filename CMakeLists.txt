cmake_minimum_required (VERSION 3.8)

# Integrate our local vcpkg install into CMake
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake")
if (WIN32)
    # the default is x86-windows -- 32Bit, really?
    set(VCPKG_TARGET_TRIPLET x64-windows)
endif()

project("arc_evo" C CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Let CMake make the call
# $> git submodule update --init --remote --recursive
# so the user doesn't have to.
#
# Optionally, you can/should add these lines to all the sub-folder CMakeLists.txt.
# That would take care about recursive submodules.
find_package(Git REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --remote --recursive
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                RESULT_VARIABLE GIT_SUBMOD_RESULT)
if (NOT GIT_SUBMOD_RESULT EQUAL "0")
    message(FATAL_ERROR "git submodule update --init --remote --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
endif()

find_package(TBB CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

if(UNIX)
set(CMAKE_BUILD_TYPE "Release")
endif()


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


target_link_libraries(hello_vcpkg PUBLIC
    TBB::tbb
    nlohmann_json::nlohmann_json
)
