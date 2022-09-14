find_package(OpenMP)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(cxxopts CONFIG REQUIRED)
find_package(range-v3 CONFIG REQUIRED)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  add_compile_options(-march=native)
  # add_compile_options(-Weffc++)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # MSVC
  add_compile_options(/GS- /bigobj)     # /GS-: no security cookie
endif()

if(MINGW)
  add_compile_options(-Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Werror)
  add_compile_options(-Wa,-mbig-obj)
  add_compile_options("$<$<OR:$<CONFIG:Release>,$<CONFIG:MinSizeRel>>:-Ofast>")
  # add_compile_options(-Weffc++)
endif()


file(GLOB src
  "${PROJECT_SOURCE_DIR}/src/*.cpp"
  "${PROJECT_SOURCE_DIR}/src/*.h"
  "${PROJECT_SOURCE_DIR}/src/*.hpp"
  )

add_executable(arc_evo ${src})

target_link_libraries(arc_evo PUBLIC
  nlohmann_json::nlohmann_json
  cxxopts::cxxopts
  range-v3
  )

if(OpenMP_CXX_FOUND)
    target_link_libraries(arc_evo PUBLIC OpenMP::OpenMP_CXX)
endif()
