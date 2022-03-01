
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
     "${PROJECT_SOURCE_DIR}/src/*.cpp"
     "${PROJECT_SOURCE_DIR}/src/*.h"
     "${PROJECT_SOURCE_DIR}/src/*.hpp"
)

add_executable(arc_evo ${src})
set_target_properties(arc_evo PROPERTIES LINKER_LANGUAGE CXX)


target_link_libraries(arc_evo PUBLIC
    TBB::tbb
    nlohmann_json::nlohmann_json
)
