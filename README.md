# ArchitectureEvolution

## CMakeLists.txt

```cmake
cmake_minimum_required (VERSION 3.8)

project("test" C CXX)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

file(GLOB src
     "${PROJECT_SOURCE_DIR}/*.cpp"
     "${PROJECT_SOURCE_DIR}/*.hpp"
     "${PROJECT_SOURCE_DIR}/*.h"
)

add_executable("test" ${src})
```

## Peregrine

```bash
module load CMake
module load binutils
```

## All platforms
```
git clone git@github.com:swom/ArchitectureEvolution.git
cd ArchitectureEvolution
git checkout # your branch
mkdir build
cd build
cmake ..
cmake --build .
```

This will generate the binary 'test' in 'build'.


## TlDr
In this project we will study how different architecture evolution algorithms affect the outcome of an evolutionary process
