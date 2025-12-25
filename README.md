# gpf (mesh)

Header-only C++23 mesh data structures under the `gpf` namespace.

**CMake target:** `gpf::mesh` (INTERFACE library)

## Requirements

- CMake 3.20+
- A C++23 compiler

## Use in Your Project

### Option A: Add as a subdirectory

```cmake
add_subdirectory(path/to/gpf)
target_link_libraries(your_target PRIVATE gpf::mesh)
```

Tests/install are disabled by default when `gpf` is not the top-level project.

### Option B: Install + `find_package`

Build and install:

```sh
cmake -S . -B build -DGPF_MESH_INSTALL=ON
cmake --build build -j
cmake --install build --prefix <prefix>
```

Consume from another project:

```cmake
find_package(gpf CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE gpf::mesh)
```

If needed, point CMake at the install prefix:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=<prefix>
```

## Headers

Public headers live in `include/` and are included like:

```cpp
#include <gpf/mesh.hpp>
// or specific modules:
// #include <gpf/surface_mesh.hpp>
// #include <gpf/manifold_mesh.hpp>
// #include <gpf/property.hpp>
```

## Build and Run Tests

```sh
cmake -S . -B build -DGPF_MESH_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## CMake Options

- `GPF_MESH_BUILD_TESTS` (default: ON if top-level, otherwise OFF)
- `GPF_MESH_INSTALL` (default: ON if top-level, otherwise OFF)
