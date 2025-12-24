cmake_minimum_required(VERSION 3.20)

project(mesh1_cpp LANGUAGES CXX)

add_library(mesh1 INTERFACE)
target_include_directories(mesh1 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(mesh1 INTERFACE cxx_std_23)

option(MESH1_BUILD_TESTS "Build mesh1 tests" ON)
if(MESH1_BUILD_TESTS)
  enable_testing()
  add_executable(
    mesh1_tests
    tests/main.cpp
    tests/test_surface_mesh.cpp
    tests/test_manifold_mesh.cpp
    tests/test_property.cpp
  )
  target_link_libraries(mesh1_tests PRIVATE mesh1)
  add_test(NAME mesh1_tests COMMAND mesh1_tests)
endif()
