cmake_minimum_required(VERSION 3.20)

project(gpf VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)

set(GPF_IS_TOP_LEVEL OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(GPF_IS_TOP_LEVEL ON)
endif()

if(GPF_IS_TOP_LEVEL)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

add_library(gpf_mesh INTERFACE)
add_library(gpf::mesh ALIAS gpf_mesh)
set_target_properties(gpf_mesh PROPERTIES EXPORT_NAME mesh)

target_include_directories(
  gpf_mesh
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(gpf_mesh INTERFACE cxx_std_23)

option(GPF_MESH_BUILD_TESTS "Build mesh tests" ${GPF_IS_TOP_LEVEL})
option(GPF_MESH_INSTALL "Install mesh targets" ${GPF_IS_TOP_LEVEL})
if(GPF_MESH_BUILD_TESTS)
  enable_testing()
  add_executable(
    gpf_mesh_tests
    tests/main.cpp
    tests/test_surface_mesh.cpp
    tests/test_manifold_mesh.cpp
    tests/test_property.cpp
  )
  target_link_libraries(gpf_mesh_tests PRIVATE gpf::mesh)
  add_test(NAME gpf_mesh_tests COMMAND gpf_mesh_tests)
endif()

if(GPF_MESH_INSTALL)
  include(CMakePackageConfigHelpers)

  install(TARGETS gpf_mesh EXPORT gpfTargets INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

  set(gpf_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

  install(
    EXPORT gpfTargets
    NAMESPACE gpf::
    FILE gpfTargets.cmake
    DESTINATION "${gpf_INSTALL_CMAKEDIR}"
  )

  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/gpfConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/gpfConfig.cmake"
    INSTALL_DESTINATION "${gpf_INSTALL_CMAKEDIR}"
  )

  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/gpfConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
  )

  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/gpfConfig.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/gpfConfigVersion.cmake"
    DESTINATION "${gpf_INSTALL_CMAKEDIR}"
  )
endif()
