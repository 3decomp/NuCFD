## CMakeLists.txt
#
## Description
#
## LICENSE
#
# SPDX-License-Identifier: BSD-3-Clause
#

## Basic configuration
cmake_minimum_required(VERSION 3.16)
project(NuCFD)
enable_language(Fortran)

include(GNUInstallDirs)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/${CMAKE_INSTALL_INCLUDEDIR})

## Set compiler flags using interface libraries
set(gcc_like_fc "$<COMPILE_LANG_AND_ID:Fortran,GNU>")
set(icc_like_fc "$<COMPILE_LANG_AND_ID:Fortran,Intel>")
add_library(nucfd_compiler_flags INTERFACE)
target_compile_options(nucfd_compiler_flags INTERFACE
  "$<${gcc_like_fc}:-std=f2008;-Wall;-Wextra;-Wpedantic;-Warray-bounds;-Wimplicit-interface;-Wimplicit-procedure;-fimplicit-none>"
  "$<${icc_like_fc}:-std08;-warn;all>")
add_library(nucfd_debug_flags INTERFACE)
target_compile_options(nucfd_debug_flags INTERFACE
  "$<${gcc_like_fc}:-g;-Og;-ffpe-trap=invalid,zero;-fcheck=bounds;-fbacktrace>"
  "$<${icc_like_fc}:-g;-O0;-check;all;-fpe0;-traceback;-init:snan,arrays>")
add_library(nucfd_strict_flags INTERFACE)
target_compile_options(nucfd_strict_flags INTERFACE
  "$<${gcc_like_fc}:-Werror>")

## Build nucfd library
add_library(nucfd
  src/nucfd_trid_solver_mod.f90
  src/nucfd_types_mod.f90
  src/nucfd_coeffs_mod.f90)
target_link_libraries(nucfd nucfd_compiler_flags)
if (CMAKE_BUILD_TYPE MATCHES "Debug")
  target_link_libraries(nucfd nucfd_debug_flags)
endif()

## Testing
enable_testing()

add_library(nucfd_test_flags INTERFACE)
target_link_libraries(nucfd_test_flags INTERFACE
  nucfd_compiler_flags
  nucfd_debug_flags
  nucfd_strict_flags)

add_library(nucfd_test_framework tests/nucfd_tests_mod.f90)
target_link_libraries(nucfd_test_framework nucfd_test_flags)

function(define_test suite test_name)
  add_executable(${test_name} tests/${suite}/${test_name}.f90)
  target_link_libraries(${test_name} nucfd)
  target_link_libraries(${test_name} nucfd_test_framework)
  target_link_libraries(${test_name} nucfd_test_flags)
  add_test(NAME ${suite}:${test_name} COMMAND ${test_name})
endfunction()

function(add_test_libs suite test_name lib)
  target_link_libraries(${test_name} ${lib})
endfunction()

# Tidiagonal solver test suite
add_library(trid_test_utils tests/tridsolver/tridsol_test_utils_mod.f90)
target_link_libraries(trid_test_utils nucfd_test_flags)
define_test(tridsolver system_11_symm nucfd_test_framework)
define_test(tridsolver system_11_anti-symm)
add_test_libs(tridsolver system_11_symm trid_test_utils)
add_test_libs(tridsolver system_11_anti-symm trid_test_utils)
define_test(tridsolver system_00_symm nucfd_test_framework)
define_test(tridsolver system_00_anti-symm)
add_test_libs(tridsolver system_00_symm trid_test_utils)
add_test_libs(tridsolver system_00_anti-symm trid_test_utils)

# Finite difference test suite
define_test(fd-schemes verify_coeffs)

## Examples
add_executable(compute_derivatives examples/compute-derivatives.f90)
target_link_libraries(compute_derivatives nucfd)

## Documentation
add_custom_target(doc ford
  ${CMAKE_SOURCE_DIR}/NuCFD.md)
