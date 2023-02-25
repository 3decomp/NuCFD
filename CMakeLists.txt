## CMakeLists.txt
#
## Description
#
## LICENSE
#
# SPDX-License-Identifier: BSD-3-Clause
#

# Basic configuration
cmake_minimum_required(VERSION 3.16)
project(NuCFD)
enable_language(Fortran)

# Set compiler flags using interface libraries
set(gcc_like_fc "$<COMPILE_LANG_AND_ID:Fortran,GNU>")
add_library(nucfd_compiler_flags INTERFACE)
target_compile_options(nucfd_compiler_flags INTERFACE
  "$<${gcc_like_fc}:-std=f2008;-Wall;-Wextra;-Warray-bounds;-Wimplicit-interface;-Wimplicit-procedure;-fimplicit-none>")
add_library(nucfd_debug_flags INTERFACE)
target_compile_options(nucfd_debug_flags INTERFACE
  "$<${gcc_like_fc}:-ffpe-trap=invalid,zero;-fcheck=bounds>")
add_library(nucfd_strict_flags INTERFACE)
target_compile_options(nucfd_strict_flags INTERFACE
  "$<${gcc_like_fc}:-Werror>")

# Build nucfd library
add_library(nucfd
  src/nucfd_trid_solver_mod.f90
  src/nucfd_coeffs_mod.f90)
target_link_libraries(nucfd nucfd_compiler_flags)

# Testing
enable_testing()

add_library(nucfd_test_flags INTERFACE)
target_link_libraries(nucfd_test_flags INTERFACE
  nucfd_compiler_flags
  nucfd_debug_flags
  nucfd_strict_flags)

add_executable(tridsolver_suite tests/tridsolver.f90)
target_link_libraries(tridsolver_suite nucfd_test_flags)
target_link_libraries(tridsolver_suite nucfd)
add_test(NAME test_tridsolver COMMAND tridsolver_suite)
