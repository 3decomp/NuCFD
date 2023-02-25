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

# Set compiler flags using interface library
set(gcc_like_fc "$<COMPILE_LANG_AND_ID:Fortran,GNU>")
add_library(nucfd_compiler_flags INTERFACE)
target_compile_options(nucfd_compiler_flags INTERFACE
  "$<${gcc_like_fc}:-Wall;-Wextra;-Wimplicit-interface;-Wimplicit-procedure>")
add_library(nucfd_strict_flags INTERFACE)
target_compile_options(nucfd_strict_flags INTERFACE
  "$<${gcc_like_fc}:-Werror>")

# Testing
enable_testing()

add_library(nucfd_test_flags INTERFACE)
target_link_libraries(nucfd_test_flags INTERFACE
  nucfd_compiler_flags
  nucfd_strict_flags)

add_executable(tridsolver_suite tests/tridsolver.f90)
target_link_libraries(tridsolver_suite nucfd_test_flags)
add_test(NAME test_tridsolver COMMAND tridsolver_suite)