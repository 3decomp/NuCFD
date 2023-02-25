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

# Testing
enable_testing()

add_executable(tridsolver_suite tests/tridsolver.f90)
add_test(NAME test_tridsolver COMMAND tridsolver_suite)
