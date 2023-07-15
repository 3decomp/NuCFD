## tests/CMakeLists.txt
#
## Description
#
# CMakeLists file controlling the tests for the NuCFD project.
#
## LICENSE
#
# SPDX-License-Identifier: BSD-3-Clause
#

add_library(nucfd_test_flags INTERFACE)
target_link_libraries(nucfd_test_flags INTERFACE
  nucfd_compiler_flags
  nucfd_debug_flags
  nucfd_strict_flags)

add_library(nucfd_test_framework OBJECT
  nucfd_tests_mod.f90)
target_link_libraries(nucfd_test_framework nucfd_test_flags)

function(define_test suite test_name)
  add_executable(${test_name} ${test_name}.f90)
  target_link_libraries(${test_name} nucfd)
  target_link_libraries(${test_name} $<TARGET_OBJECTS:nucfd_test_framework>)
  target_link_libraries(${test_name} nucfd_test_flags)
  add_test(NAME ${suite}:${test_name} COMMAND ${test_name})
endfunction()

function(add_test_libs suite test_name lib)
  target_link_libraries(${test_name} ${lib})
endfunction()

# Tridiagonal solver test suite
add_subdirectory(tridsolver)

# Finite difference test suite
add_subdirectory(fd-schemes)
