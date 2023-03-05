## src/CMakeLists.txt
#
## Description
#
# CMakeLists file for building the NuCFD library.
#
## LICENSE
#
# SPDX-License-Identifier: BSD-3-Clause
#

add_library(nucfd_types OBJECT
  nucfd_types_mod.f90
)
target_link_libraries(nucfd_types PRIVATE nucfd_compiler_flags)

add_subdirectory(coeffs)

add_library(nucfd
  nucfd_trid_solver_mod.f90
  nucfd_deriv_mod.f90)
target_link_libraries(nucfd PRIVATE $<TARGET_OBJECTS:nucfd_types>)
target_link_libraries(nucfd PRIVATE $<TARGET_OBJECTS:nucfd_coeffs>)
target_link_libraries(nucfd PUBLIC nucfd_compiler_flags)
