## src/coeffs/CMakeLists.txt
#
## Description
#
# CMakeLists file for building the coefficients module of NuCFD.
#
## LICENSE
#
# SPDX-License-Identifier: BSD-3-Clause
#

add_library(nucfd_coeffs OBJECT
  nucfd_coeffs_mod.f90
  nucfd_coeffs_a.f90
  nucfd_coeffs_b.f90
  nucfd_coeffs_c.f90
  nucfd_coeffs_d.f90
  nucfd_coeffs_e.f90
)
target_link_libraries(nucfd_coeffs PRIVATE $<TARGET_OBJECTS:nucfd_types>)
target_link_libraries(nucfd_coeffs PRIVATE nucfd_compiler_flags)
