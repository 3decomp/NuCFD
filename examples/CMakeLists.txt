## examples/CMakeLists.txt
#
## Description
#
# CMakeLists file for building the NuCFD project examples.
#
## LICENSE
#
# SPDX-License-Identifier: BSD-3-Clause
#

add_executable(compute_derivatives compute-derivatives.f90)
target_link_libraries(compute_derivatives nucfd)
