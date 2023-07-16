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

project(examples)

add_executable(compute_derivatives compute-derivatives.f90)
target_link_libraries(compute_derivatives nucfd)

add_custom_target(${PROJECT_NAME})
add_dependencies(${PROJECT_NAME} compute_derivatives)
