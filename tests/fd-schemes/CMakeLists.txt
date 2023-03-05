## tests/fd-schemes/CMakeLists.txt
#
## Description
#
# CMakeLists file for the finite difference schemes test suite.
#
## LICENSE
#
# SPDX-License-Identifier: BSD-3-Clause
#

define_test(fd-schemes verify_coeffs)

add_subdirectory(differentiation_rules)
