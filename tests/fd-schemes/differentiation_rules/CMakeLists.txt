## tests/fd-schemes/differentiation_rules/CMakeLists.txt
#
## Description
#
# CMakeLists file for the finite difference schemes differentiation rules test sub-suite.
#
## LICENSE
#
# SPDX-License-Identifier: BSD-3-Clause
#

add_library(diff_rules_utils OBJECT diff_rules_utils.f90)

define_test(fd-schemes constant_function)
add_test_libs(fd-schemes constant_function diff_rules_utils)
define_test(fd-schemes linear_rising_function)
add_test_libs(fd-schemes linear_rising_function diff_rules_utils)
define_test(fd-schemes linear_falling_function)
add_test_libs(fd-schemes linear_falling_function diff_rules_utils)
