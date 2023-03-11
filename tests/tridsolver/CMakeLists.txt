## tests/tridsolver/CMakeLists.txt
#
## Description
#
# CMakeLists file for the tri-diagonal solver test suite.
#
## LICENSE
#
# SPDX-License-Identifier: BSD-3-Clause
#

add_library(trid_test_utils OBJECT
  tridsol_test_utils_mod.f90)
target_link_libraries(trid_test_utils nucfd_test_flags)

define_test(tridsolver system_11_symm nucfd_test_framework)
define_test(tridsolver system_11_anti-symm)
add_test_libs(tridsolver system_11_symm $<TARGET_OBJECTS:trid_test_utils>)
add_test_libs(tridsolver system_11_anti-symm $<TARGET_OBJECTS:trid_test_utils>)

define_test(tridsolver system_00_symm nucfd_test_framework)
define_test(tridsolver system_00_anti-symm)
add_test_libs(tridsolver system_00_symm $<TARGET_OBJECTS:trid_test_utils>)
add_test_libs(tridsolver system_00_anti-symm $<TARGET_OBJECTS:trid_test_utils>)
