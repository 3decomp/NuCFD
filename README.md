`NuCFD`, for Non-uniform Compact Finite Differences, is a library implementing compact finite
difference schemes on non-uniform grids.

## Building

The `NuCFD` build system is generated by `cmake`, a default configuration can be created by running
```
cmake -B build
```
from the root directory.
Once the build system has been generated running
```
make -C build
```
will build the library.
The build can be configured using the `ccmake` tool, however currently only `gfortran` is supported.

Note that this project uses Fortran2008 submodules, support for these in CMake requires CMake v3.25.2
or above.
If you already have an older CMake it can be relatively easily updated for the local user by cloning
the repository and checking out tag `v3.25.2` (or later), CMake iteslf can then be built using a
fairly standard process:
```
cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/install
# Any configuration required
make -C build && make -C build install
export PATH=/path/to/install/bin:${PATH}
```
you should then be able to use your new CMake to configure the `NuCFD` build.

## Testing

After building the library it can be tested using `ctest`.
From the root directory run
```
make -C build test
```
to launch the tests.
If any test fail more detailed output can be shown by running `ctest` verbosely
```
cd build
ctest --verbose
```
which will print any output from the tests, add `--rerun-failed` to only repeat failing tests.

## Documentation

Documentation can be generated using the FORD tool by running
```
make -C build/ doc
```
