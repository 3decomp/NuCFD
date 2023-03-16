#!/bin/sh
for f in $(find . -type f -name "*.f90")
do
	flint lint $f -r utils/fortran_rc_default.yaml
done
