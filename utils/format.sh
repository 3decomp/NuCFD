#!/usr/bin/env sh

for f in $(find ${PWD} -type f -name "*.f90")
do
    emacs $f --batch --eval="(delete-trailing-whitespace)" -f save-buffer
    emacs $f --batch --eval="(mark-whole-buffer) (indent-region)"
done
