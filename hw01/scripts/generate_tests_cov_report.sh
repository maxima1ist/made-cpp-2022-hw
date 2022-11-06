#!/bin/bash

cd build
lcov -t "tests/tests" -o coverage.info -c -d tag/
genhtml -o report coverage.info
