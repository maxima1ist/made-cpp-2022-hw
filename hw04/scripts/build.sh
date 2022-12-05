#/bin/bash

PARALLEL="OFF"

print_usage() {
  printf "Usage: bash build.sh [-p]\n"
}

while getopts "p" flag; do
  case "${flag}" in
    p) PARALLEL="ON" ;;
    *) print_usage
       exit 1 ;;
  esac
done

rm -rf build
mkdir build
cd build
cmake .. -DPARALLEL=${PARALLEL}
make
