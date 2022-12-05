#!/bin/bash

cd build/tests/

sum=0
for i in {1..5}
do
  start=$(date +%s%3N)
  ctest
  end=$(date +%s%3N)
  sum=$(($sum + $end - $start))
done

echo "average time: $(($sum / 5)) milliseconds"
