#!/bin/bash

echo -n "Testing, please wait... ";

#../../debug/src/simulate.x foot foot.beam.conf foot.xgpc.txt #> output.log 2>&1
#../../debug/src/simulate.x pdf pdf.conf foot.xgpc.txt model.BahcallSoneira.conf sky.pdf.bin #>> output.log 2>&1
#../../debug/src/simulate.x catalog catalog.conf sky.pdf.bin sky.cat.txt #- | awk '{if (count++%1000==0) print $0;}' > sky.cat.txt 
../../debug/src/simulate.x postprocess postprocess.conf sky.cat.txt sky.obs.txt # >> output.log 2>&1

#cmp sky.obs.txt result/sky.obs.txt || (echo "Error, difference in result.")
#cmp sky.obs.txt result/sky.obs.txt || (echo "Error, difference in result.")
#echo "OK.";

#rm -f output.log sky.cat.txt sky.obs.txt sky.pdf.bin foot.xgpc.txt 
