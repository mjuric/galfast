#!/bin/tcsh

#
# Duplicate runset directory structure from <from> to <to>
# but link the large parts
#

if($#argv != 2) then
	echo "$0 <from> <to>"
	exit -1
endif

set from = $argv[1]
set to = $argv[2]

if(! -d $from) then
	echo "Directory $from does not exist"
	exit -1
endif

if(-d $to) then
	echo "Directory $to exists"
	exit -1
endif

mkdir $to

cp $from/bins.txt $to
cp $from/conf.sdss $to
cp $from/fits.txt $to
cp $from/run_geometry.txt $to
cp -r $from/vars $to
ln -s ../$from/dmm $to
ln -s ../$from/maps $to
ln -s ../$from/merged.bin.gz $to
ln -s ../$from/volumes $to
mkdir $to/cache
mkdir $to/outputs
mkdir $to/txt
