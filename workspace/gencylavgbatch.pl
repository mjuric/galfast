#!/usr/bin/perl

scalar(@ARGV) == 1 or die("Usage: genselectorbatch.pl <runset_directory>.\n");

$runset = shift @ARGV;
$plx = "3.2 13.30 -11.50 5.40 -0.65";

open(BINS, "$runset/bins.txt") or die("Cannot open $runset/bins.txt");

while($_ = <BINS>)
{
	next if /^#.*/;
	($ri0, $ri1, $dx0, $ndx) = split /\s+/;
	$ri0 = sprintf("%5.3f", $ri0);
	$ri1 = sprintf("%5.3f", $ri1);

	print "./median3d.pl mean pure $runset/txt/cyl.$ri0-$ri1.txt > $runset/txt/avgcyl.$ri0-$ri1.mean.txt\n";
	print "./median3d.pl mean cleaned $runset/txt/cyl.$ri0-$ri1.txt > $runset/txt/avgcyl.$ri0-$ri1.mean.cleaned.txt\n";
}
