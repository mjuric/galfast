#!/usr/bin/perl

scalar(@ARGV) == 1 or die("Usage: genselectorbatch.pl <runset_directory>.\n");

$runset = shift @ARGV;
$plx = "3.2 13.30 -11.50 5.40 -0.65";

open(BINS, "$runset/bins.txt") or die("Cannot open $runset/bins.txt");

while($_ = <BINS>)
{
	next if /^#.*/;
	($ri0, $ri1, $dx0, $ndx, $r0, $r1) = split /\s+/;
	$ri0 = sprintf("%5.3f", $ri0);
	$ri1 = sprintf("%5.3f", $ri1);
	$r0 = sprintf("%5.3f", $r0);
	$r1 = sprintf("%5.3f", $r1);

	print "./median3d.pl median pure $runset/txt/cyl.$r0-$r1.$ri0-$ri1.txt > $runset/txt/avgcyl.$ri0-$ri1.median.txt\n";
	print "./median3d.pl median cleaned $runset/txt/cyl.$r0-$r1.$ri0-$ri1.txt > $runset/txt/avgcyl.$ri0-$ri1.median.cleaned.txt\n";
	print "./median3d.pl mean pure $runset/txt/cyl.$r0-$r1.$ri0-$ri1.txt > $runset/txt/avgcyl.$ri0-$ri1.mean.txt\n";
	print "./median3d.pl mean cleaned $runset/txt/cyl.$r0-$r1.$ri0-$ri1.txt > $runset/txt/avgcyl.$ri0-$ri1.mean.cleaned.txt\n";
	print "./median3d_ng.pl $runset/txt/cyl.$r0-$r1.$ri0-$ri1.txt > $runset/txt/cyl.$r0-$r1.$ri0-$ri1.cleaned.txt\n";
}
