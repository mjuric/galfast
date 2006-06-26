#!/usr/bin/perl

scalar(@ARGV) == 1 or die("Usage: genselectorbatch.pl <runset_directory>.\n");

$runset = shift @ARGV;

open(BINS, "grep -v '^#.*' $runset/bins.txt|") or die("Cannot open $runset/bins.txt");

$map = 0;
while($_ = <BINS>)
{
	$map++;
	print "# map $map: $_";

	$name="selxy_$map";
	print "qsub -l walltime=04:00:00 -o hydra.astro.princeton.edu:/scr0/mjuric/galaxy/workspace/$runset/outputs/$name.output -N $name ".
		" -v RUNSET=$runset,CMD=selxy,MAP=$map driver_ng\n";

	$name="selcyl_$map";
	print "qsub -l walltime=04:00:00 -o hydra.astro.princeton.edu:/scr0/mjuric/galaxy/workspace/$runset/outputs/$name.output -N $name ".
		" -v RUNSET=$runset,CMD=selcyl,MAP=$map driver_ng\n";
}
