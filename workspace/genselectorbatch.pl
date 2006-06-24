#!/usr/bin/perl

scalar(@ARGV) == 1 or die("Usage: genselectorbatch.pl <runset_directory>.\n");

$runset = shift @ARGV;

open(BINS, "$runset/bins.txt") or die("Cannot open $runset/bins.txt");

while($_ = <BINS>)
{
	next if /^#.*/;
	($ri0, $ri1, $dx0, $ndx) = split /\s+/;

	$name="dcar${ri0}_$ri1";
	print "qsub -l walltime=04:00:00 -o hydra.astro.princeton.edu:/scr0/mjuric/galaxy/workspace/$runset/outputs/$name.output -N $name ".
		" -v RUNSET=$runset,RI0=$ri0,RI1=$ri1,DX0=$dx0,NDX=$ndx driver_sel_xy\n";

	$name="dcyl${ri0}_$ri1";
	print "qsub -l walltime=04:00:00 -o hydra.astro.princeton.edu:/scr0/mjuric/galaxy/workspace/$runset/outputs/$name.output -N $name ".
		" -v RUNSET=$runset,RI0=$ri0,RI1=$ri1,DX0=$dx0,NDX=$ndx driver_sel_cyl\n";
}
