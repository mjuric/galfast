#!/usr/bin/perl

scalar(@ARGV) == 1 or die("Usage: genbatch.pl <runset_directory>.\n");

$runset = shift @ARGV;

open(GEOM, "/scr0/mjuric/galaxy/workspace/$runset/run_geometry.txt") or die("Cannot open run_geometry file");
while(<GEOM>)
{
	($nada, $run) = split / +/;
	$runs{$run} = 1;
}

open(RUNS, "/scr0/mjuric/galaxy/workspace/$runset/runs.txt") or die("Cannot open runs.txt file");
foreach $run (<RUNS>)
{
	$run = $run + 0;
	next if($run == 0);

	$runs{$run} or die("Run $run is not in run_geometry file");

	print "qsub -l walltime=01:00:00 -o hydra.astro.princeton.edu:/scr0/mjuric/galaxy/workspace/$runset/outputs/ray$run.output -N ray$run -v RUNSET=$runset,RUN=$run driver\n";
}
