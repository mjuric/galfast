#!/usr/bin/perl

sub usage { die("./median3d.pl <mean|median> <pure|cleaned> cyl.1.00-1.10.txt > rz.txt"); }

$filter_mon = $ENV{'FILTER_MON'};

scalar(@ARGV == 3) or usage();
($method, $type, $file) = @ARGV;
$method eq "mean" or $method eq "median"  or usage();
$type   eq "pure" or $type   eq "cleaned" or usage();
open(IN, $file) or die("Cannot open $file");

print "# input file          : $file\n";
print "# phi averaging method: $method\n";

sub trim(\$)
{
	my $string = shift;
	for ($$string) {
		s/^\s+//;
		s/\s+$//;
	}
}

sub filter
{
	my $PI = 3.14159;

	my ($r, $rphi, $z) = @_;
	my $phi = ($r != 0 ? $rphi/$r : 0) - $PI;
	my $x = $r*cos($phi);
	my $y = $r*sin($phi);
	my $D = sqrt(($x-8000)*($x-8000) + $y*$y + $z*$z);

	# Filter out the top half where the Virgo overdensity is
	my $theta = $PI/6;
	my $plane1 = sin($theta)*($x - 8000) - cos($theta)*$y;
	my $plane2 = $x + $y;
	if($D > 2500 && $plane1 < 0 && $plane2 > 0) { return 0; }

	# Monoceros stream
	# This used to be the filter for astro-ph version of photom. paralax
	if($filter_mon) {
		##if(13000 < $r && $r < 19000 && 0 < $z && $z < 5000) { return 0; }
		##if(16000 < $r && $r < 19000 && 0 < $z && $z < 7000) { return 0; }
		if(14000 < $r && $r < 24000 && 0 < $z && $z < 7000) { return 0; }
		if(16000 < $r && $r < 24000 && 0 < $z && $z < 10000) { return 0; }
	}

	# Close by (just for testing, the fitter should apply this cut)
	#if(abs($z) < 1000) { return 0; }

	# General cleanup
	$phi2 = atan2($z, $r - 8000)*180/$PI;

	if(-30 < $phi2 && $phi2 < 15) { return 0; }
	if(153 < $phi2 && $phi2 < 180) { return 0; }
	if(-180 < $phi2 && $phi2 < -150) { return 0; }

	return 1;
}

while($_ = <IN>)
{
	trim $_;
	/^#.*$/ and next;	# comments

	($x, $y, $z, $obs, $ovol, $n, $v) = split /\s+/;
	$n == 0 and next;

	#print "($x, $y, $z, $obs, $ovol, $n, $v)\n";
	if($type eq "cleaned") { filter($x, $y, $z) or next; }

	# BUGFIX/HACK: There are some bins that got rounding screwed up
	$x = sprintf("%.0f", $x);
	$y = sprintf("%.0f", $y);
	$z = sprintf("%.0f", $z);

	if($v == 0) { print STDERR "$n $v\n"; }
	push @{$den{$x}{$z}}, {'den' => $n/$v, 'N' => $n, 'V' => $v};
}

foreach $x (sort { $a <=> $b } keys %den)
{
	foreach $y (sort { $a <=> $b } keys %{$den{$x}})
	{
		# Find median
		@a = sort { $a->{'den'} <=> $b->{'den'} } @{$den{$x}{$y}};

		if($method eq "mean") {

			# Mean calculation
			$h{'N'} = $h{'V'} = 0;
			foreach $i (@a) {
				$h{'N'} += $i->{'N'};
				$h{'V'} += $i->{'V'};
			}
			$h{'den'} = $h{'N'}/$h{'V'};

		} elsif($method eq "median") {

			# Median calculation
			$l = scalar(@a);
			$i = int($l / 2);
			%h = %{$a[$i]};
			if($l % 2 == 0) {
				%h2 = %{$a[$i+1]};
				$h{'den'} = ($h{'den'} + $h2{'den'}) / 2.;
				$h{'N'} = ($h{'N'} + $h2{'N'}) / 2.;
				$h{'V'} = ($h{'V'} + $h2{'V'}) / 2.;
			}
		}

		print "$x $y 0 0 0 $h{N} $h{V} $h{den}\n";
	}
}
