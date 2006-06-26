#!/usr/bin/perl

# One-liner which converts ./gen*.pl generated PBS submission scripts to
# a series of commands for local execution

while($_ = <>) { next if /^#.*$/; ($l, $c) = /.*-v (.*) (.*)$/; $l =~ s/,/ /g; print "env $l ./$c\n"}
