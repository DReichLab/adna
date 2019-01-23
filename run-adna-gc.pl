#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (
	P=>'broad-hominid-dev', # project name
	R=>'results',           # write results to this directory
	Z=>'us-central1-b',     # zone
	M=>'n1-standard-4',     # machine type
	I=>'lh3dev',            # development image
	S=>'200GB'              # startup disk size
);
getopts('p:P:B:m:', \%opts);
die(qq(Usage: run-adna-gc.pl [options] <gs://unmapped.bam> <gs://barcodes.txt> <gs://results.dir>
Options:
  -p STR      prefix [auto]
  -P STR      gcloud project name [$opts{P}]
  -m STR      machine name [auto]
)) if @ARGV < 3;

my $n_threads = $opts{M} =~ /-(\d+)$/? $1 : 4;

my $prefix;
if (defined $opts{p}) {
	$prefix = $opts{p};
} else {
	$prefix = $ARGV[0] =~ /\/([^\s\/]+)$/? $1 : $ARGV[0];
	$prefix =~ s/\.bam$//;
	$prefix =~ s/\.unmapped$//;
}
my $machine;
if (defined $opts{m}) {
	$machine = $opts{m};
} else {
	$machine = "$prefix-ins";
	$machine =~ tr/._:()[]{}|*@/-/;
	$machine = lc($machine);
}

system(qq/gcloud compute instances create $machine --zone $opts{Z} --boot-disk-type pd-standard --boot-disk-size 200GB --disk name=$opts{I},mode=ro --machine-type $opts{M} --scopes storage-rw/);
open(FH, qq/| gcloud compute ssh $machine --zone $opts{Z} --command "cat > run.sh"/) || die;
print FH qq(mkdir -p lh3dev && sudo mount -o discard,defaults,ro /dev/sdb lh3dev
sudo apt-get update && sudo apt-get -q -y install perl
gsutil -u $opts{P} cp $ARGV[1] $prefix.bc
lh3dev/adna.kit/run-adna -b $prefix.bc -t $n_threads -p $prefix lh3dev/bwadb/hs37d5.fa "gsutil cp $ARGV[0] -" | sh
gsutil -u $opts{P} cp $prefix.* $ARGV[2]
);
close(FH);
system(qq/gcloud compute ssh $machine --zone $opts{Z} --command "sh run.sh"/);
system(qq/gcloud compute -q instances delete $machine --zone $opts{Z}/);
