#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($indir,$outdir,$help);

GetOptions(
    "indir=s" => \$indir,
    "outputdir=s" => \$outputdir,
	"help!" => \$help,
);

my @subnetfiles = `find $indir -name "*subnet.txt"`;
foreach my $subnetfile (@subnetfiles){
	chomp $subnetfile;
    $subnetfile =~ /\/(DB.*subnet\.txt)/;
    my $outfile = $1;
    $outfile =~ s/subnet/module/;
	system("mcl $subnetfile --abc -o $outdir/$outfile")
}
}