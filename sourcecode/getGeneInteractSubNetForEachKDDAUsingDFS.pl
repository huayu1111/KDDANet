#!/usr/bin/perl -w
##########################################################################################################################
#	File Name: getGeneInteractSubNetForEachKDDAUsingDFS.pl
#	This program gets the gene interaction subnetwork for each KDDA using depth-first searching
# 	Wroted By Lu Lu and Hua Yu
##########################################################################################################################

use warnings;
use strict;
use Getopt::Long;

my ($infile,$outdir,$help);

GetOptions(
	"infile=s" => \$infile,
	"outdir=s" => \$outdir,
	"help!" => \$help,
);

my %edgeFlows;
open(IN,"<$infile") or die "$!\n";
while(<IN>){
    my $line = $_;
    chomp $line;
    my @fieldValues = split /\t/,$line;
    $edgeFlows{$fieldValues[0]}{$fieldValues[2]}{$fieldValues[1]} = $fieldValues[3];
}
close IN;


foreach my $drugid (keys %edgeFlows){
	my @diss =keys %{$edgeFlows{$drugid}{"T"}};
    	my %diss = map {$_ => 1} @diss;
	foreach my $dis (@diss){
        my (%subnet,%markHash);
        my @firstgenes = keys %{$edgeFlows{$drugid}{$dis}};
        foreach my $firstgene (@firstgenes){
            my @secondgenes = keys %{$edgeFlows{$drugid}{$firstgene}};
            foreach my $secondgene (@secondgenes){
		if($secondgene ne $drugid){
			$subnet{$firstgene}{$secondgene} = $edgeFlows{$drugid}{$firstgene}{$secondgene};
			$markHash{$firstgene}{$secondgene} = 1;
			$markHash{$secondgene}{$firstgene} = 1;
			getPaths($drugid,$secondgene,\%subnet,\%edgeFlows,\%markHash,\%diss);
		}
            }
	}
        open(OUT,">$outdir/$drugid\_$dis\_subnet.txt") or die "$!\n";
        foreach my $from (keys %subnet){
            foreach my $to (keys %{$subnet{$from}}){
                print OUT "$from\t$to\t$subnet{$from}{$to}\n";
            }
        }
        close OUT;
	}
}

sub getPaths{
	my ($drugid,$to,$subnet,$edgeFlows,$markHash,$diss) = @_;
	my $from = $to;
	my @tos = keys %{$edgeFlows -> {$drugid}{$from}};
	foreach my $to (@tos){
            next if($to eq $drugid);
            next if(exists $diss -> {$to});
            if(!exists $markHash -> {$from}{$to}){
                $subnet -> {$from}{$to} = $edgeFlows{$drugid}{$from}{$to};
                getPaths($drugid,$to,$subnet,$edgeFlows,$markHash,$diss);
        }
    }
}

