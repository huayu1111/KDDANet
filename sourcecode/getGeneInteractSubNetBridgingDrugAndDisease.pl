#!/usr/bin/perl -w
##########################################################################################################################
#	File Name: getGeneInteractSubNetBridgingDrugAndDisease.pl
#	This program gets the gene interaction subnetwork bridging drugs and diseases
# 	Wroted By Lu Lu and Hua Yu
##########################################################################################################################

use warnings;
use strict;
use Getopt::Long;

my ($infile,$outDir,$help);

GetOptions(
	"infile=s" => \$infile,
	"outDir=s" => \$outDir,
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
	if($drugid eq "DB03147" && $dis eq "114480"){
		print join("\t",@firstgenes)."\n";
	}
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
        open(OUT,">$outDir/$drugid\_$dis\_subnet.txt") or die "$!\n";
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

# sub writeNetwork{
	# my($file,$nodes,$edges) = @_;
	# open(NET,">$file.sif")|| die "can't not open the file : $!";
	# open(NODE,">$file.NODE")|| die "can't not open the file : $!";
	# open(EDGE,">$file.EDGE")|| die "can't not open the file : $!";	
	
	# my $threshold = 0.00001;
	# foreach my $from(keys %{$edges}){
		# foreach my $to(keys %{$edges->{$from}}){
			# my $type = $edges->{$from}{$to}{'TYPE'};

			# my $flow = abs($edges->{$from}{$to}{'FLOW'});
			# print NET "$from gg $to\n" if($flow > $threshold);		##Jialiang Huang 2010-11-28 '--|->' changed to 'gg'
		# }
	# }
	
	# my %nodeInNet = ();
	# foreach my $from(keys %{$edges}){
		# foreach my $to(keys %{$edges->{$from}}){
			# $nodeInNet{$from} = 1;
			# $nodeInNet{$to} = 1;
		# }
	# }
	# print NODE "nodeID\tnodeFlow\n";
	# foreach my $gene(sort {$a<=>$b}keys %nodeInNet){
		# if(defined $nodes->{$gene}{'FLOW'}){
			# print NODE "$gene\t$nodes->{$gene}{'FLOW'}\n" if($nodes->{$gene}{'FLOW'}>$threshold)
		# }		
	# }
	
	# print EDGE "edgeID\tedgeType\tedgeFlow\n";
	# foreach my $from(keys %{$edges}){
		# foreach my $to(keys %{$edges->{$from}}){
			# my $type = $edges->{$from}{$to}{'TYPE'};
			# my $flow = abs($edges->{$from}{$to}{'FLOW'});
			# print EDGE "$from (gg) $to\t$type\t$flow\n" if($flow > $threshold);	##Jialiang Huang 2010-11-28 '--|->' changed to 'gg'
		# }
	# }
	
	# close(NET);
	# close(NODE);
	# close(EDGE);	
# }

# perl getGeneInteractSubNetBridgingDrugAndDisease.pl --infile /public/ZhangJinLab/project_pharpath/PNet/Result/SDrTDi/10/opti_edge_flow.txt --outDir /public/ZhangJinLab/project_pharpath/PNet/Result/SDrTDi/10
