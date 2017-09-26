#!/usr/bin/perl -w
##########################################################################################################################
#	File Name: getPharPaths.pl
#	This program predicts the associations between drugs and diseases and gets
#   pharmacological pathways bridging drugs and diseases using greedy search algorithm
# 	Wroted By Lu Lu and Hua Yu
##########################################################################################################################

use warnings;
use strict;
use Getopt::Long;

my ($inDir,$netFile,$outDir,$help);

GetOptions(
	"inDir=s" => \$inDir,
	"netFile=s" => \$netFile,
	"outDir=s" => \$outDir,
	"help!" => \$help,
);

my %netHash;
open(NET,"<$netFile") or die "$!\n";
while(<NET>){
	my $line = $_;
	chomp $line;
	my @fieldValues = split /\t/,$line;
	$netHash{$fieldValues[0]}{$fieldValues[1]} = $fieldValues[2];
	$netHash{$fieldValues[1]}{$fieldValues[0]} = $fieldValues[2];
}

my @resFiles = `find $inDir -name "prim.simp.result*txt"`;


my (%flowHash, %pathHash); 
foreach my $resfile (@resFiles){
	chomp $resfile;
	$resfile =~ /prim\.simp\.result\.(DB.*)\.\d+\.txt/;
	my $drugid = $1;
	open(RES,"<$resfile") or die "$!\n";
	while(<RES>){
		my $line = $_;
		chomp $line;
		if($line =~ /^V_/){
			my @fieldValues = split /\s+/,$line;
			if($fieldValues[1] > 0){
				print "$fieldValues[1]\n";
				$fieldValues[0] =~ s/V\_//;
				my ($from,$to) = split /_/,$fieldValues[0];
				if($from =~ /^DB/){
					$flowHash{$drugid}{$from}{$to} = $fieldValues[1];
				}elsif($to eq "T"){
					$flowHash{$drugid}{$from}{$to} = $fieldValues[1];
				}else{
					$flowHash{$drugid}{$from}{$to} = $fieldValues[1];
					$flowHash{$drugid}{$to}{$from} = $fieldValues[1];
				}
			}
		}
	}
}

foreach my $drugid (keys %flowHash){
	my @firstgenes = keys %{$flowHash{$drugid}{$drugid}};
	foreach my $firstgene (@firstgenes){
		if($firstgene ne "S"){
			my @secondgenes = keys %{$flowHash{$drugid}{$firstgene}};
			foreach my $secondgene (@secondgenes){
				if($secondgene ne $firstgene){
					my @path = ();
					push @path, $firstgene;
					RankingCandidatesAndGetPharPaths($drugid,$firstgene,$secondgene,\@path,\%flowHash,\%pathHash);
				}
			}
		}
	}
}

sub RankingCandidatesAndGetPharPaths{
	my ($drugid,$from,$to,$path,$flowHash,$pathHash) = @_;
	my %pathelehash = map {$_ => 1} @{$path};
	if(!exists $pathelehash{$to}){
		if($to eq "T"){
			push @{$pathHash -> {$drugid}{$from}{"PATH"}},join("--",@{$path});
		}else{
			push @{$path},$to;
			$from = $to;
			my @tos = keys %{$flowHash -> {$drugid}{$from}};
			foreach my $to (@tos){
				RankingCandidatesAndGetPharPaths($drugid,$from,$to,$path,$flowHash,$pathHash);
			}
		}
	}
}

open(PRE,">$outDir/predicted_associations.txt") or die "$!\n";
open(PATH,">$outDir/predicted_pharpaths.txt") or die "$!\n";
foreach my $drugid (keys %pathHash){
	my @diseaseids = keys %{$pathHash{$drugid}};
	foreach my $diseaseid (@diseaseids){
		my $assoconf = 0;
		my @pharpaths = @{$pathHash{$drugid}{$diseaseid}{"PATH"}};
		foreach my $pharpath (@pharpaths){
			my @path = split /--/, $pharpath;
			delete $path[$#path];
			if(scalar(@path) >= 2){
				my $pathconf = 0;
				my $assopathconf = 0;
				for(my $iteri = 0; $iteri <= ($#path-1); $iteri++){
					for(my $iterj = ($iteri+1); $iterj <= $#path; $iterj++){
						my $firstgene = $path[$iteri];
						my $secondgene = $path[$iterj];
						$pathconf = $pathconf + ($netHash{$firstgene}{$secondgene})*($flowHash{$drugid}{$firstgene}{$secondgene});
						$assopathconf = $assopathconf + ($netHash{$firstgene}{$secondgene});
					}
				}
				$pathconf = $pathconf/(scalar(@path)-1);
				$assopathconf = $assopathconf/(scalar(@path)-1);
				# $pathconf = $pathconf**($alpha*scalar(@path));
				print PATH "$drugid\t$diseaseid\t".join(",",@path)."\t$pathconf\n";
				# $assopathconf = $assopathconf**($alpha*scalar(@path));
				$assoconf = $assoconf + $assopathconf;
			}
		}
		$assoconf = $assoconf * $flowHash{$drugid}{$diseaseid}{"T"};
		print PRE "$drugid\t$diseaseid\t$assoconf\n";
	}
}

# perl getPharPaths.pl --inDir /public/home/hyu/PNet/Result/12 --netFile /public/home/hyu/PNet/Data/HumanNet.txt --outDir /public/home/hyu/PNet/Result/12