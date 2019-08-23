#!/usr/bin/perl -w
use warnings;
use Getopt::Long;

my ($indir,$outdir,$help);

GetOptions(
	"indir=s" => \$indir,
	"outdir=s" => \$outdir,
	"help!" => \$help,
);

my (%flowEdges,%newEdges,%flowNodes);
my @resFiles = `find $indir -name "prim.simp.result*txt"`; 
foreach my $resfile (@resFiles){
	chomp $resfile;
	$resfile =~ /prim\.simp\.result\.(DB.*)\.\d+\.txt/;
	my $drugid = $1;
    getEdgeFlow($resfile,$drugid,\%flowEdges,\%newEdges);
}

open(OUT,">$outdir/opti_edge_flow.txt") or die "$!\n";
foreach my $drugid (keys %newEdges){
    foreach my $from (keys %{$newEdges{$drugid}}){
        foreach my $to (keys %{$newEdges{$drugid}{$from}}){
            print OUT "$drugid\t$from\t$to\t$newEdges{$drugid}{$from}{$to}\n";
        }
    }
}
close OUT;

getNodeFlow(\%newEdges,\%flowNodes);

open(OUT,">$outdir/opti_node_flow.txt") or die "$!\n";
foreach my $drugid (keys %flowNodes){
    foreach my $node (keys %{$flowNodes{$drugid}}){
        print OUT "$drugid\t$node\t$flowNodes{$drugid}{$node}\n";
    }
}
close OUT;
    
sub getEdgeFlow{
	my ($resfile,$drugid,$flowEdges,$newEdges) = @_;
	
    my (%edges,%newEdges);
    open(RES,"<$resfile") or die "$!\n";
	while(<RES>){
		my $line = $_;
		chomp $line;
		if($line =~ /^V_/){
			my @fieldValues = split /\s+/,$line;
			if($fieldValues[1] != 0){
				$fieldValues[0] =~ s/V\_//;
				my ($from,$to) = split /_/,$fieldValues[0];
                $edges{$drugid}{$from}{$to} = $fieldValues[1];
                $flowEdges -> {$drugid}{"$from\_$to"} = $fieldValues[1];
			}
		}
	}
    close RES;
    
    foreach my $drugid (keys %edges){
        foreach my $from (keys %{$edges{$drugid}}){
            foreach my $to (keys %{$edges{$drugid}{$from}}){
                if($from ne $drugid && $to ne "T"){
                    my $newFrom = $from;
                    my $newTo = $to;
                    my $newFlow = 0;
                    my $flowEdge1 = $from."_".$to;
                    my $flowEdge2 = $to."_".$from;
                    if(exists $flowEdges -> {$drugid}{$flowEdge2}){
                        $newFlow = $flowEdges -> {$drugid}{$flowEdge1} - $flowEdges -> {$drugid}{$flowEdge2};
                        if($newFlow < 0){					
                            $newFlow = -$newFlow;
                            $newFrom = $to;
                            $newTo = $from;    
                        }
                        $newEdges->{$drugid}{$newFrom}{$newTo} = $newFlow;
                    }elsif($flowEdges -> {$drugid}{$flowEdge1} < 0){				
                        $newEdges->{$drugid}{$newFrom}{$newTo} = -$flowEdges -> {$drugid}{$flowEdge1};
                    }else{
                        $newEdges->{$drugid}{$newFrom}{$newTo} = $flowEdges-> {$drugid}{$flowEdge1};
                    }
                }
            }
        }
    }
	
	##	edges from S and to T
	foreach my $drugid (keys %{$flowEdges}){
        foreach my $flowEdge (keys %{$flowEdges{$drugid}}){
            my ($from,$to) = split(/\_/,$flowEdge);
            my $flow = $flowEdges->{$drugid}{$flowEdge};		
            if($from eq $drugid){
                $newEdges->{$drugid}{$drugid}{$to} = $flow;
            }elsif($to eq 'T'){
                $newEdges->{$drugid}{$from}{'T'} = $flow;			
            }
        }
	}
    
}

sub getNodeFlow{
	my ($newEdges,$flowNodes) = @_;
	foreach my $drugid (keys %{$newEdges}){
        foreach my $from (keys %{$newEdges -> {$drugid}}){
            foreach my $to (keys %{$newEdges -> {$drugid}{$from}}){
                $flowNodes->{$drugid}{$from} = 0;
                $flowNodes->{$drugid}{$to} = 0;	
            }
        }
	}	
	
	foreach my $drugid (keys %{$newEdges}){
        foreach my $from (keys %{$newEdges -> {$drugid}}){
            foreach my $to (keys %{$newEdges -> {$drugid}{$from}}){
                print $newEdges -> {$drugid}{$from}{$to} if($to eq "T"); 
                $flowNodes->{$drugid}{$to} += $newEdges -> {$drugid}{$from}{$to};
            }
        }
	}
	
	foreach my $drugid (keys %{$newEdges}){
        foreach my $from (keys %{$newEdges -> {$drugid}}){
            foreach my $to (keys %{$newEdges -> {$drugid}{$from}}){
                print $newEdges -> {$drugid}{$from}{$to} if($from eq $drugid); 
                $flowNodes->{$drugid}{$drugid} += $newEdges -> {$drugid}{$from}{$to} if($from eq $drugid);
            }
        }
	}	 
}