#######################################################################################################################################
##  Linear Programming Model (LPM) to solve minimum cost flow optimization problem 
##	Wroted by Hua Yu and Lu Lu on 2018-08-09
#######################################################################################################################################

sub lpsolveModeling{
	my ($source,$sourceHash,$targetHash,$sourcemapHash,$targetmapHash,$ddHash,$netFile,$outDir,$gamma,$optimethod) = @_;
    print "$outDir\n";
	##	Read network file
	my %net = ();
	readEdgeFile($netFile,\%net);
	open(ERROR,">>$outDir/error.log") or die "$!\n";
	
	if($optimethod eq "SDrTDi"){
		if(!exists $sourcemapHash -> {$source}){
			print ERROR "Unrecognizable, $source is not included in our database\n";
			return(1)
		}else{
			my @genes = keys %{$sourcemapHash -> {$source}};
			if(scalar(grep {$net{$_}} @genes) == 0){
				print ERROR "The targets of $source are not included in gene network\n";
				return(1)
			}else{
				addNodeByEdge($source,$sourcemapHash,\%net,$optimethod,"Drug");
			}
		}
		
		my $mark = 0;
		foreach my $target (keys %{$targetHash}){
			my @genes = keys %{$targetmapHash -> {$target}};
			
			if(scalar(grep {$net{$_}} @genes) == 0){
				print ERROR "The related genes of $target are not included in gene network\n";
				print ERROR "This disease will be discarded\n";
			}else{
				addNodeByEdge($target,$targetmapHash,\%net,$optimethod,"Disease");
				$mark = 1;
			}
		}
		if($mark == 0){
			return(1);
		}
		
		# get nodes included in network
		my %nodes = ();	
		getNodes(\%net,\%nodes);
		
		# set capacities of edges
		setEdgeCapacity(\%net);
		
		# set attributes of nodes
		setNodesAttr($sourceHash,$targetHash,\%nodes);
		
		# add artificial nodes and edges
		addHypoNodes(\%nodes);
		addHypoEdges(\%net,$source,$ddHash);
		
		##	linear programming formulation
		my @flowsName = ();
		my @coefficients = ();
		if(!-e "$outDir/model.$source.$gamma.lp" || -z "$outDir/model.$source.$gamma.lp"){
			formulation("$outDir/model.$source.$gamma.lp",$gamma,\%nodes,\%net,\@flowsName,\@coefficients,$source);
			# writeNodeAttri($outDir."$source.$gamma.node.type",\%nodes);
			# writeCoefficients($outDir."$source.$gamma.object.in",\@coefficients,\@flowsName);
		}
		
		##	Solve the model by primal simplex algorithm 
		if(!-e "$outDir/prim.simp.result.$source.$gamma.txt" || -z "$outDir/prim.simp.result.$source.$gamma.txt"){
			LPSolver($outDir,$gamma,$source);
			return(0);
		}

	}
    
    if($optimethod eq "SDiTDr"){
		if(!exists $sourceHash -> {$source}){
			print ERROR "Unrecognizable, $source is not included in our database\n";
			return(1);
		}else{
			my @genes = keys %{$sourcemapHash -> {$source}};
			if(scalar(grep {$net{$_}} @genes) == 0){
				print ERROR "The related genes of $source are not included in gene network\n";
                print ERROR "This disease will be discarded\n";
				return(1);
			}else{
				addNodeByEdge($source,$sourcemapHash,\%net,$optimethod,"Disease");
			}
		}
		
		my $mark = 0;
		foreach my $target (keys %{$targetHash}){
			my @genes = keys %{$targetmapHash -> {$target}};
			
			if(scalar(grep {$net{$_}} @genes) == 0){
				print ERROR "The targets of $target are not included in gene network\n";
			}else{
				addNodeByEdge($target,$targetmapHash,\%net,$optimethod,"Drug");
				$mark = 1;
			}
		}
		if($mark == 0){
			return(1);
		}
		
		# get nodes included in integrated network
		my %nodes = ();	
		getNodes(\%net,\%nodes);
		
		# set capacities of edges
		setEdgeCapacity(\%net);
		
		# set attributes of nodes
		setNodesAttr($sourceHash,$targetHash,\%nodes);
		
		# add artificial nodes and edges
		addHypoNodes(\%nodes);
        
		addHypoEdges(\%net,$source,$ddHash);
		
		##	construct linear programming model
		my @flowsName = ();
		my @coefficients = ();
		if(!-e "$outDir/model.$source.$gamma.lp" || -z "$outDir/model.$source.$gamma.lp"){
			formulation("$outDir/model.$source.$gamma.lp",$gamma,\%nodes,\%net,\@flowsName,\@coefficients,$source);
		} 
		# writeNodeAttri($outDir."$source.$gamma.node.type",\%nodes);
		# writeCoefficients("$outDir/$source.$gamma.object.in",\@coefficients,\@flowsName);
		
		##	solve the model by primal simplex algorithm
		if(!-e "$outDir/prim.simp.result.$source.$gamma.txt" || -z "$outDir/prim.simp.result.$source.$gamma.txt"){
			LPSolver($outDir,$gamma,$source);
			return(0);
		}
	}
}

sub addNodeByEdge{
	my($node,$mapHash,$net,$method,$type) = @_;
	print "adding nodes (drugs or diseases) to network by related targets or genes\n";
	my @mappedEles = keys %{$mapHash -> {$node}};
	my @commEles = grep {$net -> {$_}} @mappedEles;
	if(($method eq "SDrTDi" && $type eq "Drug") || ($method eq "SDiTDr" && $type eq "Disease")){
		foreach my $commEle (@commEles){
			$net -> {$node}{$commEle}{'WEIGHT'} = 1/scalar(@commEles);
            $net -> {$node}{$commEle}{'CAPACITY'} = 1/scalar(@commEles);
		}
	}else{
		foreach my $commEle (@commEles){
			$net -> {$commEle}{$node}{'WEIGHT'} = 1/scalar(@commEles);
            $net -> {$commEle}{$node}{'CAPACITY'} = 1/scalar(@commEles);
		}
	}
}

sub getNodes{
	my ($edges,$nodes) = @_;	
	foreach my $from(keys %{$edges}){
		foreach my $to(keys %{$edges->{$from}}){
			$nodes->{$from}{'TYPE'} = '';
			$nodes->{$to}{'TYPE'} = '';
		}
	}
}

sub setEdgeCapacity{
	my($edges) = @_;
	foreach my $from(keys %{$edges}){
        foreach my $to(keys %{$edges->{$from}}){
            $edges->{$from}{$to}{'CAPACITY'} = 1 if(!defined $edges->{$from}{$to}{'CAPACITY'});
        }
	}
}

sub setNodesAttr{
	my($drugs,$diseases,$nodes) = @_;
	
	foreach my $node (keys %{$nodes}){
		if(defined $drugs->{$node}){
			$nodes->{$node}{'TYPE'} = 'DRUG';
		}elsif(defined $diseases->{$node}){
			$nodes->{$node}{'TYPE'} = 'DISEASE';
		}else{
			$nodes->{$node}{'TYPE'} = 'GENE';
		}								
	}
}

sub addHypoNodes{
	my($nodes) = @_;
	$nodes->{'T'}{'TYPE'} = 'T';
}

sub addHypoEdges{
	my ($net,$source,$ddHash) = @_;
    
	##	V - SINK
    foreach my $element (keys %{$$ddHash{$source}}){
        $net->{$element}{'T'}{'WEIGHT'} = 1/scalar(keys %{$$ddHash{$source}});
        $net->{$element}{'T'}{'CAPACITY'} =  1/scalar(keys %{$$ddHash{$source}});
    }
}

sub formulation{
	my($file,$gamma,$nodes,$net,$flowsName,$coefficients,$source) = @_;
	print "LP Modeling...\n";
	print "write file $file...\n";
	
	my @nodesName = ();
	my $count = 0;
	foreach my $node (sort keys %{$nodes}){
		if($node ne $source && $node ne 'T'){
			$nodesName[$count] = $node;
			$count++;
		}
	}
	
	my $count = 0;
	foreach my $from(sort keys %{$net}){
		foreach my $to (sort keys %{$net -> {$from}}){
			$flowsName->[$count] = $from."_".$to;
			$count++;
		}
	}
	
	my $count = 0;
	my %coeffHash = ();
	my %capacityHash = ();
	my $numFlow = @{$flowsName};
	for(my $i=0; $i<$numFlow; $i++){
		my($from,$to) = split(/\_/,$flowsName->[$i]);
		$capacityHash{$flowsName->[$i]} = $net->{$from}{$to}{'CAPACITY'};
		my $weight = $net->{$from}{$to}{'WEIGHT'};
		
		if($from eq $source){
			$coefficients -> [$i] = -log($weight)/log(10) - $gamma;
		}else{
			$coefficients -> [$i] = -log($weight)/log(10);
		}
	}
	for(my $i=0; $i<$numFlow; $i++){
		$coeffHash{$flowsName->[$i]} = $coefficients -> [$i];
	}
	writeLPFile($file,\%coeffHash,\%capacityHash,\@nodesName,$flowsName,$source);
}

sub writeLPFile{
	my($file,$coeffHash,$capacityHash,$nodesName,$flowsName,$source) = @_;
	
	open(OUT,">$file")|| die "can't not open the file : $!";
	print "minimize objective function...\n";
	print OUT "min: ";
	my $first = 1;
	my $numFlow = @{$flowsName};
	for(my $i=0;$i<$numFlow;$i++){
		my $flow = $flowsName->[$i];
		my $variable =  'V_'.$flow;
		my $weight = $coeffHash->{$flow};	
		if($weight >= 0){
			if($first){
				print OUT " $weight $variable ";
				$first = 0;
			}else{
				print OUT "+ $weight $variable ";
			}
		}elsif($weight < 0){
			$weight = abs($weight);
			if($first){
				print OUT "- $weight $variable ";
				$first = 0;
			}else{
				print OUT "- $weight $variable ";
			}			
		}
	}
	print OUT ";\n";

	## Subject to constraits
	my %INs = ();				
	my %OUTs = ();				
	for(my $j=0;$j<$numFlow;$j++){
		my($from,$to) = split(/\_/,$flowsName->[$j]);
		$INs{$to}{$from} = 1;
		$OUTs{$from}{$to} = -1;
	}

	for(my $i=0;$i<@{$nodesName};$i++){
		my $node = $nodesName->[$i];
		my $rhs = 0;
		my $dir = '=';
		my $first = 1;
		foreach my $from (keys %{$INs{$node}}){
			my $flow = $from.'_'.$node;
			my $variable = 'V_'.$flow;
			my $weight = 1;
			if($weight > 0){
				if($first){
					print OUT "$weight $variable ";
					$first = 0;
				}else{
					print OUT "+ $weight $variable ";
				}
			}elsif($weight < 0){
				$weight = abs($weight);
				if($first){
					print OUT "- $weight $variable ";
					$first = 0;
				}else{
					print OUT "- $weight $variable ";
				}			
			}				
		}
		foreach my $to (keys %{$OUTs{$node}}){
			if($node ne $to){
				my $flow = $node.'_'.$to;
				my $variable = 'V_'.$flow;
				my $weight = -1;
				if($weight >= 0){
					if($first){
						print OUT "$weight $variable ";
						$first = 0;
					}else{
						print OUT "+ $weight $variable ";
					}
				}elsif($weight < 0){
					$weight = abs($weight);
					if($first){
						print OUT "- $weight $variable ";
						$first = 0;
					}else{
						print OUT "- $weight $variable ";
					}			
				}
			}				
		}
		print OUT "\t$dir\t$rhs;\n";
	}

	my $rhs = 0;
	my $dir = '=';
	my $first = 1;
	foreach my $from (keys %{$INs{'T'}}){
		my $flow = $from.'_'.'T';
		my $variable = 'V_'.$flow;
		my $weight = 1;
		if($weight > 0){
			if($first){
				print OUT "$weight $variable ";
				$first = 0;
			}else{
				print OUT "+ $weight $variable ";
			}
		}elsif($weight < 0){
			$weight = abs($weight);
			if($first){
				print OUT "- $weight $variable ";
				$first = 0;
			}else{
				print OUT "- $weight $variable ";
			}			
		}				
	}	
	foreach my $to (keys %{$OUTs{$source}}){
		my $flow = "$source"."_"."$to";
		my $variable = 'V_'.$flow;
		my $weight = -1;		
		if($weight > 0){
			if($first){
				print OUT "$weight $variable ";
				$first = 0;
			}else{
				print OUT "+ $weight $variable ";
			}
		}elsif($weight < 0){
			$weight = abs($weight);
			if($first){
				print OUT "- $weight $variable ";
				$first = 0;
			}else{
				print OUT "- $weight $variable ";
			}			
		}			
	}	
	print OUT "\t$dir\t$rhs;\n";
	
	for(my $i=0;$i<$numFlow;$i++){
		my $flow = $flowsName->[$i];
		my $variable =  'V_'.$flow;
		my $capicity = $capacityHash->{$flow};	
		print OUT "0\t<=\t$variable\t<=\t$capicity;\n";		
	}
	close(OUT);
}

sub writeNodeAttri{
	my($file,$nodes) = @_;
	open(OUT,">$file")|| die "can't not open the file : $!";
	
	print OUT "nodeID\tnodeType\n";
	foreach my $node (keys %{$nodes}){
		my $type = $nodes->{$node}{'TYPE'};
		print OUT "$node\t$type\n";
	}
	close(OUT);
}

sub writeCoefficients{
	my ($file,$coefficients,$names) = @_;
	print "write file: $file...\n";
	open(OUT,">$file")|| die "can't not open the file : $!";
	print OUT "ID\tValue\n";
	for(my $i=0;$i<@{$coefficients};$i++){
		print OUT "$names->[$i]\t$coefficients->[$i]\n";
	}
	close(OUT);
}

sub LPSolver{
	my ($outDir,$gamma,$source) = @_;
	if(!-e "$outDir/prim.simp.result.$source.$gamma.txt" || -z "$outDir/prim.simp.result.$source.$gamma.txt"){
		my @commands = ();
		$commands[0] = join(' ', '/public/ZhangJinLab/project_pharpath/PNet/Proc/lpsolver','-lp',"$outDir/model.$source.$gamma.lp",'-prim','>',"$outDir/prim.simp.result.$source.$gamma.txt");
		$commands[1] = join(' ','rm -rf',"$outDir/model.$source.$gamma.lp");
		my $commands = join("\n",@commands)."\n";  
		# open(LSF,">$outDir/lpsolver_".$source."_".$gamma.".lsf") or die "$!\n";
		# print LSF <<END;
# #BSUB -J lpsolver_$source\_$gamma
# #BSUB -n 1
# #BSUB -o $outDir/lpsolver_$source\_$gamma.log
# #BSUB -e $outDir/lpsolver_$source\_$gamma.err
# #BSUB -q normal

# date
# $commands
# date	
# END
		# close LSF;
		
        # my $taskNum =`bjobs | grep RUN | grep normal | grep lpsolver | wc -l`; 
		# while($taskNum > 100){
			# print "The num of task remaining $taskNum\n";
			# sleep 30;
			# print `date`;
			# $taskNum = `bjobs | grep RUN | grep normal | grep lpsolver | wc -l`;
		# }
        # my $bsub = `bsub < $outDir/lpsolver_$source\_$gamma.lsf`;
        # if($bsub ne ""){
            # print "$bsub\n";
        # }
        open(SH,">$outDir/lpsolver_".$source."_".$gamma.".sh") or die "$!\n";
        print SH $commands;
        close SH;
        my $taskNum = `ps -aux | grep lpsolver | wc -l`;
        while($taskNum > 200){
			print "The num of task remaining $taskNum\n";
			sleep 30;
			print `date`;
			$taskNum = `ps -aux | grep lpsolver | wc -l`;
		}
        my $out = system("sh $outDir/lpsolver_$source\_$gamma.sh 1>>$outDir/lpsolver_$source\_$gamma.log 2>>$outDir/lpsolver_$source\_$gamma.err &");
        if($out==0){
            print "The task of $sample_id is successfully submitted\n";
        }
	}
}

1;



