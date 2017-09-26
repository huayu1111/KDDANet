sub readEdgeFile{
	my ($file,$edges) = @_;
	
	print "Reading the file\t$file\n";
	open(IN,"$file")|| die "can't not open the file : $!";
	while(<IN>){
		chomp();
		my ($from,$to,$weight) = split(/\s+/,$_);
		if(!$weight){
			$edges->{$from}{$to}{'WEIGHT'} = 1;	
			$edges->{$to}{$from}{'WEIGHT'} = 1;	
		}else{
			$edges->{$from}{$to}{'WEIGHT'} = abs($weight);	
			$edges->{$to}{$from}{'WEIGHT'} = abs($weight);
		}
	}
	close(IN);
}

sub readNodeFile{
	my($file,$nodes) = @_;
	print "Reading file\t$file\n";
	open(IN,"$file")|| die "can't not open the file : $!";
	while(<IN>){
		$node = $_;
		$node =~ s/^\s+|\s+//g;
		$nodes->{$node} = 1;
	}		
	close(IN);	
}

1;