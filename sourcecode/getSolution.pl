###################################################################################################################################
##	KDDANet: Main function for uncovering roadmap of gene interactions mediating KDDAs
##	Wroted by Hua Yu and Lu Lu on 2018-09-01
####################################################################################################################################

use strict;
use Getopt::Long;
require("lpsolveModeling.pl");
require("readFile.pl");

####### get options from command line
my ($drFile,$diFile,$netFile,$dtFile,$dgFile,$ddFile,$outPath,$gammaMin,$gammaMax,$gammaStep,$optimethod,$help);			
GetOptions(	
        "drug_info|dr:s" => \$drFile,	##	A tab-delimited file gives drug info (default DrugBank drugs)
	"disease_info|di:s" => \$diFile,	##	A tab-delimited file gives disease info (default OMIM diseases) 
	"network_info|n:s" => \$netFile,		##	A tab-delimited file gives interactome network info (default HumanNet)
	"drug_target_info|dt:s" => \$dtFile,	##	A tab-delimited file gives drug target info (default DrugBank targets)
	"disease_gene_info|dg:s" => \$dgFile,	##	A tab-delimited file gives disease gene info (default OMIM genes)
        "drug_disease_info|dd:s" => \$ddFile,    ##	A tab-delimited file gives drug-disease associations (default CTD database and literature)
	"out_path|o:s", => \$outPath,		##	Output path of prediction results (default ./results/)
	"gamma_min|gmin:f" => \$gammaMin, 	##	minimum of gamma values (default 4)
	"gamma_max|gmax:f", => \$gammaMax,	##	maximum of gamma values	(default 12)
	"gamma_step|gstep:f", => \$gammaStep,	## step increment of gamma (default 1)
	"opti_method|opt:s", => \$optimethod,	## optimization method (SDrTDi or SDiTDr, default Source_Drug_Target_Disease (SDrTDi))
	"help|h!" => \$help,	## display help info
) or usage();

usage() if(defined $help);

####################### USAGE;
sub usage {
    print <<"END_USAGE";
Usage:
perl KDDANet.pl -dr <DRUG_FILE> -di <DISEASE_FILE> -n [NETWORK_FILE] -dt [DRUG_TARGET_FILE] -dg [DISEASE_GENE_FILE] -o [OUT_PATH] -gmain [GAMMA_MIN] -gmax [GAMMA_MAX] -gstep [GAMMA_INCREAMENT] -h
-dr  DRUG_FILE for drug info (default DrugBank drugs)
-di  DISEASE_FILE for disease info (default OMIM diseases)
-n  NETWORK_FILE for network info (default HumanNet, a functional gene association network '../Data/HumanNet.txt')
-dt  DRUG_TARGET_FILE for drug target info (default DrugBank targets)
-dg  DISEASE_GENE_FILE for disease gene info (default OMIM genes)
-dd  DRUG_DISEASE_FILE for drug-disease associations (default drug-disease associations extracted from CTD database)
-o  OUT_PATH for output files (default '../result/')
-gmin  MIN_GAMMA for minimum value of gamma (default 4)
-gmax  MAX_GAMMA for maximum value of gamma (default 10)
-gstep  GAMMA_INCREMENT (default 1)
-opt  CONTEXT (SDrTDi or SDiTDr, default Source_Drug_Target_Disease (SDrTDi))
-h  display help information
END_USAGE
    exit;
}

if(!defined $netFile){
	my $netFile = '../inputdata/HumanNet.txt';
}

if(!defined $outPath){
	my $outPath = '../result/';
}

if(!defined $gammaMin){
	my $gammaMin = 4;
}

if(!defined $gammaMax){
	my $gammaMax = 12;
} 	
if(!defined $optimethod){
	my $optimethod = "SDrTDi";
} 
my @gammas = SettingGamma($gammaMin,$gammaMax,$gammaStep);

#########################  main code using test data ###############################
pNet($drFile,$diFile,$netFile,$dtFile,$dgFile,$ddFile,$outPath,\@gammas,$optimethod);

sub pNet{
	my($drFile,$diFile,$netFile,$dtFile,$dgFile,$ddFile,$outPath,$gammas,$optimethod) =@_;
	
	print "DrugInfoFile:\t$drFile\n";
	print "DiseaseInfoFile:\t$diFile\n";
	print "NetFile:\t$netFile\n";
	print "DrugTargetFile:\t$dtFile\n";
	print "DiseaseGeneFile:\t$dgFile\n";
        print "DrugDiseaseFile:\t$ddFile\n";
	print "OutPath:\t$outPath\n";
	
	foreach my $gamma (@{$gammas}){
		my $outDir = "$outPath/$gamma";
		mkdir($outDir);
		LPModeler($drFile,$diFile,$netFile,$dtFile,$dgFile,$ddFile,$outDir,$gamma,$optimethod);
	}
}

#######################################################################################
##	Setting a range of gamma values
##	Hua Yu and Lu Lu 2017-09-01
#######################################################################################
sub SettingGamma{
	my($min,$max,$step) = @_;
	my @gammaArray = ();	
	my $gamma = $min;
	for(my $i=0;$gamma<=$max;$i++){
		$gammaArray[$i] = $gamma;
		$gamma = $gamma+$step;
	}
	return(@gammaArray);	
}

#############################################################################################
##	Implementing minimum cost flow optimization algorithm by Linear Programming Model
#############################################################################################

sub LPModeler{
	my($drFile,$diFile,$netFile,$dtFile,$dgFile,$ddFile,$outDir,$gamma,$optimethod) = @_;
	my(%drHash,%diHash,%dtHash,%dgHash,%ddHash);
	readNodeFile($drFile,\%drHash);
	readNodeFile($diFile,\%diHash);
	readEdgeFile($dtFile,\%dtHash);
	readEdgeFile($dgFile,\%dgHash);
    readEdgeFile($ddFile,\%ddHash);
	print "pNet optimization (gamma=$gamma) ...\n";
	if($optimethod eq "SDrTDi"){
		foreach my $drug (keys %drHash){
			if(!-e "$outDir/prim.simp.result.$drug.$gamma.txt" || -z "$outDir/prim.simp.result.$drug.$gamma.txt"){
				my $stat = lpsolveModeling($drug,\%drHash,\%diHash,\%dtHash,\%dgHash,\%ddHash,$netFile,$outDir,$gamma,$optimethod);
				if($stat == 0){
					print "$drug is successfully submitted\n";
				}else{
					print "$drug is not successfully submitted\n";
				}
			}
		}
	}else{
		foreach my $disease (keys %diHash){
			if(!-e "$outDir/prim.simp.result.$disease.$gamma.txt" || -z "$outDir/prim.simp.result.$disease.$gamma.txt"){
				my $stat = lpsolveModeling($disease,\%diHash,\%drHash,\%dgHash,\%dtHash,\%ddHash,$netFile,$outDir,$gamma,$optimethod);
				if($stat == 0){
					print "$disease is successfully submitted\n";
				}else{
					print "$disease is not successfully submitted\n";
				}
			}
		}
	}
}

1;
