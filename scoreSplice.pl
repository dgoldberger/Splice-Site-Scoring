#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/perl_5.10.1/bin/perl
use warnings;
use strict;
use Bio::EnsEMBL::Registry;

use lib '/home/unix/davegold/lib/';

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'calcium.broadinstitute.org', # alternatively 'useastdb.ensembl.org'
    -user => 'mlek',
    -pass => 'ml3k'
);


open TF, "ThreePFreqs.txt" or die "can't open needed file: $!";
open ACC, "AcceptorDistuptMin3Ratios.txt" or die "can't open needed file: $!"; #n=231
open DI, "DiNuc2.txt" or die "can't open needed file: $!";
open FH, ">vcfMixed.vcf";

open NOT, "DonorNotGtDisRat.txt" or die "can't open needed file: $!"; # n=1075
open SC, "randomSNPRatios.txt" or die "can't open needed file: $!"; # n=35928 , all non splice SNP with AF>1%

open SNP, "DonorDiNucSNPRatio.txt" or die "can't open needed file: $!"; # n=942
open CRT, "spliceCreateRatios.txt" or die "can't open needed file: $!"; # n=84
open ASNP, "AcceptorSNPRatios.txt" or die "can't open needed file: $!"; #n=243


my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );

#print "#CHROM\tPOS\tREF\tALT\tINFO\tEnsembl_Exon_ID\tSeq\n";

my %negStrand= (A=>"T", C=>"G", T=>"A", G=>"C"); # used to convert bases from neg strand later

#thresholds to determine if files are annotated or not
my $donorThresh=.6; #only annotate donor sites if they score above
my $donCreateThresh=.4; #only annotate donor creations if they score above 
my $accThresh=0; # only annotate acceptor sites if they score above 


#make array to hold all values of splice creating variants' ratios
my $p=0;
my @create=();
while (<CRT>){
    chomp $_;
    $create[$p]=$_;
    $p++;
}

#array holds values of common acceptor site variants ratios
my $z=0;
my @accSNP=();
while (<ASNP>){
    chomp $_;
    $accSNP[$z]=$_;
    $z++;
}

#sliding window ratios for non splice site SNP
my $o=0;
my @randSNP=();
while (<SC>){
    chomp $_;
    $randSNP[$o]=$_;
    $o++;
}

#make arrays to keep all ratio values of non GT disrupting mutations and common splice SNP

my $k=0;
my @notGt=();
while (<NOT>){
    chomp $_;
    $notGt[$k]=$_;
    $k++;
}


#same for set of common SNP in donor splice sites
my $n=0;
my @snp=();
while (<SNP>){
    chomp $_;
    $snp[$n]=$_;
    $n++;
    }

#array of hashes for 3' single nucleotide frequencies
my @threeFreq=();
for (my $i=0; $i<5; $i++){
    $threeFreq[$i]={};
}
my $track=0;
while(<TF>){
    chomp $_;
    if (length($_)<2){
	$track++;
	next;
    }
    else{
	my @info=split(/\s+/, $_);
	chomp $info[0];
	chomp $info[1];
	$threeFreq[$track]{$info[0]}=$info[1];
    }
}


#array of ratios of known three prime disrupting variants
my @threeRatios=();
my $b=0;
while (<ACC>){
    chomp $_;
    $threeRatios[$b]=$_;
    $b++;
}



#make array of hashes for dinucleotide frequencies, to be used in scoring method
my @freq=();
for (my $i=0; $i<7;$i++){
    $freq[$i]={};
}
my $num=0;
while(<DI>){
    chomp $_;
    if (length($_)<3){
	$num++;
	next;
    }
    else{
	my @info=split(/\s+/, $_);
	chomp $info[0];
	chomp $info[1];
	$freq[$num]{$info[0]}=$info[1];
    }
}

#PWM scoring for acceptor site
sub scoreThreeP {
    my $seq=shift;
    chomp $seq;
    $seq=uc($seq);
    my $total=0;
    for (my $k=0; $k<length($seq); $k++){
        my $nt=substr($seq,$k,1);
        my $value=$threeFreq[$k]{$nt};
	my $mod= (log($value/.25))/log(2);
        $total=$total+ $mod;
    }
    $total=$total+31.06098; # min is -31.06098, to make all positive, i shift
    return $total;

}

# takes in acceptor site with variant like AA[C/G]GT and finds ratios between scores
sub threePRatio{
    my $seq=shift;
    $seq=uc($seq);
    chomp $seq;
    my $start=index($seq, '[');
    my $end=index($seq, ']');
    my $first;
    if ($start==0){
        $first="";
    }
    else{
        $first=substr($seq,0,$start);
    }
    if ($start==-1){
	return 0;
    }
    else{
	my $ref=$first.substr($seq, $start+1,1).substr($seq,$end+1,length($seq)-($end+1));
	my $alt= $first.substr($seq, $start+3,1).substr($seq,$end+1,length($seq)-($end+1));
	my @refs=[];
	my @alts=[];
	for (my $i=0;$i<length($alt)-4;$i++){
	    $alts[$i]=scoreThreeP(substr($alt,$i,5));
	}
	for(my $i=0;$i<length($ref)-4;$i++){
	    $refs[$i]=scoreThreeP(substr($ref,$i,5));
	}
	my $maxRef=-1000;
	foreach (@refs){
	    if ($_>$maxRef){
		$maxRef=$_;
	    }
	}
	my $maxAlt=-1000;
	foreach (@alts){
	    if ($_> $maxAlt){
		$maxAlt=$_;
	    }
	}
	return $maxAlt/$maxRef;
    }
	
}

#converts 3p ratio into a score less than one based on distribution of ratios
sub ThreeRatioToScore{
    my $ratio=shift;
    chomp $ratio;
    my $freq=0; #will be evaluated to frequency of ratios smaller than one given
#score change 0819 
  # foreach (@threeRatios){
#	chomp $_;
#        if ($_<=$ratio){
#            $freq++;
#        }
#    }
    foreach (@accSNP){
        chomp $_;
        if ($_<=$ratio){
            $freq++;
        }
    }
    $freq=$freq/243.0; # there are 243 elements in @accSNP
    my $score=1-($freq);
    return $score;
}


#use dinuc frequencies above to calculate score. for each position score is: log_2(freq/.0625)
sub scoreFiveP{
    my $seq=shift;
    chomp $seq;
    $seq=uc($seq);
    my $total=0;
    for (my $k=0; $k<length($seq)-1; $k++){
        my $nt=substr($seq,$k,2);
        my $value=$freq[$k]{$nt};
        my $mod= (log($value/.0625))/log(2);
        $total=$total+ $mod;
    }
    $total=$total+35.129672; # min value is -35.129672, makes all positive
    return $total;
}





#scores reference site and alt site using method above. sliding window around site, produces ratio
#ratios is between highest scores produced by sliding window for donor site
sub makeRatio{
    my $seq=shift;
    $seq=uc($seq);
    my $start=index($seq, '[');
    my $end=index($seq, ']');
    my $first;
    if ($start==0){
        $first="";
    }
    else{
        $first=substr($seq,0,$start);
    }
    my $ref=$first.substr($seq, $start+1,1).substr($seq,$end+1,length($seq)-($end+1));
    my $alt= $first.substr($seq, $start+3,1).substr($seq,$end+1,length($seq)-($end+1));

    my @refs=[];
    my @alts=[];
    for (my $i=0;$i<length($alt)-7;$i++){
	$alts[$i]=scoreFiveP(substr($alt,$i,8));
    }
    for(my $i=0;$i<length($ref)-7;$i++){
	$refs[$i]=scoreFiveP(substr($ref,$i,8));
    }
    my $maxRef=-1000;
    foreach (@refs){
        if ($_>$maxRef){
            $maxRef=$_;
        }
    }
    my $maxAlt=-1000;
    foreach (@alts){
        if ($_> $maxAlt){
            $maxAlt=$_;
        }
    }

    return $maxAlt/$maxRef;

}


#using know ratios of splice disrupting variants and SNP this method conversts ratio into score
#score= 1-(percent SNP ratios less than given)
sub ratioToScore{
    my $ratio=shift;
    chomp $ratio;
    my $notGt1=0; #will be evaluated to frequency of ratios smaller than one given
    my $snp1=0;
#commented out score change, 08/19
  #  foreach (@notGt){
#	chomp $_;
#	if ($_<=$ratio){
	#    $notGt1++;
#	}
 #   }
    
  #  $notGt1=$notGt1/1075.0;

    foreach (@snp){
	chomp $_;
	if ($_<=$ratio){
	    $snp1++;
	}
    }
    $snp1=$snp1/942.0; #942 donor SNPs

    my $score=1-$snp1;

    return $score;
    
}

#ratio maker for splice creation evaluation
sub makeRatioForCreate{
    my $ref=shift;
    my $alt=shift;
    $ref=uc($ref);
    $alt=uc($alt);


    my @refs=[];
    my @alts=[];
    for (my $i=0;$i<length($alt)-7;$i++){
        $alts[$i]=scoreFiveP(substr($alt,$i,8));
    }
    for(my $i=0;$i<length($ref)-7;$i++){
        $refs[$i]=scoreFiveP(substr($ref,$i,8));
    }
    my $maxRef=-1000;
    foreach (@refs){
        if ($_>$maxRef){
            $maxRef=$_;
        }
    }
    my $maxAlt=-1000;
    foreach (@alts){
        if ($_> $maxAlt){
            $maxAlt=$_;
        }
    }

    my @info=($maxAlt, $maxAlt/$maxRef);

    return @info;

}


#sub routine to test for splice creation, when exonic variant is not near a splice site
#produces score same as one above using splicing creating (not disrupting) data
sub spliceCreate {
    my $ref=shift;
    my $alt=shift;
    my @ratioInfo=makeRatioForCreate($ref,$alt);
    my $ratio=$ratioInfo[1];
    my $score=0;
#score change 0819
#    foreach(@create){
#	chomp $_;
#	if ($_<$ratio){
#	    $score++;
#	}
#    }
    foreach (@randSNP){
        chomp $_;
        if ($_<=$ratio){
            $score++;
        }
    }
    
    $score=$score/35928.0; #there are 35928 SNP in exons with AF>1%


    #if ref is really really low, and alt is twice as good, it should not be annotated.
    if ($ratioInfo[0]<31){ # alt score is low, regardless of ratio, the score should be 0
	$score=0;
    }
#    my $gt=index($alt, "GT");
 #   if ($gt==-1){
#	$score=0;
#    }
    #if (index($ref, "GT")==$gt){
#	$score=0;
#    }


    return $score;

}




#begin looking through VCF file
while(<>){
#my $seqOut;
if($_ =~ /^#/){
   print FH $_;
	next;
}

if (length($_)<3){
    next;
}
my $scoreSeq="";

#print "$_\n";
my @data = split(/\t/, $_);
my $chr = $data[0];
my $pos = $data[1];
my $ref = $data[3];
my $altNt = substr($data[4],0,1);

#some vcf have this field
chomp $data[7];
if ($data[7]=~/HGMD_CHANGE/){
     $data[7]=~/HGMD_CHANGE=(.*);/;
     $altNt=substr($1,2,1);
     $ref= substr($1,0,1);
}


my $posRef; #position relative to acceptor or donor site

#$data[7]=~/HGMD_TYPE=(\w+);/;
#my $type=$1;
#print "HGMD_TYPE=$type\n";
if ($data[7]=~/HGMD_POS=/){
    $data[7]=~/HGMD_POS=(-?\w+);/;
    $posRef=$1;
}

my $neg=0; #flag variable for negative strand transcripts for later
#flag variables for each type of site
my $acceptor=0;
my $donor=0;


#flags for variants not in splice sites to be assessed for creation
my $midExon=0;
my $midIntron=0;

my @exonPosData; # declared but only used for variants not in splice sites for splice creation evaluation
my @intronPosData=(0,-1,-1,0); # setting defaults so I can test for changes and it doesnt complain
#more flags
my $inExonSite=0; #splice site and exonic
my $inIntronSite=0; # splice site and intronic

my $inExon=0; # flag for exonic variants NOT in splice Sites. if $inExon==1 and $inExonSite==1 then its a non splice site exon variant
#inExonSite, several lines above, flags splice site variants in exons, to not be confused

my $inIntron=0; # intronic and not splice site

my $goodTrans=0; # flag to make sure that transcript contains variant

my $testSeq; # this will be tested for splice creating
#to determine if it is near a splice site and which one

#call genome slice based on VCF info
my $slice1 = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $pos, $pos );
my $genes1 = $slice1->get_all_Genes();
my $gene1;

#make sure gene I am using is protein coding
foreach (@{$genes1}){
    if ($_->biotype=~ /protein_coding/){
	$gene1 =$_;
	last;
    }
}
#if there is no protein coding gene at place (error in vcf?) then just print it out and go next
if (not $gene1){
    print FH join ("\t", @data), "\n";
    next;
}


my $transcripts1 = $gene1->get_all_Transcripts();
my $transNum=0; # keep track of transcipt number in array of transcripts for gene to be called later

foreach my $t1 (@{$transcripts1}){
    if ($t1-> biotype =~ /protein_coding/){
	my @exons1 = @{$t1->get_all_Exons()};
	my @introns1 = @{$t1->get_all_Introns()};
	
	#this checks if the variant is contained in a particular transcript by looking thru exons and introns
	for(my $i=0; $i < @exons1; $i++){
	    my $pos2=0; #posref for exon vairants not in splice sites
	    #check if variant is in exon or near start or end of exon
	    if (($exons1[$i]->start>=0 && $exons1[$i]->end<=0) || ($exons1[$i]->start<=0 && $exons1[$i]->end>=0)){
			$goodTrans++;
			$inExon++;
			if ($exons1[$i]->start<=0){
			$pos2=-1*($exons1[$i]->start)+1; # hard to understand but it works
		    }
			else{
			    $pos2=-1*($exons1[$i]->end)+1;
			}
			@exonPosData=($pos2, $i, $transNum, $t1->strand); #data array to be evaluated for splice creation only
			if ($exonPosData[0]>10){
			    last;
			}
			
		    }
	#just in case variant is intronic it needs to be flagged
	    if (abs($exons1[$i]->start)<8 || abs($exons1[$i]->end)<8){
		$goodTrans++;
		last;
	    }
	}


	    #also checking if variant is in intron corresponding to transcript 
	for(my $i=0; $i < @introns1; $i++){
	    my $pos3=0;
	    if (($introns1[$i]->start>=0 && $introns1[$i]->end<=0) || ($introns1[$i]->start<=0 && $introns1[$i]->end>=0)){
			$goodTrans++;
			$inIntron++;
			if ($introns1[$i]->start>-20 && $introns1[$i]->start<-8 && $introns1[$i]->end>=1){
			    $pos3=-1*($introns1[$i]->start)+1;
			}
			elsif ($introns1[$i]->end<20 && $introns1[$i]->end>8 && $introns1[$i]->start<=1){
			    $pos3=-1*($introns1[$i]->end)+1;
			}
			@intronPosData=($pos3, $i, $transNum, $t1->strand); # same as exonPosData above but for creation in intron regions
		    }
	    if ($intronPosData[0]!=0){
		last;
	    }

	}
	#following block determines which site and how far variant is 
	
	if ($goodTrans>0){
	    for(my $i=0; $i < @exons1; $i++){
	  	 if ($t1->strand==1){
		
	         #test if variant is within 3 bases upstream of acceptor (intron)
	    	if ($exons1[$i]->start>1 && $exons1[$i]->start<4 && $exons1[$i]->end>0){
			 $acceptor=1;
			 $posRef=($exons1[$i]->start);
			 $inIntronSite++;
			 last;
	    	}
	    	#test if variant is within 2 bp in the exon
	    	elsif ($exons1[$i]->start<=1 && $exons1[$i]->start>-5 &&$exons1[$i]->end>0){
		    $acceptor=1;
		    $posRef=-1*($exons1[$i]->start);
		    $inExonSite++;
		    $posRef++;
		    
		    last;
	    	}
	    	#test if variant is in the last 2 bp of the exon
	    	elsif($exons1[$i]->end>=1 && $exons1[$i]->end<3 && $exons1[$i]->start<0){
		    $donor=1;
		    $posRef=$exons1[$i]->end;
		    $inExonSite++;
		    #$posRef++;
		    
		    last;
	    	}
	    	#test if variant is in the 6 bp after an exon
	    	elsif($exons1[$i]->end<1 && $exons1[$i]->end>-7 && $exons1[$i]->start<0){
		    $donor=1;
		    $posRef=-1*$exons1[$i]->end;
		    $posRef++;
		    $inIntronSite++;
		    last;
            	}
	
		   
	    }
		#same as above but for negative strand so signs are switched
	    elsif($t1->strand==-1){
		
		$neg=1;
	
		if ($exons1[$i]->start>1 && $exons1[$i]->end>0 && $exons1[$i]->end>$exons1[$i]->start && $exons1[$i]->start<8){
		    $posRef=$exons1[$i]->start-1;
		    $donor=1;
		    $inIntronSite++;
		    last;
		}
		elsif ($exons1[$i]->start<1 && $exons1[$i]->end<=0 && $exons1[$i]->start<$exons1[$i]->end && abs($exons1[$i]->end<5)){
		    $posRef=$exons1[$i]->end-1;
		    $acceptor=1;
		    $inIntronSite++;
		    last;
		}

		#following block is always false, old method, left in just in case
#	    if (0){
	       #im going to flip all the start and stops
	    	#if ($exons1[$i]->end<0 && $exons1[$i]->end>-4 && $exons1[$i]->start>0){
		
#	       if ($exons1[$i]->end<0 && $exons1[$i]->end>-1*(length($exons1[$i]->seq->seq)+4) && $exons1[$i]->start>0){      
#		      $acceptor=1;
		
#		      $posRef=-1*($exons1[$i]->end);
#		      $inExonSite++;
#		      last;
#	    	}
	    	#test if variant is within 2 bp in the exon
	    	#elsif ($exons1[$i]->end>=0 && $exons1[$i]->end<5 && $exons1[$i]->start>0){
#		elsif ($exons1[$i]->end>=0 && $exons1[$i]->end<length($exons1[$i]->seq->seq)-5 && $exons1[$i]->start>0){ 
#	            $acceptor=1;
#		    print "acc2\n";
#		    print "exon start: ", $exons1[$i]->end, "\n";
#		    $posRef=($exons1[$i]->end)-1;
#		    $inIntronSite++;
#		    last;
	 #   	}
	    	#test if variant is in the last 2 bp of the exon
	 #   	elsif($exons1[$i]->start<0 && $exons1[$i]->start>-3 && $exons1[$i]->end<0){
	#	    $donor=1;
	#	    print "don1\n";
	#	    print "exon end: ", $exons1[$i]->start, "\n";
	#	    $posRef=-1*(abs($exons1[$i]->start)+1);
	#	    $inIntronSite++;
	#	    last;
	#    	}
	    	#test if variant is in the 6 bp after an exon
	 #   	elsif($exons1[$i]->start>=0 && $exons1[$i]->start<7 && $exons1[$i]->end<0){
	#	    $donor=1;
	#	    print "don2\n";
	#	    print "exon end: ", $exons1[$i]->start, "\n";
	#	    $posRef=$exons1[$i]->start+1;
		   
	#	    $inExonSite++;
	#	    last;
     #       	}
	 #
	   #}
		
	    }
	}

    }
      if ($goodTrans>0){ # if transcript contains variant, stop looking
	    last;
	}

    
    }

    $transNum++; #keeping track of which transcript I am using
}

#print "goodtran: $goodTrans   inIntronSite: $inIntronSite   inExonSite: $inExonSite  inIntron: $inIntron  intronposdata: $intronPosData[0]\n";

#corrector so intron variants are not flagged as in a site
if (abs($intronPosData[0])>6){
    $inIntronSite=0;
}

#determine position of non splice site variants in exon
if ($goodTrans>0 && $inIntronSite==0 && $inExonSite==0 && $inExon>0 && abs($exonPosData[0])>10){ # trans has it, not intron site, not exon site, in exon, more than 10 from SS
    $midExon++; #flag variable
    my @trans=@{$transcripts1};
    my $tran=$trans[$exonPosData[2]];
    my @exs=@{$tran->get_all_Exons()};
    my $ex=$exs[$exonPosData[1]];
    #rare case where exons arent defined for gene
    if (not $ex){
	print FH join ("\t", @data), "\n";
	next;
    }
    $testSeq=substr($ex->seq->seq, $exonPosData[0]-7,16);
    if ($exonPosData[3]==-1){
		$altNt=$negStrand{$altNt}; # exonposData[3] is strand, so i need to provide compliment
    }
#    print "exon create\n";
    my $altSeq= substr($testSeq,0,7).$altNt.substr($testSeq,8,8);
    my $createScore=spliceCreate($testSeq, $altSeq);
    if ($createScore>=$donCreateThresh){
	print FH join ("\t", @data), ";DONOR_CREATE=$createScore\n";
    }
    else{
	print FH join ("\t", @data), "\n";
    }
    next;
}
    

#checks for creation for intron variants
elsif ($goodTrans>0 && $inIntronSite==0 && $inExonSite==0 && $inIntron>0 && $intronPosData[0]!=0){
    $midIntron++;
    my @trans=@{$transcripts1};
    my $tran=$trans[$intronPosData[2]];
    my @ins=@{$tran->get_all_Introns()};
    my $in=$ins[$intronPosData[1]];

    if (not $in){
	print FH join ("\t", @data), "\n";
	next;
    }
    $testSeq=substr($in->seq, $intronPosData[0]-7,16);
    if ($intronPosData[3]==-1){
        $altNt=$negStrand{$altNt};
    }
 #   print "intron create\n";
    my $altSeq= substr($testSeq,0,7).$altNt.substr($testSeq,8,8);
    my $createScore=spliceCreate($testSeq, $altSeq);
    if ($createScore>$donCreateThresh){
	#print SC $createScore, "\n";
        print FH join ("\t", @data), ";DONOR_CREATE=$createScore\n";
    }
    else{
	print FH join ("\t", @data), "\n";
    }
    next;
}


#if not in splice site or middle of exon write to file and next
#this case would only happen for intronic variants more than 20 nt away from ss
if ($inIntronSite==0 && $inExonSite==0 && $midExon<1 && $midIntron<1){ # midExon and midIntron are flagged in two blocks above
    print FH join ("\t", @data), "\n";
    next;
}

#evaluate exonic splice variants
if ($inExonSite>0 && $inIntronSite==0 && abs($posRef)<8){
    #find protein coding gene at that position
    foreach my $gene (@{$genes1} ) {
	if (($gene->biotype)!~ /protein_coding/){
	    next;
	}
#go thru each trans until you find one that has an exon containing variant
    foreach my $t (@{$transcripts1} ) {
	if($t -> biotype =~ /protein_coding/){
		my @exons = @{$t->get_all_Exons()};
		my @introns = @{$t->get_all_Introns()};
	    
		my $exonNum = -1; # indicates variant is not found in particular exon
		
		for(my $i=0; $i < @exons; $i++){
		    	if($exons[$i]->start <= 0 && $exons[$i]->end >= 0){
				$exonNum = $i;
				last;
			}
			
			elsif($exons[$i]->start >=0 && $exons[$i]->end <= 0){
				$exonNum = $i;
				last;
			    }
			
		    }
		#if found...
		if($exonNum != -1){
			my $prevIntron = ".";
			my $nextIntron = ".";
			#bracket and lower case the variant to be detected later
			my @exon_bases = split(//,$exons[$exonNum]->seq->seq);
			if($t->strand == 1){
			       #$exon_bases[abs($exons[$exonNum]->start)+1] = "[".lc($exon_bases[abs($exons[$exonNum]->start)+1])."]";
			       if ($donor){
				     $exon_bases[-1*$posRef] = "[".lc($exon_bases[-1*$posRef])."]";
				
				     }
			       elsif($acceptor){
				   $exon_bases[$posRef-1] = "[".lc($exon_bases[$posRef-1])."]";
				    }
			   }
			  else{
				#$exon_bases[$#exon_bases-abs($exons[$exonNum]->start)-1] = "[".lc($exon_bases[$#exon_bases-abs($exons[$exonNum]->start)-1])."]";
				$exon_bases[$#exon_bases-abs($exons[$exonNum]->start)-1] = "[".lc($exon_bases[$#exon_bases-abs($exons[$exonNum]->start)-1])."]";
				
				#HGMD_CHANGE field accounts for neg strand, otherwise I need to "transcribe altNt" using hash defined earlier
				if ($data[7]!~/HGMD_CHANGE/){
				     $altNt=$negStrand{$altNt};
				    }
			   }	
			
				my $exonSeq = join("",@exon_bases);
			 


				# is this some exon other than the first exon
				if($exonNum - 1 >= 0){
					if(length($introns[$exonNum-1]->seq) >= 50){
						$prevIntron = lc(substr($introns[$exonNum-1]->seq,length($introns[$exonNum-1]->seq)-50,50));
					    }
					else{
						$prevIntron = lc($introns[$exonNum-1]->seq);
						
					    }
						
					
				    }
				# is this some exon other than the last exon
				if($exonNum+1 < (scalar @introns)){
					if(length($introns[$exonNum]->seq) >= 50){
					    $nextIntron = lc(substr($introns[$exonNum]->seq,0,50));
					}
					else{
					    $nextIntron = lc($introns[$exonNum]->seq);
					}

					#print "$nextIntron\n";
				    }
		   # }
			
				#donor site uses next intron                
				if ($donor>0){
				    $nextIntron=substr($nextIntron,0,9);
				    my $open=index($exonSeq, '[');
				    my $strLen=$open-(length($exonSeq)-4);
				    my $strEnd=length($exonSeq)-($open+2);
				    $scoreSeq=substr($exonSeq, length($exonSeq)-4, $strLen+2)."/$altNt".substr($exonSeq, $open+2,$strEnd).$nextIntron;

				}
			        #acceptor site uses previous intron
			  	elsif ($acceptor>0){
				    my $open=index($exonSeq, '[');
				    $prevIntron=substr($prevIntron,length($prevIntron)-4,4);
			  	    #$exon_introns=$prevIntron.$exonSeq;

			  	    $scoreSeq=$prevIntron.substr($exonSeq, 0, $open+2)."/$altNt".substr($exonSeq, $open+2,4);

			  	}		  
				
               
		    }

		if ($exonNum!=-1){
		    last; #if exon that contains variant is found, it stops looping through transcripts
			}
	    } 
	else{
	    next; #if a transcript is not protein coding, it goes on to the next one
	}
    }
    last; # looks at first gene that has variant
}
}


#when splice site variant is in intron
if ($inIntronSite>0 && $inExonSite<1 && abs($posRef)<8){
    foreach my $gene (@{$genes1} ) {
	if ($gene->biotype!~/protein_coding/){
	    next;
	}

	my $transcripts = $gene->get_all_Transcripts();

	foreach my $t (@{$transcripts} ) {
		my $intronNum=-1;
	    	if($t -> biotype =~ /protein_coding/){
		

		    my @exons = @{$t->get_all_Exons()};
		    my @introns = @{$t->get_all_Introns()};
		
		    for(my $i=0; $i < @introns; $i++){
		   
			#to check if exon contains the nt
			if($introns[$i]->start <= 1 && $introns[$i]->end >= 0){
				$intronNum = $i;
				last;
			    }
			elsif($introns[$i]->start >= 0 && $introns[$i]->end <= 0){
				$intronNum = $i;
				last;
			    }

		    }

		
		if ($intronNum!=-1){
		    my $prevExon=".";
		    my $nextExon=".";
		  #  $varPos=$posRef;
		    my @intron_Bases=split(//, $introns[$intronNum]->seq);
		    if ($t->strand == 1){
			#$intron_Bases[abs($introns[$intronNum]->start)+1]="[".lc($intron_Bases[abs($introns[$intronNum]->start)+1])."]"; 
			if ($acceptor){
			
			    $intron_Bases[-1*($posRef)+1]="[".lc($intron_Bases[-1*($posRef)+1])."]";
		    }
			elsif($donor){
			    $intron_Bases[$posRef-1]="[".lc($intron_Bases[$posRef-1])."]";
			}
		    }
		    elsif($t->strand==-1) {
			if ($acceptor){  
			    $intron_Bases[$#intron_Bases-(abs($posRef)-1)] = "[".lc($intron_Bases[$#intron_Bases-(abs($posRef)-1)])."]";
		    }
			elsif($donor){
			    #$intron_Bases[$#intron_Bases-(abs($posRef))] = "[".lc($intron_Bases[$#intron_Bases-(abs($posRef))])."]";
			    $intron_Bases[(abs($posRef)-1)] = "[".lc($intron_Bases[(abs($posRef)-1)])."]";
			}
			if ($data[7]!~/HGMD_CHANGE/){
			    $altNt=$negStrand{$altNt};

			}
			
		    }
		   
		    my $intronSeq=lc(join("", @intron_Bases));
		  
		    my $intronLen= length($intronSeq);
		    my $midString;
		    
		    if ($acceptor==1){
		    	$intronSeq=substr($intronSeq,length($intronSeq)-10,10);	
		    }
		    
		    elsif ($donor==1){
		    	$intronSeq=substr($intronSeq,0,14);	
		    }
		  
			
		    #$intronSeq= substr($intronSeq,0,10).$midString.substr($intronSeq,$intronLen-30,30);
			
		    if ($acceptor==1){
			$nextExon=substr($exons[$intronNum+1]->seq->seq, 0, 6);
			
			
			my $open=index($intronSeq, '[');
			my $len= $open-(length($intronSeq)-6);
			my $lenEnd=length($intronSeq)-($open+2);
			$scoreSeq=substr($intronSeq, length($intronSeq)-6, $len+2)."/$altNt".substr($intronSeq, $open+2,$lenEnd).$nextExon;
		    }
		    elsif ($donor==1){
				$prevExon=substr($exons[$intronNum]->seq->seq, length($exons[$intronNum]->seq->seq)-4, 4);
				my $open=index($intronSeq, '[');
				$scoreSeq=$prevExon.substr($intronSeq, 0, $open+2)."/$altNt".substr($intronSeq, $open+2,8);
			}
		}
		if ($intronNum!=-1){ #if the protein coding transcript contains the variant, stop 
			last;
		}	
	
	}	
		else{
			next;
		}
	    }
	
		last; #only one gene
    }
	
}

#catches some intronic edge cases, not within 20 nt of exon, limitation of exome data
#else{
#    print "elsed out\n";
 #   print "inIntron: $inIntron\n";
  #  print "scoreSeq: $scoreSeq\n";
  #  print "intronposdata: $intronPosData[0]\n";
#}



#to remove portions of sequence not adjecent to splice site
#make sure this corresponds to what I think it is


# to make sure I'm actually at a junction
#two caps ensures I am two nt into the exon
my $capCount=0;
for (my $j=0;$j<length($scoreSeq); $j++){
    if (substr($scoreSeq, $j,1) eq uc(substr($scoreSeq, $j,1))){
	$capCount++;
    }
}

#exons are capitalized and introns are lower case so this makes sure there are enough cap letters
#brackets and slashes are also "capital" so we need at least 5

if ($capCount<5){
    print FH join ("\t", @data), "\n";
    next;
}



if ($donor==1){
    my $ratio=makeRatio($scoreSeq);

    my $score=ratioToScore($ratio);
    if ($score>$donorThresh){
	print FH join("\t", @data).";DONOR_DISRUPT=$score\n";
    }

    else{
	print FH join("\t", @data)."\n";
    }

}
elsif($acceptor==1){
    if ($scoreSeq!~/\[/ || $scoreSeq!~/\]/){
	print FH join("\t", @data)."\n";
    }
    else{
	$scoreSeq=substr($scoreSeq, 0, length($scoreSeq)-3);
	my $ratio=threePRatio($scoreSeq);
  
	my $score=ThreeRatioToScore($ratio);
	if ($score>$accThresh){
	    print FH join ("\t", @data).";ACCEPTOR_DISRUPT=$score\n";
	    
	}
	else{
	    print FH join("\t", @data)."\n";
	}
    }
}


}

print "yabba dabba doo\n"; # just because...

close(TF);
close(FH);
close(SNP);
close(NOT);
close(DI);
close(CRT);
close(SC);
close(ACC);
