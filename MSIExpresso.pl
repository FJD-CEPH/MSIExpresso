#!/usr/bin/perl
# MSI Search Tools V0.8 2018
# Emmanuel Tubacher
# Fondation Jean Dausset-CEPH
# Read the configuration file and store all parameters in hash table and in variables for easier access
my %parameters;
readConfigFile();
checkParameters();
my $outpath="$parameters{'output dir'}/$parameters{'output files root'}";
print STDERR "Using output directory $outpath for all files\n";
if ($parameters{'output files root'} ne "") {
	$outpath=$outpath.".";
}
if ($ARGV[1] eq "updateStatus") {
	rerunMSIStatus();
	exit(0);
}
my $samout=$parameters{'bam files'};
my $readfiles=$parameters{'read info files'};

# Search for coding microsatellites

if ($parameters{'coding'} eq "yes") {
	searchCodingMicrosats($outpath,$samout,$readfiles);
}

# Search for intronic microsatellites
my %MSIintronicMicrosatellites;
if ($parameters{'intronic'} eq "yes") {
        %MSIintronicMicrosatellites=searchIntronicMicrosat($outpath,$samout,$readfiles);
}

# Search for exon skipping

if ($parameters{'skipping'} eq "yes") {
	searchExonSkipping($outpath,$samout,$readfiles,%MSIintronicMicrosatellites);
}

#if ($parameters{'check NMD'} eq "yes") {
#	searchNMD($outpath,$samout,$readfiles);
#}

sub rerunMSIStatus {
	reassignStatus("CODING");
	my $codingFile="$parameters{'output dir'}/$parameters{'output files root'}.CODING.counts.csv";
	getMsiStatus($codingFile);
	reassignStatus("INTRONIC");
	reassignSkippingStatus();
}

sub reassignStatus {
	my ($type)=@_;
	my $codingFile="$parameters{'output dir'}/$parameters{'output files root'}.$type.counts.csv";
	print STDERR "Computing MSI Status for $codingFile\n";
	open(COUNTS,$codingFile) || die "No $codingFile, unable to set MSIstatus\n";
	my $tmpFile=$codingFile;
	$tmpFile=~s/csv/tmp/;
	open(TMP,"> $tmpFile") ||die "Can't write $tmpFile\n";
	my $header=<COUNTS>;
	print TMP $header;
	while (my $line=<COUNTS>) {
		chomp $line;
		my ($id,$base,$length,$gene,$stat,$sco,@genos)=split(/\t/,$line);
		my ($status,$score)=getMicrosatelliteStatus($length,@genos);
		print TMP "$id\t$base\t$length\t$gene\t$status\t$score";
		foreach my $geno (@genos) {
			print TMP "\t$geno";
		}
		print TMP "\n";
	}
	close COUNTS;
	close TMP;
	system("mv $tmpFile $codingFile");
	print STDERR "Done\n";
}

sub getStatusScore {
	my ($msi,$tot)=@_;
	my $msiStatus="NORMAL";
	my $score = -10;
	if ($tot >= $parameters{'minimum coverage'}){
		if  ($msi >0) {
			$score=log($msi*$msi/$tot)/log(10);
		}
	}
	else {
		$msiStatus="LOW EXPRESSION";
	}
	if ($score >= $parameters{'MSI score cutoff'}) {
		$msiStatus="MSI";
	}
	if ($score <  $parameters{'MSI score cutoff'} && $score > $parameters{'MSI score low cutoff'}) {
		$msiStatus="MSI LOW EXPRESSION";
	}
	return ($msiStatus,$score);
}

sub reassignSkippingStatus {
	my $skippingFile="$parameters{'output dir'}/$parameters{'output files root'}.SKIPPING.counts.csv";
	my $intronicFile="$parameters{'output dir'}/$parameters{'output files root'}.INTRONIC.counts.csv";
	open(INTRONIC,"cut -f1,5,6 $intronicFile|") || die "Can't open $inotrnicFile\n";
	my $hd=<INTRONIC>;
	my %intronStatus;
	my %intronScore;
	while (my $intron=<INTRONIC>) {
		chomp $intron;
		my ($i,$st,$sc)=split(/\t/,$intron);
		$intronStatus{$id}=$st;
		$intronScore{$id}=$sc;
	}
	close INTRONIC;
        open(COUNTS,$skippingFile) || die "No $skippingFile, unable to set MSIstatus\n";
        my $tmpFile=$skippingFile;
        $tmpFile=~s/csv/tmp/;
        open(TMP,"> $tmpFile") ||die "Can't write $tmpFile\n";
        my $header=<COUNTS>;
        print TMP $header;
        while (my $line=<COUNTS>) {
                my ($id,$base,$length,$gene,$nm,$ex,$dist,$totex,$skip,$c5,$c3,$stat,$sco,$istat,$isco)=split(/\t/,$line);
		if ($intronStatus{$id} eq "") {
			$intronStatus{$id}="NE";
			$intronScore{$id}=-10;
		}
		my ($msiStatus,$score)=getStatusScore($skip*5,$skip+$c5+$c3);
		print TMP "$id\t$base\t$length\t$gene\t$nm\t$ex\t$dist\t$totex\t$skip\t$c5\t$c3\t$msiStatus\t$score\t".$intronStatus{$id}."\t".$intronScore{$id}."\n";
	}
	close COUNTS;
	close TMP;
	system("mv $tmpFile $skippingFile");
}
sub searchNMD {
	my ($out,$sam,$read)= @_;
	my $minln=$parameters{'minimum overlap length'}+int($parameters{'database sequence length'}/2);
	open (NMD,"> ${out}NMDcounts.csv") || die "Can't create ${out}NMDcounts.csv\n";
	if ($sam eq "yes") {
	open(SAMOUT,"|$parameters{'path to samtools'} view -Sb -t $parameters{'MMR database'}.fai - > ${out}NMD.bam") || die "Can't create ${out}NMD.bam\n";
	}
	if ($read eq 'yes') {
		open(READS,"> ${out}readsNMD.csv")||die "Can't create ${out}readsNMD.csv\n";
	}
	my %countNMD;
        my $bwaNumThreads=$parameters{'bwa threads'};
        if ($bwaNumThreads ne "") { $bwaNumThreads=" -t $bwaNumThreads";}
	my $bwaNMDCmd="$parameters{'path to bwa'} mem $bwaNumThreads $parameters{'MMR database'} $parameters{'fastq read1'}  $parameters{'fastq read2'} | $parameters{'path to samtools'} view -S -F4 -q 30 - ";
        open(IN,"$bwaNMDCmd|") ||die "Command $bwaNMDCmd failed\n";
        while (my $line = <IN>) {
                chomp $line;
                my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,$options)=split(/\t/,$line);
		if (skippingMatch($cigar,$minln)) {
			$countNMD{$rname}++;
			if ($read eq 'yes') {
				print READS "$qname\t$rname\n";
			}
			if ($samout eq "yes") {
				print SAMOUT "$line\n";
			}
		}
	}
	close IN;
	if ($read eq 'yes') {
		close READS;
	}
	if ($samout eq "yes") {
		close SAMOUT;
	}
	foreach my $nmd (sort keys %countNMD) {
		print NMD "$nmd\t$countNMD{$nmd}\n";
	}
	close NMD;	
}	

sub searchCodingMicrosats {
	my ($out,$sam,$read)= @_;
	# Count of WT reads for each sequence containing a coding microsatellite
	my %countNL;
	# Count of reads for each sequence containing a coding microsatellite shorter than the reference (deletion)
	my %countDel;
	# Count of reads for each sequence containing a coding microsatellite longer than the reference (inserttion)
	my %countIns;
	#Hash table storing the fasta description of matching references sequences
	my %names;
	# Saving reads id and matching sequence to a file (optional)
	if ($read eq 'yes') {
		open(READS,"> ${out}readsCoding.csv")||die "Can't create ${out}readsCoding.csv\n";
	}
	# Saving bwa alignment for each matching read in sam format (optional
	if ($sam eq 'yes') {
		open(SAMOUT,"|$parameters{'path to samtools'} view -Sb -t $parameters{'database dir'}/CODING.fa.fai - > ${out}coding.bam")||die "Can't create ${out}coding.bam\n";
	}
	# bwa command line filtering out reads with no match and with a matching quality < 30
	my $bwaNumThreads=$parameters{'bwa threads'};
	if ($bwaNumThreads ne "") { $bwaNumThreads=" -t $bwaNumThreads";}
	my $bwaCodingCmd="$parameters{'path to bwa'} mem $bwaNumThreads $parameters{'database dir'}/CODING.fa $parameters{'fastq read1'}  $parameters{'fastq read2'} | $parameters{'path to samtools'} view -S -F4 -q 30 - ";
	open(IN,"$bwaCodingCmd|") ||die "Command $bwaCodingCmd failed\n";
	while (my $line = <IN>) {
		chomp $line;
		# sam line parsing
		my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,$options)=split(/\t/,$line);
		# fasta description parsing
		#my %microsats=parseFastaDescription($rname);
		# reference length of microsatellite
		my $msln=$microsat;
		$cigar=~s/H/S/g;
		$msln=~s/[A-Z]//g;
		# Minimum alignment length according to microsatellite length and position in reference sequence
		my $minmatch=$msln+$posInSeq+1;
		# Determining read type (WT, insertion, deletion of microsatellite) according to the CIGAR string
		my $score=getTypeFromCigar($cigar,$minmatch);
		# Empty score indicates a match shorter than minimum length
		if ($score eq '') {next;}
		# Score equals to 0 indicates a WT read
		if ($score ==0) {
			$countNL{$rname}++;
			$names{$rname}=1;
		# Negative score indicates a microsatellite deletion and its relative size e.g : -1 for a single base deletion
		} elsif ($score < 0) {
			$countDel{$rname}=$countDel{$rname}."#$score";
			$names{$rname}=1;
		# Positive score indicates a microsatellite insertion and its relative size e.g : 1 for a single base insettion
		} elsif ($score eq 'I') {
			$countIns{$rname}=$countIns{$rname}."#$score";
			$names{$rname}=1;
		}
		# Detailed read information saved to the reads file
		if ($read eq "yes" && $score ne '') {
			print READS "$qname\t$rname\t$score\n";
		}
		# Alignment saved to sam file
		if ($sam eq 'yes') {
			print SAMOUT "$line\n";
		}
	}
	close IN;
	# Closing optional output files
	if ($read eq "yes") { close READS;}
	if ($sam eq 'yes') {close SAMOUT;}
	# Exporting counts for each reference sequence and each type of read : WT, insertion, deletion detailed by sizes
	open (OUT,"> ${out}CODING.counts.csv ") ||die "Unable to write file : ${out}CODING.counts.csv\n";
	print OUT "Microsatellite id (ch.position)\tbase\tlength\tAnnotation\tMSI Status\tScore\tSizes distribution\n";
	foreach my $name (sort keys %names) {
		my $desc=getCodingMicrosatDesc($name);
		my ($id,$base,$wtlen,@reste)=split(/\t/,$desc);
		if ($countNL{$name} eq "") { $countNL{$name}=0;}
		my @delgeno=split(/#/,$countDel{$name});
		my @insgeno=split(/#/,$countIns{$name});
		my @numdelgeno;
		my @numinsgeno;
		my $maxdel=0;
		my $mindel=0;
		# Maximum  insertion and deletion sizes
		for (my $i=1;$i<=$#delgeno;$i++) {
			$delgeno[$i]=~s/\-//;
			$numdelgeno[$delgeno[$i]]++;
			if ($delgeno[$i] > $maxdel && $delgeno[$i] < $wtlen/2) { $maxdel=$delgeno[$i];}
		}
		for (my $i=1;$i<=$#insgeno;$i++) {
			$numinsgeno[$insgeno[$i]]++;	
			if ($insgeno[$i] > $maxins && $insgeno[$i] < $wtlen*2) { $maxins=$insgeno[$i];}
		}
		if ($maxins ==0 && $maxdel==0) {next;}
		# Making human readable description of reference sequence
		print OUT "$desc";
		# Printing reads count for each deletion size found
		my $genos="";
		my @genosizes;
		for (my $i=$maxdel;$i>0;$i--) {
			if ($numdelgeno[$i] > 0) { 
				$genos=$genos."\t-$i $numdelgeno[$i]";
				push(@genosizes,"-$i $numdelgeno[$i]");
			}
		}
		# Printing reads count for WT reads
		$genos=$genos."\t0 $countNL{$name}";
		push(@genosizes,"0 $countNL{$name}");
		# Printing reads count for each insertion size found
		for (my $i=1;$i<=$maxins;$i++) {
			if ($numinsgeno[$i] > 0) { 
				$genos=$genos."\t+$i $numinsgeno[$i]";
				push(@genosizes,"+$i $numinsgeno[$i]");
			}
		}
		my($status,$score)=getMicrosatelliteStatus($wtlen,@genosizes);
		print OUT "\t$status\t$score$genos\n";
	}
	close OUT;
	getMsiStatus("${out}CODING.counts.csv");
}

sub getMicrosatelliteStatus {
	my ($ln,@gn)=@_;
        my %sizecount;
        my $sc=0;
        my $totindel=0;
        my $totsize=0;
	# Cutoff for microsatellite indel size filtering
        my $cutoff=15;
        foreach my $g (@gn) {
                my ($size,$cnt)=split(/\s/,$g);
                $sizecount{"$size"}=$cnt;
                if ($size != 0) {
                        $totindel=$totindel+$cnt;
                }
        }
        if ($totindel > $parameters{"minimum indel count"})  {
                foreach my $size (sort keys(%sizecount)) {
                        if ($sizecount{"$size"}/$totindel >= .1 && abs($size) > 0) {#$ln/$cutoff) {
                                $totsize=$totsize+abs($size)*$sizecount{"$size"};
                        }
                }
        }
	my ($st,$sc)=getStatusScore($totsize,$sizecount{"0"}+$totindel);
	return ($st,$sc);
}

sub getMsiStatus {
	my ($codingFile) = @_;
	my $totalExpressed=0;
	my $totalMSI=0;
	my $list;
	print STDERR "Computing sample MSI status on $codingFile\n";
	open(IDS,$parameters{'UTR panel'}) || die "No UTR panel found\n";
	while (my $idline = <IDS> ) {
		chomp $idline;
		my($id,@reste)=split(/\t/,$idline);
		open(IN," grep -w $id $codingFile |") || die "No $codingFile found.\n";
		while (my $line = <IN>) {
			my ($id,$base,$ln,$gene,$stat,$score,@genos)=split(/\t/,$line);
			if ($stat eq "MSI") {
				$list=$list.$line;
				$totalMSI++;
				$totalExpressed++;
			}
			elsif ($stat eq "NORMAL") {
				$totalExpressed++;
			}
		}
		close IN;
	}
	close IDS;
	my $outFile=$codingFile;
	$outFile=~ s/CODING/MSIStatus/;
	my $ratio=0;
	print STDERR "Writing $outFile\n";
	open(OUT,"> $outFile") || die "Unable to write $outFile\n";
	if ($totalExpressed > 0) {$ratio=100*$totalMSI/$totalExpressed;}
	if ($totalExpressed < $parameters{'minimal expressed UTR'}) {
		print OUT "Status : NA Not enough microsatellites expressed $totalExpressed\n$list";
	}
	elsif ($ratio >= $parameters{'MSI-H cutoff'}) {
		print OUT "Status : MSI-H $ratio% ($totalMSI/$totalExpressed)\n$list";
	}
	elsif ($ratio >= $parameters{'MSI-L cutoff'}) {
		print OUT "Status : MSI-L $ratio% ($totalMSI/$totalExpressed)\n$list";
	}
	elsif ($ratio < $parameters{'MSI-L cutoff'}) {
                print OUT "Status : MSS $ratio% ($totalMSI/$totalExpressed)\n$list";
        }
	close OUT;
}

sub getCodingMicrosatDesc {
	my ($ds)=@_;
	my $msdesc;
	my ($ms,$d1,@desc)=split(/\//,$ds);
	my %types=('C' => 'coding', 'U' => "5'utr", 'D' => "3'utr");
	my ($id,$base,$ln,$pos)=split(/#/,$ms);
	my ($nm,$gn,$exnum,$totex,$type)=split(/#/,$d1);
	my $msdesc="$id\t$base\t$ln\t$nm-$gn-$exnum/$totex $types{$type}";
	foreach my $d (@desc) {
		my ($nm,$gn,$exnum,$totex,$type)=split(/#/,$d);
		$msdesc=$msdesc."#$nm-$gn-$exnum/$totex $types{$type}";
	}
	return $msdesc;
}

sub getIntronicMicrosatDesc {
	my ($ds)=@_;
	my $msdesc;
	my ($d1,$d2,@desc)=split(/\//,$ds);
	my ($id,$base,$ln,$pos)=split(/#/,$d1);
	my ($nm,$gn,$intnum,$totex,$dist)=split(/#/,$d2);
	$msdesc="$id\t$base\t$ln\t$nm-$gn-$intnum$dist/$totex";
	foreach my $d (@desc) {
		my ($nm,$gn,$intnum,$totex,$dist)=split(/#/,$d);
		$msdesc=$msdesc."#$nm-$gn-$intnum$dist/$totex";
	}
	return $msdesc;
}
sub searchExonSkipping {
	my ($out,$sam,$read,%msintronic)=@_;
	#minimum match size to make sure a read is overlaping the junction
	my $minln=$parameters{'minimum overlap length'}+int($parameters{'database sequence length'}/2);
	# Saving reads id and matching sequence to a file (optional)
	if ($read eq "yes") {
		open(READS,"> ${out}readsSkipping.csv")||die "Can't create ${out}readsSkipping.csv\n";
	}
	# Saving bwa alignment for each matching read in sam format (optional)
	if ($sam eq "yes") {
	open (SAMOUT,"|$parameters{'path to samtools'} view -Sb -t $parameters{'database dir'}/SKIPPING.fa.fai - > ${out}skipping.bam")||die "Can't create ${out}skipping.bam\n";
	}
	my %skippingCounts;
	# bwa command line filtering out reads with no match and with a matching quality < 30
	my $bwaNumThreads=$parameters{'bwa threads'};
        if ($bwaNumThreads ne "") { $bwaNumThreads=" -t $bwaNumThreads";}
	my $bwaSkippingCmd="$parameters{'path to bwa'} mem $bwaNumThreads $parameters{'database dir'}/SKIPPING.fa $parameters{'fastq read1'} $parameters{'fastq read2'} 2> /dev/null | $parameters{'path to samtools'} view -S -F4 -q 30 - ";
	open(IN,"$bwaSkippingCmd |") ||die "Command $bwaSkippingCmd failed\n";
	my %nmSkip;
	while (my $line = <IN>) {
        	chomp $line;
		# sam line parsing
        	my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,$options)=split(/\t/,$line);
		# fasta description parsing
		my (@descs)=split(/\|/,$rname);
		my @keys;
		foreach my $ds (@descs) {
			my ($id,$base,$ln,$nm,$gn,$junction,$totex,$dist)=split(/[#,\/]/,$ds);
			my $key="$nm#$gn#$totex#$junction#$dist";
			push(@keys,$key);
			if ($nmSkip{$key} eq "") {
				$nmSkip{$key}="$id#$base#$ln";
			}
			elsif ($nmSkip{$key} !~ "$id" ) {
				$nmSkip{$key}=$nmSkip{$key}."-$id#$base#$ln";
			}
		}
		# Determining if read is presenting an exon skipping event according to the CIGAR string
		if (skippingMatch($cigar,$minln)) {
			foreach my $nj (@keys) {
				$skippingCounts{$nj}++;
			}
			# Detailed read information saved to the reads file
			if ($read eq "yes") { print READS "$qname\t$rname\n"; }
			# Alignment saved to sam file
			if ($sam eq "yes") { print SAMOUT "$line\n"; }
		}
	}
	close IN;
	# Closing optional output files
	if ($read eq "yes") { close READS; }
	if ($sam  eq "yes") { close SAMOUT;}
	# searching for normal exon-exon junctions 
	my %junctionsCounts;
	# bwa command line filtering out reads with no match and with a matching quality < 30
        my $bwaNumThreads=$parameters{'bwa threads'};
        if ($bwaNumThreads ne "") { $bwaNumThreads=" -t $bwaNumThreads";}
	my $bwaJunctionsCmd="$parameters{'path to bwa'} mem $bwaNumThreads $parameters{'database dir'}/JUNCTIONS.fa $parameters{'fastq read1'} $parameters{'fastq read2'} 2> /dev/null | $parameters{'path to samtools'} view -S -F4 -q 30 - ";
	open(IN,"$bwaJunctionsCmd |") ||die "Command $bwaJunctionsCmd failed\n";
	if ($read eq "yes") {
		open (READS,"> ${out}readsJunctions.csv")||die "Can't create ${out}readsJunctions.csv\n";
	}
	if ($sam eq "yes") {
        	open (SAMOUT,"|$parameters{'path to samtools'} view -Sb -t $parameters{'database dir'}/JUNCTIONS.fa.fai - >  ${out}junctions.bam") ||die "Can't create ${out}junctions.bam\n";
	}
	while (my $line = <IN>) {
        	chomp $line;
		# sam line parsing
        	my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,$options)=split(/\t/,$line);
		# fasta description parsing
        	my (@descs)=split(/\|/,$rname);
		my %nmj;
		foreach my $ds (@descs) {
			my ($id,$base,$ln,$nm,$gn,$totex,$junction)=split(/[#,\/]/,$ds);
			$nmj{"$nm#$gn#$totex#$junction"}=1;
		}
		# Determining if read is overlaping an WT exon-exon junction according to the CIGAR string
        	if (skippingMatch($cigar,$minln)) {
			foreach my $nj (keys %nmj) {
                		$junctionsCounts{"$nj"}++;
			}
			# Detailed read information saved to the reads file
                	if ($read eq "yes") { print READS "$qname\t$rname\n"; }
			# Alignment saved to sam file
              		if ($sam eq "yes") { print SAMOUT "$line\n";}
        	}
	}
	close IN;
	# Closing optional output files
	if ($read eq "yes") { close READS; }
	if ($sam  eq "yes") { close SAMOUT;}
	# Exporting counts for each reference sequence and each type of read (when at least one skipping event is found) :
	# normal #n-1/#n junction count, exon #n skipping count, normal #n/#n+1 junction count
	open (OUT,"> ${out}SKIPPING.counts.csv ") ||die "Unable to write file : ${out}SKIPPING.counts.csv\n";
	my %intronicStatus=getIntronicStatus($out);
	print OUT "microsatellite id (chr.position)\tbase\tlength\tGene name\tRefSeq id\tskipped exon\tposition of microsatellite\tnumber of exons in transcript\tnumber of skipping reads\tnumber of normal 5' junction reads\tnumber of 3' normal junction reads\tMSI Status\tScore\tintronic Status\tIntronic Score\n";
	foreach my $name (sort keys %skippingCounts) {
		my ($nm,$gn,$totex,$ex,$dist)=split(/#/,$name);
		my $prevnum=$ex-1;
		my $nexnum=$ex+1;
		my $prevname="$nm#$gn#$totex#$prevnum-$ex";
		my $nexname="$nm#$gn#$totex#$ex-$nexnum";
		if ($junctionsCounts{$prevname} eq "") { $junctionsCounts{$prevname}=0;}
		if ($junctionsCounts{$nexname} eq "") { $junctionsCounts{$nexname}=0;}
		# Making human readable description of reference sequence
		my $totreads=$junctionsCounts{$prevname}+$junctionsCounts{$nexname}+$skippingCounts{$name};
		my @ids=split(/\-/,$nmSkip{$name});
		my ($msiStatus,$score)=getStatusScore($skippingCounts{$name},$totreads);
		foreach my $idbl (@ids) {
			my ($id,$base,$length)=split(/#/,$idbl);
			if ($intronicStatus{$id} eq "") { $intronicStatus{$id}="NOT EXPRESSED\t-10";}
			print OUT "$id\t$base\t$length\t$gn\t$nm\t$ex\t$dist\t$totex\t$skippingCounts{$name}\t$junctionsCounts{$prevname}\t$junctionsCounts{$nexname}\t$msiStatus\t$score\t".$intronicStatus{$id}."\n";
		}
	}
	close OUT;
}

sub getIntronicStatus {
	my ($fileroot)=@_;
	my %res;
	open(INTRON,"cut -f1,5,6 ${fileroot}INTRONIC.counts.csv|") ||die;
	my $head=<INTRON>;
	while (my $line=<INTRON>) {
		chomp $line;
		my ($i,$st,$sc)=split(/\t/,$line);
		$res{$i}="$st\t$sc";
	}
	close INTRON;
	return %res;
}

sub searchIntronicMicrosat {
	my($out,$sam,$read)=@_;	
	my %MSI;
	if ($read eq "yes") {
		open(READS,"> ${out}readsIntronic.csv")||die "Can't create ${out}readsIntronic.csv\n";
	}
	if ($sam eq "yes") {
		open(SAMOUT,"|$parameters{'path to samtools'} view -Sb -t $parameters{'database dir'}/INTRONIC.fa.fai - >${out}intronic.bam")||die "Can't create ${out}intronic.bam\n";
	}
	my %countRTNL;
	my %countRTIns;
	my %countRTDel;
        my $bwaNumThreads=$parameters{'bwa threads'};
        if ($bwaNumThreads ne "") { $bwaNumThreads=" -t $bwaNumThreads";}
	my $bwaIntronicCmd="$parameters{'path to bwa'} mem $bwaNumThreads $parameters{'database dir'}/INTRONIC.fa  $parameters{'fastq read1'} $parameters{'fastq read2'} | $parameters{'path to samtools'} view -S -F4 -q 30 -";
	open(IN,"$bwaIntronicCmd |") ||die "Command $bwaIntronicCmd failed\n";
	while (my $line = <IN>) {
        	chomp $line;
        	my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,$options)=split(/\t/,$line);
		my $desc=getIntronicMicrosatDesc($rname);
        	my $minmatch=int($parameters{'database sequence length'}/2)+$parameters{'minimum overlap length'};
        	my $score=getTypeFromCigar($cigar,$minmatch);
		if ($score eq '') { next;}
        	if ($score == 0) {
                	$countRTNL{$desc}++;
                	$names{$desc}=1;
        	} elsif ($score < 0) {
                	$countRTDel{$desc}=$countRTDel{$desc}."#$score";
                	$names{$desc}=1;
        	} elsif ($score > 0) {
                	$countRTIns{$desc}=$countRTIns{$desc}."#$score";
                	$names{$desc}=1;
        	}
		if ($read eq "yes") { print READS "$qname\t$desc\t$score\n"; }
		if ($sam  eq "yes") { print SAMOUT "$line\n"; }
	}
	close IN;
	if ($read eq "yes") {
		close READS;
	}
	if ($sam eq "yes") {
		close SAMOUT;
	}
	open (OUT,"> ${out}INTRONIC.counts.csv ") ||die "Unable to wrtie file ${out}INTRONIC.counts.csv\n";
	print OUT "Microsatellite id (ch.position)\tbase\tlength\trefseq annotation (NM, gene name, exon number,relative distance/total number of exons)\tMSI Status\tScore\tSizes distribution\n";
	foreach my $name (sort keys %names) {
		my ($i,$b,$l,@reste)=split(/\t/,$name);
		if ($countRTNL{$name} eq "") { 
			$countRTNL{$name}=0;
		}
        	my @delgeno=split(/#/,$countRTDel{$name});
       	 	my @insgeno=split(/#/,$countRTIns{$name});
        	my @numdelgeno;
        	my @numinsgeno;
        	my $maxdel=0;
        	my $mindel=0;
        	for (my $i=1;$i<=$#delgeno;$i++) {
                	$delgeno[$i]=~s/\-//;
                	$numdelgeno[$delgeno[$i]]++;
                	if ($delgeno[$i] > $maxdel && $delgeno[$i] < $l/2) { $maxdel=$delgeno[$i]; }
        	}
        	for (my $i=1;$i<=$#insgeno;$i++) {
                	$numinsgeno[$insgeno[$i]]++;
               		if ($insgeno[$i] > $maxins && $insgeno[$i] < $l*2) { $maxins=$insgeno[$i]; }
        	}
		if ($maxdel==0 && $maxins==0) {next;}
        	print OUT "$name";
                my $genos="";
                my @genosizes;
                for (my $i=$maxdel;$i>0;$i--) {
                        if ($numdelgeno[$i] > 0) {
                                $genos=$genos."\t-$i $numdelgeno[$i]";
                                push(@genosizes,"-$i $numdelgeno[$i]");
                        }
                }
                $genos=$genos."\t0 $countRTNL{$name}";
                push(@genosizes,"0 $countRTNL{$name}");
                for (my $i=1;$i<=$maxins;$i++) {
                        if ($numinsgeno[$i] > 0) {
                                $genos=$genos."\t+$i $numinsgeno[$i]";
                                push(@genosizes,"+$i $numinsgeno[$i]");
                        }
                }
                my($status,$score)=getMicrosatelliteStatus($wtlen,@genosizes);
		if ($status eq "MSI") { $MSI{$name}=$score;}
                print OUT "\t$status\t$score$genos\n";

	}
	close OUT;
	return %MSI;
}
sub skippingMatch {
	my ($cig,$min)= @_;
	my @values=split(/[M,I,D,S,H]/,$cig);
        my @tags=split(/[0-9][0-9]*/,$cig);	
	if (($cig=~ /(^[0-9]*[H,S][0-9]*M[0-9]*S$)/ || $cig=~ /(^[0-9]*[H,S][0-9]*M$)/) && $values[1] >= $min) {
		return 1;
	} elsif ($cig=~ /(^[0-9]*M[0-9]*[H,S]$)/ && $values[0] >= $min) {
		return 1;
	}
	return 0;
}
sub getTypeFromCigar {
	my ($cig,$min)=@_;
	my $matchln=0;
	my $type='0';
	my @values=split(/[M,I,D,S,H]/,$cig);
	my @tags=split(/[0-9][0-9]*/,$cig);
	if (($cig=~ /(^[0-9]*[H,S][0-9]*M[0-9]*[H,S]$)/ || $cig=~ /(^[0-9]*[H,S][0-9]*M$)/) && $values[1] >= $min) {
		return 0;
	} elsif ($cig=~ /(^[0-9]*M[0-9]*[H,S]$)/ && $values[0] >= $min) {
		return 0;
	} else {
		for (my $i=1; $i <=$#tags ; $i++) {
			if ($tags[$i] eq 'M') {$matchln=$matchln+$values[$i-1];}
			if ($tags[$i] eq 'D') {$matchln=$matchln+$values[$i-1]; $type=0-$values[$i-1];}
			if ($tags[$i] eq 'I') {$matchln=$matchln+$values[$i-1]-1; $type=$values[$i-1];}
		}
	}
	if ($matchln < $min) {$type='';} 
	return $type;
}

sub readConfigFile {
	print STDERR "Reading configuration file\n";
	if ($ARGV[0] eq "") {die "Usage : ./MSIExpresso.pl <configuration file>\n";}
        open(IN,"grep -v '#' $ARGV[0]|") ||die "No configuration file $ARGV[0]!\n";
        while (my $line = <IN>) {
                chomp $line;
                if ($line =~ "#" || $line eq "") { next;}
                my ($key,$value)=split(/\=/,$line);
                $parameters{"$key"}=$value;
		print STDERR "$key=$value\n";
        }
        close IN;
}
sub checkParameters {
if ($ARGV[1] ne "updateStatus") {
	if (! -f $parameters{'fastq read1'}) {
		die "$parameters{'fastq read1'} fastq file not found\n";
	}
	if ($parameters{'fastq read2'} ne "") {
		if (! -f $parameters{'fastq read2'}) {
			die "$parameters{'fastq read2'} fastq file not found\n";
		}
	}
	if ($parameters{'minimum overlap length'} eq "") {
		$parameters{'minimum overlap length'}=10;
		print STDERR "minimum overlap length parameter no set, default value of 10bp will be used\n";
	}
}
	if ($parameters{"MSI-H cutoff"} eq "") {
		$parameters{"MSI-H cutoff"}=30;
		print STDERR "MSI-H cutoff parameter no set, default value of 30% will be used\n";
	}
	if ($parameters{"MSI-L cutoff"} eq "") {
		$parameters{"MSI-H cutoff"}=10;
		print STDERR "MSI-L cutoff parameter no set, default value of 10% will be used\n";
	}
	if ($parameters{"minimal expressed UTR"} eq "") {
		$parameters{"minimal expressed UTR"}=30;
		print STDERR "minimal expressed UTR parameter no set, default value of 30 will be used\n";
	}
	if ($parameters{"minimum coverage"} eq "") {
		$parameters{"minimum coverage"}=15;
		print STDERR "minimum coverage parameter no set, default value of 15 will be used\n";
	}
	if ($parameters{"database sequence length"} eq "") {
		$parameters{"database sequence length"}=60;
	}
}

1;
