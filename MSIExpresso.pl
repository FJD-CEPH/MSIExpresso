#!/usr/bin/perl
# MSI Search Tools V2.0 2025
# Emmanuel Tubacher
# Fondation Jean Dausset-CEPH
# Read the configuration file and store all parameters in hash table and in variables for easier access
my %parameters;
readConfigFile();
checkParameters();
my $outpath="$parameters{'output dir'}/$parameters{'output files root'}";
print STDERR "Using output directory $parameters{'output dir'}, output all files names starting with $parameters{'output files root'}\n";
if ($parameters{'output files root'} ne "") {
	$outpath=$outpath.".";
}

if ($parameters{'MMR status'} eq 'y' || $parameters{'MMR status'} eq 'yes') {
	MMRStatus();
}
if ($ARGV[1] eq "panelOnly") {
	my $panelDB=makePanelDB();
	%MSIintronicMicrosatellites=searchCodingMicrosats($outpath,$samout,$readfiles);	
	rmPanelDB();
	rerunMSIStatus();
	exit(0);
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

sub rerunMSIStatus {
	reassignStatus("CODING");
	my $codingFile="$parameters{'output dir'}/$parameters{'output files root'}.CODING.counts.csv";
	getMsiStatus($codingFile);
	if ($ARGV[1] ne 'panelOnly') {
		reassignStatus("INTRONIC");
		reassignSkippingStatus();
	}
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
	open(INTRONIC,"cut -f1,5,6 $intronicFile|") || die "Can't open $intronicFile\n";
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
if ($ARGV[1] eq "panelOnly") {
                print STDOUT "Computing MSI Status only on $parameters{'output dir'}/$parameters{'panel DB'}\n";
                $bwaCodingCmd="$parameters{'path to bwa'} mem $bwaNumThreads $parameters{'output dir'}/$parameters{'panel DB'}  $parameters{'fastq read1'} $parameters{'fastq read2'} | $parameters{'path to samtools'} view -S -F4 -q 30 -";
print STDOUT "$bwaCodingCmd";
        }

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
	print STDERR "Reading sample configuration file\n";
	if ($ARGV[0] eq "") {die "Usage : ./MSInspector.pl <sample configuration file>\n";}
	open(IN,"grep 'global configuration' $ARGV[0]| cut -f2 -d'='|");
	my $line=<IN>;
	if ($line ne '') {
		chomp $line;
		print STDERR "Reading global configuration file: $line\n";
		open(GLOBAL,"grep -v '#' $line|") ||die "No global configuration file $line!\n";
		while (my $line=<GLOBAL>) {
			chomp $line;
			my ($key,$value)=split('=',$line);
                        $parameters{$key}=$value;
			print STDERR "$key=$value\n";
                }
                close GLOBAL;			
	}	
	close IN;
        open(IN,"grep -v '#' $ARGV[0]|") ||die "No configuration file $ARGV[0]!\n";
	print STDERR "Reading sample configuration file $ARGV[0]!\nDuplicates entries from global configuration file will be overriden\n"; 
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
			print STDERR "minimum overlap length parameter not set, default value of 10bp will be used\n";
		}
	}
	open(BWA,"which $parameters{'path to bwa'}|") ||die ("bwa not found, please set and check path to bwa in configuration file\n");
	my $bwa=<BWA>;
	chomp $bwa;
	if ($bwa eq'') {
		die "bwa not found, please set and check path to bwa in configuration file\n";
	}
	close BWA;
	open(SAM,"which $parameters{'path to samtools'}|") ||die "samtools not found, please set and check path to samtools in configuration file\n";
	my $sam=<SAM>;
	chomp $sam;
	if ($sam eq'') {
		die "samtools not found, please set and check path to samtools in configuration file\n";
	}
	close SAM;
	open(BCF,"which $parameters{'path to bcftools'}|") ||die "bcftools not found, please set and check path to bcftools in configuration file\n";
	my $bcf=<BCF>;
	chomp $bcf;
	if ($bcf eq'') {
		die "bcftools not found, please set and check path to bcftools in configuration file\n";
	}
	if ($parameters{"MSI-H cutoff"} eq "") {
		$parameters{"MSI-H cutoff"}=30;
		print STDERR "MSI-H cutoff parameter not set, default value of 30% will be used\n";
	}
	if ($parameters{"MSI-L cutoff"} eq "") {
		$parameters{"MSI-H cutoff"}=10;
		print STDERR "MSI-L cutoff parameter not set, default value of 10% will be used\n";
	}
	if ($parameters{"minimal expressed UTR"} eq "") {
		$parameters{"minimal expressed UTR"}=30;
		print STDERR "minimal expressed UTR parameter not set, default value of 30 will be used\n";
	}
	if ($parameters{"minimum coverage"} eq "") {
		$parameters{"minimum coverage"}=15;
		print STDERR "minimum coverage parameter not set, default value of 15 will be used\n";
	}
	if ($parameters{"database sequence length"} eq "") {
		$parameters{"database sequence length"}=60;
	}
}

sub makePanelDB {
	%panelMS;
	open(PN,$parameters{"UTR panel"}) ||die "UTR Panel $parameters{'UTR panel'} not found\n";;
	my $dbname="paneldb$$.fa";
	while (my $line=<PN>) {
		chomp $line;
		my @pn=split(/\t/,$line);
		$panelMS{">$pn[0]"}=1;
	}
	close PN;
	my $MSCount=0;
	open(DB,"$parameters{'database dir'}/CODING.fa") || die;
	open(DBPANEL,"> $parameters{'output dir'}/$dbname") || die;
	while (my $line=<DB>) {
		chomp $line;
		my @data=split(/#/,$line);
		if ($panelMS{$data[0]}==1) {
			print DBPANEL "$line\n";
			my $seq=<DB>;
			print DBPANEL $seq;
			$panelMS{$data[0]}=0;
			$MSCount++;
		}
	}
	close DB;
	if ($MSCount < $#panelMS) {
		open(DB,"$parameters{'database dir'}/INTRONIC.fa") || die;
		while (my $line=<DB>) {
			chomp $line;
			my @data=split(/#/,$line);
			if ($panelMS{$data[0]}==1) {
				print DBPANEL "$line\n";
				my $seq=<DB>;
				$panelMS{$data[0]}=0;
				$MSCount++;
			}
		}
	}
	close DBPANEL;
	if ($MSCount < $#panelMS) {
		foreach my $id (sort keys %panelMS) {
			if ($panelMS{$id}==1) {print STDERR "$id not in INTRONIC.fa nor CODING.fa sequences databases\n";}
		}
	}
	system("$parameters{'path to bwa'} index $parameters{'output dir'}/$dbname");
	print STDERR "PANEL database $parameters{'output dir'}/$dbname created and bwa indexed\n";
	$parameters{'panel DB'}=$dbname;
	return $dbname;
}
sub rmPanelDB {
	unlink($arameters{"output dir"}."/".$parameters{'panel DB'});
}
sub MMRStatus {
	my %cdsLength;
	open(CDS, "$parameters{'database dir'}/MMR.cds.fa") ||die;
	        my @lines;
        my @heads;
        my @gns;
        my $nline=0;
        while (my $line=<CDS>) {
                chomp $line;
                if ($line =~'>') {
                        $line=~s/\>//;
                        push (@gns,$line);
                        push (@heads,$nline);
                }
                push (@lines,$line);
                $nline++;

        }
        push(@heads,$nline);
        close CDS;
        for (my $i=0; $i<= $#gns; $i++) {
                my $seq='';
                for (my $j=$heads[$i]+1; $j < $heads[$i+1]; $j++) {
                        $seq=$seq.$lines[$j];
                }
                $cdsLength{$gns[$i]}=length($seq);
        }
	my $outBam="$parameters{'output dir'}/$parameters{'output files root'}_MMR.bam";
	my $outVCF="$parameters{'output dir'}/$parameters{'output files root'}_MMR.vcf";
	print STDERR "Writing $outBam and $outVCF\n";
	my $mmrCommand="$parameters{'path to bwa'} mem $bwaNumThreads $parameters{'database dir'}/MMR.cds.fa $parameters{'fastq read1'}  $parameters{'fastq read2'} |$parameters{'path to samtools'} view -S -h -F4 -q 30 - |$parameters{'path to samtools'} sort -o $outBam ";
	system($mmrCommand);
	my $mmrCommand2="$parameters{'path to samtools'} view $outBam |cut -f3|sort |uniq -c|awk '{print \$2\" \"\$1}'";
	open(COUNTS,"$mmrCommand2|") ||die;
	my $m6=0;
	my $m2=0;
	my $m1=0;
	my $tot=0;
	my $p6=0;
	my $p2=0;
	my $p1=0;
	while (my $line=<COUNTS>) {
		chomp $line;
        	my @m=split(/\s/,$line);
        	if ($m[0] eq 'MSH6') {
			$m6=$m[1]; }
        	if ($m[0] eq 'MSH2') {
        		$m2=$m[1];
		}
        	if ($m[0] eq 'MLH1') {
       			$m1=$m[1];
		}
	}
	close COUNTS;
	$m6=$m6*1000/$cdsLength{'MSH6'};
	$m2=$m2*1000/$cdsLength{'MSH2'};
	$m1=$m1*1000/$cdsLength{'MLH1'};
	$tot=$m1+$m2+$m6;
        $p6=100*$m6/$tot;
        $p2=100*$m2/$tot;
        $p1=100*$m1/$tot;
	open(STATS,"> $parameters{'output dir'}/$parameters{'output files root'}_MMRStatus.txt")||die "Open $parameters{'output dir'}/$parameters{'output files root'}_MMRStatus.txt\n";
        if ($tot < 300) {print STAT "Not enough reads to estimate MMR deficiency accurately !\n";}
	print STATS "Relative expression of MLH1, MSH2 and MSH6\n";
	print STATS "MLH1\t$p1\t";
	my $st=1;
	if ($p1 <16) {print STATS "MLH1 low expression"; $st=0;} else { print STATS "OK";}
	print STATS "\nMSH2\t$p2\t";
	if ($p2 <16) {print STATS "MSH2 low expression"; $st=0;} else { print STATS "OK";}
	print STATS "\nMSH6\t$p6\t";
	if ($p6 <11) {print STATS "MSH6 low expression"; $st=0;} else { print STATS "OK";}
	print STATS "\n";
	if ($st==1) {print STATS "MMR genes normal expression\n";} else { print STATS "MMR genes loss of expression\n";}
	my $variantCmd="$parameters{'path to bcftools'} mpileup -f $parameters{'database dir'}/MMR.cds.fa $outBam | bcftools call -mv -Ov -o $outVCF";
	system($variantCmd);
	my @variants=searchFunctionalVariants($outVCF);	 
	if ($#variants <0) {
		print STATS "No functional variant found in MMR genes\n";
	}
	else {
		print STATS "List of functional variants in MMR genes\nGene name\tCDS position\tVaraitn description\n";
		foreach my $v (@vars) {
			$v=~s/#/\t/g;
			print STATS "$v\n";
		}
	}
	close STATS;
	print STDOUT "MMR status available: $parameters{'output dir'}/$parameters{'output files root'}_MMRStatus.txt\n";
}
sub searchFunctionalVariants {
	my ($vcf)=@_;
	my %gseq;
	open(SEQ,"$parameters{'database dir'}/MMR.cds.fa") ||die "No MMMR database found\n";
	my $sq='';
	my $gname='';	
	while (my $line=<SEQ>) {
		chomp $line;
		if ($line =~'>') {
			if ($gname ne '' && $sq ne '') {
				$gseq{$gname}=$sq;
				$gname=$line;
				$gname=~s/>//g;
				$sq='';
			}
		}
		else {
			$sq=$sq.$line;
		}
	}
	close SEQ;
	$gseq{$gname}=$sq;
	my (%codons_table) = ( 'TTT' => 'F', 'TTC' => 'F', 'TTA' => 'L', 'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
    'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', 'ATG' => 'M', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V', 'TCT' => 'S', 'TCC' => 'S',
    'TCA' => 'S', 'TCG' => 'S', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P', 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
    'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TAT' => 'Y', 'TAC' => 'Y', 'TAA' => '.', 'TAG' => '.', 'CAT' => 'H', 'CAC' => 'H',
    'CAA' => 'Q', 'CAG' => 'Q', 'AAT' => 'N', 'AAC' => 'N', 'AAA' => 'K', 'AAG' => 'K', 'GAT' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E',
    'TGT' => 'C', 'TGC' => 'C', 'TGA' => '.', 'TGG' => 'W', 'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGT' => 'S', 'AGC' => 'S',
    'AGA' => 'R', 'AGG' => 'R', 'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G',);
	my @vars;
	open(VCF,"cut -f1,2,4,5,8 $vcf|grep -v '#'|") || die;
	while (my $line=<VCF>) {
        	chomp $line;
        	my ($gene,$pos,$all1,$all2,$dp)=split(/\t/,$line);
       		my @vals=split(';',$dp);
        	my $nba1=0;
       		my $nba2=0;
        	foreach my $val (@vals) {
                	if ($val =~'DP4') {	
				$val=~s/DP4=//;
                        	my ($f1,$r1,$f2,$r2)=split(',',$val);
                        	$nba1=$f1+$r1;
                        	$nba2=$f2+$r2;
                        	if ($nba2 > $nba1) {
                                	my $ta=$all1;
                                	$all1=$all2;
                                	$all2=$ta;
                        	}
                	}
        	}
		if ($nba1+$nba2 < 20) {next;}
        	my $gnseq="";
        	my $utr3="";
        	my $cds="";
        	my $utr5="";
                $gnseq=$gseq{$gene};
        	($utr5,$cdutr)=split('atg',$gnseq);
        	$cdutr='ATG'.$cdutr;
        	if ($gnseq =~ 'taa') {
                	($cds,$utr3)=split('taa',$cdutr);
                	$cds=$cds.'TAA';
        	}
        	if ($gnseq =~ 'tag') {
                	($cds,$utr3)=split('tag',$cdutr);
                	$cds=$cds.'TAG';
        	}
        	if ($gnseq =~ 'tga') {
                	($cds,$utr3)=split('tga',$cdutr);
                	$cds=$cds.'TGA';
        	}
        	my $u5l=length($utr5);
        	my $u3l=length($utr3);
        	$cds=$cds.uc(substr($gnseq,3,$u5l+3+length($cds)));
        	my $cdl=length($cds);
        	if ($pos <= $u5l || $pos >= $u5l+$cdl) {next;}
		my $indel=length($all1)-length($all2);
		my $idl=abs($indel);
		if ($indel !=0) {
			my $deltype='insertion';
			if ($indel >0 ) {$deltype='deletion';}
			if ($indel%3 !=0) {
				push(@vars,"$gene#$pos#frameshift $deltype of $idl bases");
			}
			else {
				push(@vars,"$gene#$pos#codon(s) $deltype of $idl bases");
			}
			next;
		}
print "VAR $#vars $vars[$#vars]\n";
        	my $base=substr($cds,$pos-$u5l-1,length($all1));
        	my $env=substr($cds,$pos-$u5l-3,length($all1)+4);
        	my $genv=substr($gnseq,$pos-3,length($all1)+4);
        	if ($base ne $all1 && $base ne $all2) {print "$gene $pos $base ne $all1\n";next;}
        	else {
                	my $picodon=($pos-$u5l-1)%3;
                	my $oricodon=substr($cds,($pos-$u5l-1)-$picodon,3);
                	my $altcodon="";
                	for (my $i=0;$i<3; $i++) {
                        	if ($i==$picodon) {
                                	$altcodon=$altcodon.$all2;
                        	}
                        	else {
                                	$altcodon=$altcodon.substr($oricodon,$i,1);
                        	}
                	}
                	my $protpos=int(($pos-$u5l)/3)+1;
                	if ($codons_table{$altcodon} ne $codons_table{$oricodon}) {
				if ($codons_table{$altcodon} eq '.') {
					push(@vars,"$gene#$pos#".$codons_table{$oricodon}."$protpos -> STOP CODON ($oricodon / $altcodon)\n");
				}
				else {
                			push(@vars,"$gene#$pos#".$codons_table{$oricodon}."$protpos".$codons_table{$altcodon}."($oricodon / $altcodon)\n");
				}
			}
        	}
	
	}
	return @vars;
}
1;
