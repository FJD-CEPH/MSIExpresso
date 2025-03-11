#!/usr/bin/perl
my %parameters;
open(IN,$ARGV[0]) ||die "No configuration file found\n";
while (my $line=<IN>) {
	chomp $line;
	my($key,$val)=split(/\=/,$line);
	$parameters{$key}=$val;
}
close IN;
my %codingList=getMSISpecificMicrosatellites('CODING');
my %intronicList=getMSISpecificMicrosatellites('INTRONIC');
my %skippingList=getMSISpecificSkippingEvents();
sub getSamplesList {
	my ($msi,$mss);
	open(LS,"ls $parameters{'output directory'}/*MSIStatus.counts.csv|") ||die;
	while (my $samplepath=<LS>) {
		chomp $samplepath;
		my @pth=split(/\//,$samplepath);
		my $sample=$pth[-1];
		open(HEAD,"head -1 $samplepath|cut -f3 -d' '|") ||die "Can't read status in $samplepath\n";
		my $status=<HEAD>;
		close HEAD;
		$sample=~s/.MSIStatus.counts.csv//;
		if ($status =~ "MSI") {
			$msi="$msi|$sample";
		}
		elsif ($status =~"MSS") {
			$mss="$mss|$sample";
		}
	}
	close LS;
	$mss=~s/\|//;
	$msi=~s/\|//;
	return ($mss,$msi);
}

sub getMSISpecificSkippingEvents {
	open(LS,"ls $parameters{'output directory'}/*MSIStatus.counts.csv|") ||die;
	my $msiSamplesCount=0;
	my $mssSamplesCount=0;
	my $mssSamples="";
	my $msiSamples="";
        while (my $samplepath=<LS>) {
                chomp $samplepath;
                open(HEAD,"head -1 $samplepath|cut -f3 -d' '|") ||die "Can't read status in $samplepath\n";
                my $status=<HEAD>;
                close HEAD;
                $samplepath=~s/MSIStatus/SKIPPING/;
                if ($status =~ "MSI") {
                        $msiSamples="$msiSamples $samplepath";
                        $msiSamplesCount++;
                }
                elsif ($status =~"MSS") {
                        $mssSamples="$mssSamples $samplepath";
                        $mssSamplesCount++;
                }
        }
	open(IDS,"cut -f1,4,5,6,12 $msiSamples | grep -w MSI | grep -v Status|") ||die "Can't select MSI id from $type files\n";	
	my %skipId;
	my %skipCountMsi;
	my %multiId;
	my $prevgnm="";
	while (my $line=<IDS>) {
		chomp $line;
		my ($id,$gn,$nm,$numex,$status)=split(/\t/,$line);
		if ($status=~"MSI" && $prevgnm ne "$gn,$nm,$numex") {
			$skipCountMsi{"$gn,$nm,$numex"}++;
			$skipId{"$gn,$nm,$numex"}=$id;
			$prevgnm="$gn,$nm,$numex";
		}
	}
	close IDS;
	my %multiId;
	foreach my $gnm (keys %skipId) {
		my $MSSCount=0;
		if ($skipCountMsi{$gnm} >= $msiSamplesCount*$parameters{'minimal MSI percentage'}/100) {
			my $id=$skipId{$gnm};
			if ($multiId{$id} == 1) { next;}
			$multiId{$id}=1;
			open(MSS,"grep -w $id $mssSamples | cut -f12 | grep -w MSI| wc -l|") ||die "Can't read status from $mssSamples\n";
			while (my $line=<MSS>) {
				chomp $line;
				$MSSCount=$line;
			}
			close MSS;
			if ($MSSCount/$mssSamplesCount <= $parameters{'maximal MSS percentage'}/100) {
				$gn=$gnm;
				$gn=~s/\,/ /g;
				print "$skipId{$gnm}\t$gn\t$skipCountMsi{$gnm}\t$MSSCount\tSKIPPING\n";
			}
		}
	}
}
sub getMSISpecificMicrosatellites {
	my ($type)=@_;
	my %counts;
	my @ids;
	my $msiSamplesCount=0;
	my $mssSamplesCount=0;
	my %msiCounts;
	my %msids;
	my $totalSamples;
	my $msiSamples="";
	my $mssSamples="";
	open(LS,"ls $parameters{'output directory'}/*MSIStatus.counts.csv|") ||die;
	while (my $samplepath=<LS>) {
		chomp $samplepath;
		open(HEAD,"head -1 $samplepath|cut -f3 -d' '|") ||die "Can't read status in $samplepath\n";
		my $status=<HEAD>;
		close HEAD;
		$samplepath=~s/MSIStatus/$type/;
		if ($status =~ "MSI") {
			$msiSamples="$msiSamples $samplepath";
			$msiSamplesCount++;
		}
		elsif ($status =~"MSS") {
			$mssSamples="$mssSamples $samplepath";
			$mssSamplesCount++;
		}
	}
	my $fld="1,5";
	open(IDS,"cut -f$fld $msiSamples | grep -w MSI | cut -f1 | sort| uniq -c| grep -v Status| awk '{printf \"%s\t%d\\n\",\$2,\$1}'|") ||die "Can't select MSI id from $type files\n";
	my $mssCount=0;
	while (my $line=<IDS>) {
		chomp $line;
		my ($id,$count)=split(/\t/,$line);
		if ($count >= $msiSamplesCount*$parameters{'minimal MSI percentage'}/100) {
			open(IDS2," grep -wh $id $mssSamples | cut -f1,$fld | grep -w MSI | wc -l|") ||die "Can't select MSI id from $type files\n";
			$MSScount=<IDS2>;
			chomp $MSScount;
			close IDS2;
			if ($MSScount/$mssSamplesCount <= $parameters{'maximal MSS percentage'}/100) {
				my $desc=getDescription($id,$type);
				print "$id\t$count\t$MSScount\t$type\t$desc\n";
			} 
		}
	}
	close IDS;
}
sub getDescription {
	my ($i,$t)=@_;
	my $fl="2-4";
	if ($t eq "SKIPPING") {
		$fl="2-7";
	}
	open(DESC,"grep -w $i $parameters{'output directory'}/*$t.counts.csv|cut -f$fl| head -1|") || die;
	my $desc=<DESC>;
	close DESC;
	chomp $desc;
	return $desc;
	
}
1;
