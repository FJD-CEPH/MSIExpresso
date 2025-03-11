#!/usr/bin/perl
my %parameters=checkParameters($ARGV[0]);

my @chrs;
my %msType;
open(IN,"cut -f1 -d'.' $parameters{'circos microsatellite list'}| sort -u|") ||die;
while (my $line = <IN>) {
	chomp $line;
	push(@chrs,$line);
}
close IN;
print STDERR "Found ".($#chrs+1)." chromosomes in $parameters{'circos microsatellite list'} file\n";
my @int;
my $chromint="";
my $size=$parameters{'microsatellite display size'};

my $plotCircos="$parameters{'circos data directory'}/$parameters{'circos plot name'}plot-circos.conf";
open(CIRCOS,"> $plotCircos") || die "Can't write $plotCircos\n";
print CIRCOS "karyotype = $parameters{'circos templates directory'}/karyotype.human.txt\n";
foreach my $chr (@chrs) {
	my $pint1=0;
	my $pint2=0;
	my $minpos=0;
	my $maxpos=0;
	my @intchr;
	open(CHR,"cut -f1 $parameters{'circos microsatellite list'} | grep -w $chr| sed s/'\\.'/' '/ |sort -nk2|") ||die "Unable to read  $parameters{'circos microsatellite list'} microsatellites list\n";
	print "Preparing data for chromosome $chr\n";
	while (my $line=<CHR>) {
		chomp $line;
		my ($ch,$pos)=split(/\s/,$line);
		if ($ch ne $chr) { next;}
		my $id="$ch.$pos";
		my $type=getMsType($id);
		if ($type eq "") {print STDERR "Skipping microsatellite $id not found in data files\n"; next;}
		if ($pos < $size/2) {
			$pint1=1;
			$pint2=$size+$pint1;
		}
		else {
			$pint1=$pos-$size/2;
			$pint2=$pint1+$size;
		}
		if ($#intchr < 0) {
			push(@intchr,"$id $pint1 $pint2 $type");
			$minpos=int(($pint1-$parameters{'microsatellite minimum distance'})/1000000);
			if ($minpos < 0) {$minpos=0;}
			$maxpos=int(($pint1+$parameters{'microsatellite minimum distance'})/1000000);
		}
		else {
			my ($pid,$p1,$p2)=split(/\s/,$int[$#int]);
			my ($chr,$pp)=split(/\./,$pid);
			if ($p2 > $pint1 && $ch eq $chr) {
				$pint1=$p2+$parameters{'microsatellite minimum distance'};
				$pint2=$pint1+$size;
			}
			push(@intchr,"$id $pint1 $pint2 $type");
			$maxpos=int($pint2/1000000)+1;
		}
	}
	close CHR;
	$chromint=$chromint."hs$chr:$minpos-$maxpos;";
	push(@int,@intchr);
}
print CIRCOS "chromosomes_display_default = no\nchromosomes = $chromint\nchromosomes_units = 1000000\n";
open(TEMPLATE,"$parameters{'circos templates directory'}/template.txt") ||die "No $parameters{'circos templates directory'}/template.txt file\n";
while (my $line=<TEMPLATE>) {print CIRCOS $line;}
close TEMPLATE;
my $geneTracks="$parameters{'circos data directory'}/$parameters{'circos plot name'}genes.track";
open (TRACK,"> $geneTracks") ||die "Can't write $geneTracks\n";
foreach my $in (@int) {
	my ($id,$p1,$p2,$type)=split(/\s/,$in);
	my ($chr,$pos)=split(/\./,$id);
	my $pos1=$p1+$size/2-20;
	my $pos2=$pos1+40;
	my $geneMS=getMS($id,$type);
	if ($geneMS eq "") {$geneMS="UNK";}
	print TRACK "hs$chr  $pos1  $pos2  $geneMS\n";  
}
close TRACK;

my @sampleList;
my %MSIpercent;
if ($parameters{'sample list file'} ne "") {
	open(XDATA,"$parameters{'sample list file'}") ||die;
	while (my $sample=<XDATA>) {
	    	chomp $sample;
		if ($sample ne 'blank') {
		open(IN, "head -1 $parameters{'output directory'}/$sample.MSIStatus.counts.csv |") ||print STDERR "Skipping sample $sample, MSI Status file not found: $parameters{'output directory'}/$sample.MSIStatus.counts.csv\n";
		while (my $line=<IN>) {
			chomp $line;
			my ($s,$s,$s,$msiPerc,$ratio)=split(/\s/,$line);
			if ($msiPerc =~ '%') {
				my $msiPerc=~s/\%//;
				$MSIpercent{$sample}=$msiPerc;
			}
			else {
				print STDERR "Skipping sample $sample, MSI Status NA\n";
			}
		}
		close IN;
		}
		push(@sampleList,$sample);
	}
	close XDATA;
}
elsif ($parameters{'sample list order'} eq "MSI in" || $parameters{'sample list order'} eq "MSI out") {
	open(XDATA,"head -1 $parameters{'output directory'}/*.MSIStatus.counts.csv | grep 'Status'|") ||die "Unable to find status files in $parameters{'output directory'} directory\n";
	while (my $line = <XDATA>) {
		chomp $line;
		my @path=split(/\//,$line);
		my $statusFile=$path[-1];
		my ($sample,@remainder)=split(/\./,$statusFile);
		my $line=<XDATA>;
		my ($s,$t,$u,$msiPerc,$ratio)=split(/\s/,$line);
		if ($msiPerc =~ '\%') {
			$msiPerc=~s/\%//;
			$MSIpercent{$sample}=int($msiPerc);
		}
		else {
			print STDERR "Skipping sample $sample, MSI Status NA\n";
		}
	}
	close XDATA;
	if ($parameters{'sample list order'} eq "MSI in") {
		@sampleList = sort { $MSIpercent{$a} <=> $MSIpercent{$b} } keys(%MSIpercent);
	}
	else {
		@sampleList = sort { $MSIpercent{$b} <=> $MSIpercent{$a} } keys(%MSIpercent);
	}
}
my $r0=0.9;
my $step=($r0-$parameters{'circos inner circle'})/(2+$#sampleList);
print CIRCOS "<image>
dir = $parameters{'output directory'}\n";
if ($parameters{'circos plot name'} ne "") {
	my $name=substr($parameters{'circos plot name'},0,length($parameters{'circos plot name'})-1);
	print CIRCOS "file* = $name\n";
}
my $radius=30*($#sampleList+1);
if ($radius < 1500) {
	$radius=1500;
}
print CIRCOS "radius* = ${radius}p
</image>\n";
print CIRCOS "<plots>\n
type            = tile
margin      = 0.01u
orientation = out
layers      = 16
thickness = 30
stroke_thickness = 1
stroke_color     = black


<backgrounds>
<background>
color     = white
y0        = 0.0001
</background>
<background>
color = white
y1        = 0.0001
</background>
</backgrounds>

<plot>
file = $geneTracks
type  = text
color = black
r0 = .91r
r1 = .91r+100p
label_size = 10
label_font = condensed
show_links     = yes
link_dims      = 0p,2p,6p,2p,5p
link_thickness = 2p
link_color     = black
label_snuggle        = yes
snuggle_refine       = yes
</plot>\n";
my $prevsample="";
my $blank=-1;
foreach my $sample (@sampleList) {
	my $xdataFile= "$parameters{'circos data directory'}/xdata.$parameters{'circos plot name'}$sample";
	open(OUT,"> $xdataFile") ||die "Can't write $xdataFile\n";
	print STDERR "Writing $xdataFile phenotype file for CIRCOS\n";
	if ($r0 <=0) {
		print STDERR "Too many samples to display, reaching graph center, computing CIRCOS graph until current sample $sample.\n";
		print CIRCOS "</plots>\n";
		close CIRCOS;
		last;
	}
	my $r1=$r0;
	$r0=$r0-$step;
	print CIRCOS "<plot>\nfile = $xdataFile\nr1 = ${r1}r\nr0 = ${r0}r\n<<include $parameters{'circos templates directory'}/rules.conf>>\n</plot>\n";
	$r0=$r0-0.001;
	if ($sample eq "blank") {
                        open(BLANK,"> $parameters{'circos data directory'}/xdata.blank$blank") ||die;
                        foreach my $in (@int) {
                                my ($id,$p1,$p2,$type)=split(/\s/,$in);
				my ($chr,$p)=split(/\./,$id);
                                print BLANK "hs$chr  $p1  $p2  id=blank,pheno=5\n";
                        }
                        close BLANK;
			$r1=$r0;
                        $r0=$r0-$step;
                        print CIRCOS "<plot>\nfile = $parameters{'circos data directory'}/xdata.blank$blank\nr1 = ${r1}r\nr0 = ${r0}r\n<<include $parameters{'circos templates directory'}/rules.conf>>\n</plot>\n";
                       $r0=$r0-0.001;
                        
	}
	elsif ($prevsample ne "" && $parameters{'sample list order'} ne "") {
		if (($MSIpercent{$prevsample} <= $parameters{'MSI-L cutoff'} && $MSIpercent{$sample} > $parameters{'MSI-L cutoff'}) || ($MSIpercent{$prevsample} > $parameters{'MSI-L cutoff'} && $MSIpercent{$sample} <= $parameters{'MSI-L cutoff'})) {
			$blank=0;
		}
		elsif (($MSIpercent{$prevsample} <= $parameters{'MSI-H cutoff'} && $MSIpercent{$sample} > $parameters{'MSI-H cutoff'}) || ($MSIpercent{$prevsample} > $parameters{'MSI-H cutoff'} && $MSIpercent{$sample} <= $parameters{'MSI-H cutoff'})) {
			$blank=1;
}
		if ($blank != -1) {
			open(BLANK,"> $parameters{'circos data directory'}/xdata.blank$blank") ||die; 
			foreach my $in (@int) {
				my ($id,$p1,$p2,$type)=split(/\s/,$in);
				my ($chr,$p)=split(/\./,$id);
                        	print BLANK "hs$chr  $p1  $p2  id=blank,pheno=5\n";
                	}
			close BLANK;
			$r1=$r0;
			$r0=$r0-$step;
			print STDERR "Writing blank to separate MSI-H/MSI-L/MSS threshold\n";
			print CIRCOS "<plot>\nfile = $parameters{'circos data directory'}/xdata.blank$blank\nr1 = ${r1}r\nr0 = ${r0}r\n<<include $parameters{'circos templates directory'}/rules.conf>>\n</plot>\n";
        		$r0=$r0-0.001;
			$blank=-1;
		}
	}
	if ($sample ne 'blank') {
	$prevsample=$sample;
	my $codingFile=1;
	my $skippingFile=1;
	my $intronicFile=1;
	if (! -f "$parameters{'output directory'}/$sample.CODING.counts.csv") {
		$codingFile=0;
		print STDERR "No coding file $parameters{'output directory'}/$sample.CODING.counts.csv for sample, coding phenotype will be NA for all coding/UTR microsatellites\n";
	}
	if (! -f "$parameters{'output directory'}/$sample.SKIPPING.counts.csv") {
                $skippingFile=0;
                print STDERR "No skipping file $parameters{'output directory'}/$sample.SKIPPING.counts.csv for sample, skipping phenotype will be NA for all coding/UTR microsatellites\n";
        }
	if (! -f "$parameters{'output directory'}/$sample.INTRONIC.counts.csv") {
		$intronicFile=0;
		print STDERR "No intronic file $parameters{'output directory'}/$sample.INTRONIC.counts.csv for sample, intronic phenotype will be NA for all INTRONIC microsatellites\n";
	}
	foreach my $in (@int) {
		my $pheno=0;
		my ($id,$p1,$p2,$type)=split(/\s/,$in);
		my ($chr,$pos)=split(/\./,$id);
		if ($type eq "panel" ||$type eq "coding") {
			if ($codingFile!=0) {
				$pheno=getCodingPheno($sample,$id,$type);
			}
		}
		elsif ($type eq "skipping") {
			if ($skippingFile!=0) {
				$pheno=getSkippingPheno($sample,$id,$type);
			}
		}
		elsif ($type eq "intronic") {
			if ($intronicFile!=0) {
				$pheno=getCodingPheno($sample,$id,$type);
			}
		}
		print OUT "hs$chr  $p1  $p2  id=$sample,pheno=$pheno\n";
	}
	}
	close OUT;
}
print CIRCOS "</plots>\n";
close CIRCOS;
print STDERR "Done writing CIRCOS configuration file $parameters{'circos data directory'}/$plotCircos\n";
print STDERR "Executing $parameters{'circos binary'} -conf $plotCircos\n";
system("$parameters{'circos binary'} -conf $plotCircos");
sub getMsType {
	my ($i)=@_;
	my $ret="";
	open(PANEL,"grep -w $i $parameters{'UTR panel'} | wc -l|") ||die "No $parameters{'UTR panel'} found\n";
	my $count=<PANEL>;
	chomp $count;
	close PANEL;
	if ($count > 0) {
		$ret="panel";
	}
	else {
		open(TYPE,"grep -w $i $parameters{'output directory'}/*CODING.counts.csv|wc -l|") ||die "No $parameters{'output directory'}/*.csv files \n";
		my $coding=<TYPE>;
		chomp $coding;
		close TYPE;
		if ($coding > 0) {
			$ret="coding";
		}
		else {
			open(TYPE,"grep -w $i $parameters{'output directory'}/*SKIPPING.counts.csv|wc -l|") ||die "No $parameters{'output directory'}/*.csv files \n";
			my $skipping=<TYPE>;
			chomp $skipping;
			close TYPE;
			if ($skipping > 0) {
				$ret="skipping";

			}
			else {
				open(TYPE,"grep -w $i $parameters{'output directory'}/*INTRONIC.counts.csv|wc -l|") ||die "No $parameters{'output directory'}/*.csv files \n";
				my $intronic=<TYPE>;
				chomp $intronic;
				close TYPE;
				if ($intronic > 0) {
					$ret="intronic";
				}
			}
		
		}
	}
	return $ret;
}

sub getCodingPheno {
	my ($bc,$i,$t)=@_;
	my %panelPheno=("NORMAL"=> 3, "MSI" => 2, "MSI LOW EXPRESSION" => 1, "LOW EXPRESSION" => 0);
	my %codingPheno=("NORMAL"=> 3, "MSI" => 8, "MSI LOW EXPRESSION" => 9, "LOW EXPRESSION" => 0);
	my %intronicPheno=("NORMAL"=>3, "MSS"=>3, "MSI" => 11, "MSI LOW EXPRESSION" => 12, "LOW EXPRESSION" => 0);
	my $ph=0;
	my $file="INTRONIC";
	if ($t eq "coding" ||$t eq "panel") { 
		$file="CODING";
	}
	open(GREP,"grep -w $i $parameters{'output directory'}/$bc.$file.counts.csv|cut -f5|") ||die "Can't read $parameters{'output directory'}/$bc.CODING.counts.csv\n";
	while (my $p=<GREP>) {
		chomp $p;
		if ($t eq "coding") {
			$ph=$codingPheno{$p};
		}
		if ($t eq "panel") {
			$ph=$panelPheno{$p};
		}
		if ($t eq "intronic") {
			$ph=$intronicPheno{$p};
		}
		if ($ph eq "") {$ph=0;}
	}
	close GREP;
	return $ph;
}

sub getSkippingPheno {
	my ($bc,$i)=@_;
	my $ph=0;
	my %skippingPheno=("NORMAL"=> 3, "MSI" => 6, "MSI LOW EXPRESSION" => 10, "LOW EXPRESSION" => 0);
	open(GREP,"grep -w $i $parameters{'output directory'}/$bc.SKIPPING.counts.csv|cut -f12|") ||die "Can't read $parameters{'output directory'}/$bc.SKIPPING.counts.csv file\n";
	while (my $p=<GREP>) {
		chomp $p;
		$ph=$skippingPheno{$p};
	}
	close GREP;
	return $ph;
}
sub getMS {
	my ($i,$t)=@_;
	my $ms="";
	if ($t eq "panel" ||$t eq "coding") {
		open(GREP,"grep -w $i $parameters{'output directory'}/*CODING.counts.csv|head -1|cut -f2,3,4|") ||die "Can't read $parameters{'output directory'}/*CODING.counts.csv\n";
		while (my $ln=<GREP>) {
			chomp $ln;
			my ($base,$length,$desc)=split(/\t/,$ln);
			my ($n,$gn,@reste)=split(/\-/,$desc);
			$ms="$gn-($base$length)";
		}
		close GREP;
		return $ms;
	}
	if ($t eq "skipping") {
		open(GREP,"grep -w $i $parameters{'output directory'}/*SKIPPING.counts.csv|head -1|cut -f2,3,4|") ||die "Can't read $parameters{'output directory'}/*SKIPPING.counts.csv\n";
		while (my $ln=<GREP>) {
			chomp $ln;
			my ($base,$length,$desc)=split(/\t/,$ln);
			$ms="$desc-($base$length)";
		}
		close GREP;
		return $ms;
	}
	if ($t eq "intronic") {
		open(GREP,"grep -w $i $parameters{'output directory'}/*INTRONIC.counts.csv|head -1|cut -f2,3,4|") ||die "Can't read $parameters{'output directory'}/*INTRONIC.counts.csv\n";
                while (my $ln=<GREP>) {
                        chomp $ln;
                        my ($base,$length,$desc)=split(/\t/,$ln);
                        my ($n,$gn,@reste)=split(/\-/,$desc);
                        $ms="$gn-($base$length)";
                }
                close GREP;
                return $ms;
	}
}
sub checkParameters {
	my ($config)=@_;
	my %params;
	open(CONFIG,$config) ||die "Can't read configuration file $config\n";
	while (my $line=<CONFIG>) {
		if ($line=~ "=") {
			chomp $line;
			my ($key,$value)=split(/\=/,$line);
			$params{$key}=$value;
		}
	}
	close CONFIG;
	if (! -f $params{'circos microsatellite list'}) {
		die $params{'circos microsatellite list'}." file not found\n";
	}
	if (! -d $params{'circos data directory'}){
		print STDERR "Creating data directory $params{'circos data directory'}\n";
		mkdir ($params{'circos data directory'});
	}
	if (! -d $params{'circos templates directory'}) {
		die "CIRCOS templates directory $params{'circos templates directory'} not found.\n";
	}
	if ($params{'sample list order'} eq "file"){
		if (! -f $params{'sample list file'}) {
			die "List of samples $params{'sample list file'} not found\n";
		}
	}
	if (! -d $params{'output directory'}) {
		die "MSInspector data directory $params{'output directory'} not found\n";
	}
	if ($params{'microsatellite minimum distance'} eq "") {
		$params{'microsatellite minimum distance'}=500000;
	}
	else {
		$params{'microsatellite minimum distance'}=$params{'microsatellite minimum distance'}*1000000;
	}
	if ($params{'microsatellite display size'} eq "") {
		$params{'microsatellite display size'}=4000000;
	}
	else {
		$params{'microsatellite display size'}=$params{'microsatellite display size'}*1000000;
	}
	if (! -f $params{'UTR panel'}) {
		die "UTR panel $params{'UTR panel'} file not found\n";
	}
	if ($params{'MSI high cutoff'} eq "") {
		$params{'MSI high cutoff'}=40;
	}
	if ($params{'MSI low cutoff'} eq "") {
		$params{'MSI low cutoff'}=30;
	}
	if ($params{'circos binary'} eq "") {
		open(WHICH,"which circos|") ||die "Unable to find circos binary\n";
		my $circos=<WHICH>;
		chomp $circos;
		close WHICH;
		if ($circos eq "") { die "Unable to test circos binary\n";}
		else {$params{'circos binary'}=$circos;}
	} elsif (! -f $params{'circos binary'}) {
		die "CIRCOS binary $params{'circos binary'} not found\n";
	}
	if ($params{'circos inner circle'} eq "") {
		$params{'circos inner circle'}=.3;
		print STDERR "Circos inner circle is not defined, suing default value of .3\n";
	}
	if ($params{'circos plot name'} ne "") { $params{'circos plot name'}=$params{'circos plot name'}.".";}
	return %params;
}
1;
