#! /usr/bin/perl

use strict;
use warnings;

my $argn = $#ARGV+1;
if($argn != 3){
        print "\n3_stat_alignment.pl aln.phy window_size step_size\n\n";
        exit(0);
}

my @msa;
my $window = $ARGV[1];
my $step = $ARGV[2];
my $position = $window/2;
my %pval;
my $aln_file = $ARGV[0];

open CSV, ">CSV";

print CSV "Position\tNSS\tChi\tPHI\n";

while ($position < 1044459) {

	&extract($position, $window);
	&run_phi($position);
	print "Position: ".$position."\n";
	print CSV $position."\t".$pval{$position}{NSS}."\t".$pval{$position}{Chi}."\t".$pval{$position}{Phi}."\n";
	$position = $position+$step;

}

close CSV;




############################# SUBROUTINES ##########################



############## extract from aln ########

sub extract
{
	my ($pos, $w_size) = @_;
	my $from = ($pos-($w_size/2))+1;
	my $to = $w_size;
	open MSA, "$aln_file" or die"Can't open MSA";
	@msa = <MSA>;
	close MSA;
	
	open (T, ">temp");
	my $i=1;
	while ($i < $#msa+1) {
		my @a = split(/\s+/,$msa[$i]);
		print T ">".$a[0]."\n".substr($a[1], $from, $to)."\n";
		$i = $i+1;
	}
	close T;
}


######## Running Phi #########

sub run_phi
{
	my $pos = $_[0];
	open PIPE, "/home/anna/bijagosgenomes2012/fastq95/PhiPack/Phi -f temp -o -p 1000 | tail -5 | head -3 |";
	my @LaPhi = <PIPE>;
		my @NSS = split(/\s+/, $LaPhi[0]);
		$pval{$pos}{NSS} = $NSS[1];
		my @Chi = split(/\s+/, $LaPhi[1]);
		$pval{$pos}{Chi} = $Chi[2];
		my @PHI = split(/\s+/, $LaPhi[2]);
		$pval{$pos}{Phi} = $PHI[2];
	close PIPE;
}

