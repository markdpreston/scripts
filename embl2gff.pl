#!/usr/bin/perl -w
# Take EMBL file(s) and make glimmerHMM training files
#

use strict;
use Bio::SeqIO;

#my $id_type = 'locus_tag';
my $id_type = 'systematic_id';

# Pass series of files or union file on command line
my @files = @ARGV;

my $exon_counter = 1;

foreach my $filename (@files)
{
	my $contig_name = $filename;

	open(F, "<$filename") or die "$!";
	
	while(<F>)
	{
		if(/ID\s+(\S+)/)
		{
			$contig_name = $1;
			last;
		}
	}
	close F;
	
	my $stream = Bio::SeqIO->new(-file => $filename, -format => 'EMBL');

	$filename =~ s/\.embl//;

	while ( (my $seq = $stream->next_seq()) ) 
	{
		my @features = $seq->get_SeqFeatures();

		my $current_id = '';
		my $current_locus_tag = '';
		my $current_string = 0;

		foreach my $feat (@features)
		{
			# Ignore unless CDS
			next unless $feat->primary_tag eq 'CDS';

			$current_id = join(' ',$feat->get_tag_values($id_type));

			my $current_strand = $feat->strand;

			my $feat_string = $feat->location->to_FTstring();

			my @exons = split /,/, $feat_string;

			my %exons = ();

			my $gene_start = $feat->start;
			my $gene_stop = $feat->end;

			my $str = $current_strand == 1 ? '+' : '-';

			# gene_id "CUFF.10"; transcript_id "CUFF.10.1"; exon_number "1"

			#print "$filename\tDEFAULT\tgene\t$gene_start\t$gene_stop\t\.\t$str\t\.\tID=Gene:$current_id\;\n";
			print "$contig_name\tDEFAULT\ttranscript\t$gene_start\t$gene_stop\t\.\t$str\t\.\tgene_id \"$current_id\"; transcript_id \"$current_id\.1\";\n";

			foreach my $exon (@exons)
			{
				$exon =~ /(\d+)\.\.(\d+)/;

				next if !defined $2;

				print "$contig_name\tDEFAULT\texon\t$1\t$2\t\.\t$str\t\.\tgene_id \"$current_id\"; transcript_id \"$current_id\.1\"; exon_number \"$exon_counter\";\n";

				$exon_counter++;

				$exons{$1} = $2;
			}

			my @sorted_exons = ();

			if($current_strand == 1)
			{
				@sorted_exons = sort {$a<=>$b} keys %exons;
			}
			else
			{
				@sorted_exons = sort {$b<=>$a} keys %exons;
			}

			foreach my $start (@sorted_exons)
			{
			}
		}
	}
}
