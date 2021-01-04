#!/usr/bin/perl -w

open (f1, "peptide.fasta"); #Input: Peptide sequence file in ".fasta" format. The header should only have the gene name, no gene description. eg >YAL007:Mitochondrial gene regulating the water level(Not allowed).>YAL007(Allowed)
@line1=<f1>; chomp(@line1);
close(f1);
$n1 = @line1;
$l1 = ($n1 - 1);

open (f2, "codontable.txt"); #Input: The codon table. You can input the desired number of synonymous codons for corresponding amino-acids. See SAMPLE file attached (codontable.txt).
@line2=<f2>; chomp(@line2);
close(f2);
$n2 = @line2;
$l2 = ($n2 - 1);

@seq = ();

open (f3, ">Randomized_CDS.fasta"); #Output: Name of the (nucleic-acid) output file having all the sequences as in the input peptide file.


foreach $i(0..$l1)
	{
	$j = ($i + 1);
	if ($line1[$i] =~ /^>/)
		{
		print f3 "$line1[$i]\n";
		print "$line1[$i]\n";
		}
		
	elsif ($line1[$i] !~ /^>/ && $line1[$j] =~ /^>/ || $i == $l1)
		{
		push (@seq, $line1[$i]);
		$fullseq = join("", @seq);
		@seq_amino = split(//, $fullseq);
		
		foreach $_(@seq_amino)
			{
			foreach $k(0..$l2)
				{
				@new1 = split(/ = /, $line2[$k]); 
				if ("$_" eq "$new1[0]")
					{
					@new2 = split(/, /, $new1[1]);
					$n3 = @new2;
					$r1 = int(rand($n3));
					print f3 "$new2[$r1]";
					}
				else
					{
					}
				}
			}
		@end1 = split(/ = /, $line2[$l2]);
		@end2 = split(/, /, $end1[1]);
		$n4 = @end2;
		$r2 = int(rand($n4));
		print f3 "$end2[$r2]\n";
		@seq = ();
		}
		
	else
		{
		push (@seq, $line1[$i]);
		}
	}