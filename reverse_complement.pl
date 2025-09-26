#!usr/bin/env perl
# Doing reverse complement based on the results of TransDecoder or BLASTP hits to Swiss-Prot.
# Usage: perl reverse_complement.pl lpi.okay.fasta.transdecoder_dir/transdecoder_direction.list lpi_filtered_long.fasta > lpi_filtered_long_plus.fasta
# reverse_complement.pl
# YJ Luo | 2019.06.10

use strict;
use warnings;

open (LIST, $ARGV[0]);
open (IN, $ARGV[1]);

my %hash_list = ();
while (<LIST>) {
        chomp;
        my @fields = split (" ", $_);
        $hash_list{$fields[0]} = $fields[1];
}

local   $/ = '>';
while (<IN>) {
        chomp;
    if (/(\S+_\S+\S+)/) {
                my $id = $1;
                $_ =~ /$id/p;
                # Look at http://perldoc.perl.org/perlreref.html for explanation of using ${^POSTMATCH}
        my $seq = ${^POSTMATCH};
        $seq =~ s/\n//g;
                if (exists $hash_list{$id} && $hash_list{$id} eq '-') {
                        print '>' . $id . "\n";
                        print reverse_complement($seq) . "\n";
                } else {
                        print '>' . $id . "\n";
                        print $seq . "\n";
                }
        }
}

sub reverse_complement {
my $dna = shift;
# reverse the DNA sequence
my $revcomp = reverse($dna);
# complement the reversed DNA sequence
$revcomp =~ tr/ACGTacgt/TGCAtgca/;
return $revcomp;
}
#--//