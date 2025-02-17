###################
#
#  gblock.pl
#
#  Converts a fasta format text input to a list of DNA gblocks with the pre-defined T7 5' adaptor sequence.
#  It converts the RNA sequence to DNA, then takes the reverse compliment, and adds the primer adaptor to the 3' end.
#
##################

use strict;
use warnings;

# Define the adaptor sequence
my $adaptor = "TATAGTGAGTCGTATTAGAATGCCC";

# Compute reverse complement of the adaptor sequence
my $rev_adaptor = reverse $adaptor;
$rev_adaptor =~ tr/ACGT/TGCA/;

# Check for correct usage
if (@ARGV != 1) {
    die "Usage: perl gblock.pl input.txt > output.txt\n";
}

# Input file name from command line argument
my $input_file = $ARGV[0];

# Open input file
open(my $in,  "<", $input_file) or die "Cannot open input file: $!";

# Store RNA sequences and DNA sequences separately
my @rna_sequences;
my @dna_sequences;

# Print the reverse complement of the adaptor as the first sequence
print ">T7_5strong.FWD\n$rev_adaptor\n\n";

my $header = "";
my $sequence = "";
my $expect_sequence = 0;  # Flag to track if the next valid line is the sequence

while (my $line = <$in>) {
    chomp $line;
    next if ($line =~ /^#/ || $line =~ /^[().\[\]]+$/); # Ignore comments and dot-bracket structures
    
    if ($line =~ /^>/) {
        # Process previous sequence if exists
        if ($header ne "" && $sequence ne "") {
            push @rna_sequences, "$header\n$sequence";
            my $processed_seq = process_sequence($sequence, $adaptor);
            push @dna_sequences, "$header\n$processed_seq";
        }
        $header = $line;  # Store new header
        $sequence = "";  # Reset sequence
        $expect_sequence = 1;  # Set flag to expect sequence next
    } elsif ($expect_sequence) {
        next if ($line =~ /[^ACGUTNacgutn]/); # Ignore invalid sequence lines
        $sequence = $line;  # Store only the first sequence line after header
        $expect_sequence = 0;  # Reset flag after reading one sequence line
    }
}

# Process the last sequence in the file
if ($header ne "" && $sequence ne "") {
    push @rna_sequences, "$header\n$sequence";
    my $processed_seq = process_sequence($sequence, $adaptor);
    push @dna_sequences, "$header\n$processed_seq";
}

close $in;

# Print RNA sequences
print "# RNA Sequences\n";
print join("\n", @rna_sequences), "\n\n";

# Print separator
print "# Processed DNA Sequences\n";

# Print DNA sequences
print join("\n", @dna_sequences), "\n";

print "Processing complete. Output written to standard output\n";

sub process_sequence {
    my ($seq, $adaptor) = @_;
    
    # Convert RNA to DNA (U -> T)
    $seq =~ tr/U/T/;
    
    # Compute reverse complement
    my $rev_complement = reverse $seq;
    $rev_complement =~ tr/ACGT/TGCA/;
    
    # Append adaptor sequence
    return $rev_complement . $adaptor;
}