#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

# Input FASTA file
my $fasta_file = shift @ARGV or die "Usage: $0 <input.fasta>\n";

# Parameters
my $bulge_penalty = 3;     # Penalty for a bulge in Tm calculation
my $min_tm = 40;           # Tm threshold
my $min_complementary = 6; # Minimum number of complementary bases to report
my $context_window = 5;    # Number of neighboring bases for context
my $top_results = 20;      # Number of top results to display for each group

# Parse FASTA
my @sequences = parse_fasta($fasta_file);
die "FASTA file must contain exactly 3 sequences\n" unless @sequences == 3;

my ($primer1, $primer2, $template) = @sequences;

# Check for RNA sequences
foreach my $seq (@sequences) {
    if ($seq->{seq} =~ /U/) {
        warn "Warning: Sequence $seq->{name} contains 'U' - possibly RNA.\n";
    }
    if ($seq->{seq} =~ /[^ATCG]/i) {
        die "Error: Sequence $seq->{name} contains invalid characters.\n";
    }
}

# Compute reverse complement of the template and reverse sequences
$template->{revcom} = reverse_complement($template->{seq});
$primer1->{rev} = reverse($primer1->{seq});
$primer2->{rev} = reverse($primer2->{seq});
$template->{rev} = reverse($template->{seq}); # Reverse the template as well

# Generate variants with a single bulge for each primer
my @primer1_variants = generate_bulged_variants($primer1->{seq});
my @primer2_variants = generate_bulged_variants($primer2->{seq});

# Debug: Parsed sequences and bulged variants
print "Debug: Parsed sequences:\n", Dumper(\@sequences);
print "Debug: Primer1 Variants:\n", Dumper(\@primer1_variants);
print "Debug: Primer2 Variants:\n", Dumper(\@primer2_variants);

# Create a smart output file name based on the input file name
my ($input_base) = $fasta_file =~ /([^\/\\]+)\.[^.]+$/;  # Extract the base name (excluding path and extension)
my $output_file = "${input_base}_output.txt";

# Open the output file
open my $out, '>', $output_file or die "Could not write to $output_file: $!\n";
print "Output will be saved to: $output_file\n";

# Write template and primers to output file
print $out ">Template\n$template->{seq}\n";
print $out ">Reverse_Complement\n$template->{revcom}\n";
print $out ">$primer1->{name}\n$primer1->{seq}\n";
print $out ">Reverse_$primer1->{name}\n$primer1->{rev}\n";
print $out ">$primer2->{name}\n$primer2->{seq}\n";
print $out ">Reverse_$primer2->{name}\n$primer2->{rev}\n";

# Initialize annotation arrays for the dsDNA
my $template_annotation = ' ' x length($template->{seq});    # Top strand annotation
my $revcom_annotation = ' ' x length($template->{revcom});  # Bottom strand annotation

# Define pairings for analysis
my @pairings = (
    [$template->{seq},         \@primer1_variants, $template->{name}, $primer1->{name}, "Template vs Forward"],
    [$template->{revcom},      \@primer1_variants, "Reverse Complement of $template->{name}", $primer1->{name}, "Reverse Complement vs Forward"],
    [$template->{seq},         \@primer2_variants, $template->{name}, $primer2->{name}, "Template vs Reverse"],
    [$template->{revcom},      \@primer2_variants, "Reverse Complement of $template->{name}", $primer2->{name}, "Reverse Complement vs Reverse"],
    [$primer1->{seq},          \@primer2_variants, $primer1->{name}, $primer2->{name}, "Fwd vs Reverse"],
    [$primer1->{seq},          \@primer1_variants, $primer1->{name}, "$primer1->{name} Self", "Fwd Self"],
    [$primer2->{seq},          \@primer2_variants, $primer2->{name}, "$primer2->{name} Self", "Rev Self"],
);

# Perform analysis
foreach my $pair (@pairings) {
    my ($template_seq, $primer_variants, $template_name, $primer_name, $label) = @$pair;
    print "Debug: Analyzing $label\n";

    my @results;
    foreach my $primer_seq (@$primer_variants) {
        # Reverse the primer sequence before testing
        my $reversed_primer_seq = reverse($primer_seq);
        push @results, find_complementary_matches($template_seq, $reversed_primer_seq);
    }

    @results = sort { $b->{tm} <=> $a->{tm} } @results;  # Sort by Tm descending
    splice(@results, $top_results) if @results > $top_results;  # Keep only top results

    # Annotate binding sites only for relevant pairings
    if ($label =~ /Template vs Forward/) {
        foreach my $result (@results) {
            annotate_binding_sites($result, \$template_annotation, 'F');  # Forward primer on top strand
        }
    }
    elsif ($label =~ /Reverse Complement vs Reverse/) {
        foreach my $result (@results) {
            annotate_binding_sites($result, \$revcom_annotation, 'R');  # Reverse primer on bottom strand
        }
    }
    elsif ($label =~ /Reverse Complement vs Forward/) {
        foreach my $result (@results) {
            annotate_binding_sites($result, \$revcom_annotation, 'F');  # Forward primer on bottom strand
        }
    }
    elsif ($label =~ /Template vs Reverse/) {
        foreach my $result (@results) {
            annotate_binding_sites($result, \$template_annotation, 'R');  # Reverse primer on top strand
        }
    }

    # Write analysis results
    print $out "\n# $label (Top $top_results Tm)\n";
    print $out "# Specific Sequences: $template_name, $primer_name\n\n";

    foreach my $result (@results) {
        # Add 5- and -3 to the sequences
        my $template_with_direction = "5-" . $result->{template_context} . "-3";
        my $primer_with_direction = "3-" . $result->{primer_context} . "-5";

        # Calculate the padding for alignment
        my $alignment_padding = ' ' x 2;  # Indentation
        my $template_padding = ' ' x 2;   # Indentation for template
        my $primer_padding = ' ' x 2;

        # Adjust the alignment to account for the `5-` marker
        my $adjusted_alignment = (' ' x 2) . (' ' x 2) . $result->{alignment};

        # Ensure proper padding between markers and sequences
        $template_with_direction = "5-" . $result->{template_context} . "-3";
        $primer_with_direction = "3-" . $result->{primer_context} . "-5";

        print $out "Tm:$result->{tm}_Offset:$result->{offset}\n";
        print $out "$template_padding$template_with_direction\n";
        print $out "$adjusted_alignment\n";
        print $out "$primer_padding$primer_with_direction\n\n";
    }
}

# Annotate the sequence with binding sites
sub annotate_binding_sites {
    my ($result, $annotation_ref, $marker) = @_;
    my $offset = $result->{offset};
    my $alignment = $result->{alignment};
    my $primer_len = length($alignment);  # Alignment length matches the primer length

    for my $i (0 .. $primer_len - 1) {
        if (substr($alignment, $i, 1) eq ':') {
            # Only annotate positions with base-pairing
            my $current_char = substr($$annotation_ref, $offset + $i, 1);
            if ($current_char eq ' ') {
                substr($$annotation_ref, $offset + $i, 1) = $marker;
            } elsif (($current_char ne $marker)) {
                substr($$annotation_ref, $offset + $i, 1) = 'X';  # Overlap as 'X'
            }
        }
    }
}

# Reverse a string while preserving spacing alignment
sub reverse_with_padding {
    my ($string) = @_;
    my $reversed = reverse($string);
    return $reversed;
}

# Output the annotated dsDNA
sub output_annotated_dsDNA {
    my ($template, $template_annotation, $revcom_annotation) = @_;

    # Reverse the reverse complement and its annotation for proper alignment
    my $reversed_revcom = reverse($template->{revcom});
    my $reversed_annotation = reverse_with_padding($revcom_annotation);

    # Create the base-pairing line
    my $base_pairing = '';
    for my $i (0 .. length($template->{seq}) - 1) {
        $base_pairing .= ':';
    }

    print $out "\n# Annotated dsDNA\n";
    print $out "$template_annotation\n";
    print $out "   5-" . $template->{seq} . "-3\n";
    print $out "     $base_pairing\n";
    print $out "   3-" . $reversed_revcom . "-5\n";
    print $out "     $reversed_annotation\n";
}

# After performing analysis, output the annotated dsDNA
output_annotated_dsDNA($template, $template_annotation, $revcom_annotation);

close $out;

# --- Subroutines ---

# Generate bulged variants of a sequence
sub generate_bulged_variants {
    my ($seq) = @_;
    my @variants;

    # Add the original sequence as the first variant
    push @variants, $seq;

    # Generate variants with a bulge
    for my $i (1 .. length($seq) - 1) { # Start at position 1, ignore end-of-sequence bulge
        my $variant = substr($seq, 0, $i) . '-' . substr($seq, $i);
        push @variants, $variant;
    }

    return @variants;
}


# Parse FASTA file
sub parse_fasta {
    my ($file) = @_;
    my @seqs;
    open my $fh, '<', $file or die "Could not open $file: $!\n";
    my ($name, $seq);
    while (<$fh>) {
        chomp;
        if (/^>(.*)/) {
            if (defined $name) {
                push @seqs, { name => $name, seq => uc $seq };
            }
            $name = $1;
            $seq = '';
        } else {
            $seq .= uc $_;
        }
    }
    push @seqs, { name => $name, seq => uc $seq } if defined $name;
    close $fh;
    return @seqs;
}

# Reverse complement
sub reverse_complement {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ATCG/TAGC/;
    return $seq;
}

# Reverse sequence (no complement)
sub reverse {
    my ($seq) = @_;
    return scalar reverse($seq);
}

# Find complementary matches
sub find_complementary_matches {
    my ($template_seq, $primer_seq) = @_;
    my $primer_len = length($primer_seq);
    my $template_len = length($template_seq);

    my @matches;

    for my $offset (0 .. $template_len - $primer_len) {
        # Get the aligned part of the template
        my $aligned_template = substr($template_seq, $offset, $primer_len);
        
        # Compute alignment and complementary stats
        my ($alignment, $complementary_count, $tm) = calculate_complementary($primer_seq, $aligned_template);

        if ($complementary_count >= $min_complementary && $tm >= $min_tm) {
            # Compute the context before and after the aligned region
            my $before_context = '';
            my $after_context = '';
            if ($offset >= $context_window) {
                $before_context = substr($template_seq, $offset - $context_window, $context_window);
            } else {
                $before_context = (' ' x ($context_window - $offset)) . substr($template_seq, 0, $offset);
            }

            if ($offset + $primer_len + $context_window <= $template_len) {
                $after_context = substr($template_seq, $offset + $primer_len, $context_window);
            } else {
                $after_context = substr($template_seq, $offset + $primer_len);
            }

            # Construct the full template and primer contexts
            my $template_context = $before_context . $aligned_template . $after_context;
            my $primer_context = (' ' x length($before_context)) . $primer_seq . (' ' x length($after_context));

            # Align the alignment string with the same padding
            my $aligned_symbols = (' ' x length($before_context)) . $alignment . (' ' x length($after_context));

            # Add the result to the matches
            push @matches, {
                tm                => $tm,
                offset            => $offset,
                alignment         => $aligned_symbols,
                template_context  => $template_context,
                primer_context    => $primer_context,
            };
        }
    }

    return @matches;
}

# Calculate complementarity and Tm
sub calculate_complementary {
    my ($primer_seq, $template_seq) = @_;
    my $alignment = '';
    my $complementary_count = 0;
    my $tm = 0;

    for my $i (0 .. length($primer_seq) - 1) {
        my $p_base = substr($primer_seq, $i, 1);
        my $t_base = substr($template_seq, $i, 1);

        if ($p_base eq '-' || $t_base eq '-') {
            # Treat '-' as a bulge
            $alignment .= ' ';
        } elsif (is_complementary($p_base, $t_base)) {
            $alignment .= ':';
            $complementary_count++;
            $tm += ($p_base eq 'G' || $p_base eq 'C') ? 4 : 2;
        } else {
            $alignment .= ' ';
        }
    }

    # Apply bulge penalty if '-' exists in the primer
    if ($primer_seq =~ /-/) {
        $tm -= $bulge_penalty;
    }

    return ($alignment, $complementary_count, $tm);
}

# Check if bases are complementary
sub is_complementary {
    my ($b1, $b2) = @_;
    return ($b1 eq 'A' && $b2 eq 'T') || ($b1 eq 'T' && $b2 eq 'A') ||
           ($b1 eq 'C' && $b2 eq 'G') || ($b1 eq 'G' && $b2 eq 'C');
}

