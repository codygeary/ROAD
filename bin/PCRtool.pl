#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

# Input FASTA file
my $fasta_file = shift @ARGV or die "Usage: $0 <input.fasta>\n";

# Parameters
my $min_tm = 4;           # Tm threshold
my $min_complementary = 4; # Minimum number of complementary bases to report
my $top_results = 10;      # Number of top results to display for each group

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

# # Debug: Parsed sequences and bulged variants
# print "Debug: Parsed sequences:\n", Dumper(\@sequences);
# print "Debug: Primer1 Variants:\n", Dumper(\@primer1_variants);
# print "Debug: Primer2 Variants:\n", Dumper(\@primer2_variants);

# Create a smart output file name based on the input file name
my ($input_base) = $fasta_file =~ /([^\/\\]+)\.[^.]+$/;  # Extract the base name (excluding path and extension)
my $output_file = "${input_base}_output.txt";

# Open the output file
open my $out, '>', $output_file or die "Could not write to $output_file: $!\n";
print "Output will be saved to: $output_file\n";

# Write template and primers to output file
print $out ">Primer1: $primer1->{name}\n$primer1->{seq}\n";
print $out ">Reverse_$primer1->{name}\n$primer1->{rev}\n\n";
print $out ">Primer2: $primer2->{name}\n$primer2->{seq}\n";
print $out ">Reverse_$primer2->{name}\n$primer2->{rev}\n\n";
print $out ">Template: $template->{name}\n$template->{seq}\n";
print $out ">Reverse_Complement\n$template->{revcom}\n\n";

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

    # Calculate Tm for each result
    foreach my $result (@results) {
        my $template_context = $result->{template_context};
        my $primer_context = reverse_complement($result->{primer_context});  # Reverse complement to original direction
        my $alignment = $result->{alignment};
        $result->{tm} = calculate_tm($template_context, $primer_context, $alignment);  # Compute Tm
    }
    
    # Sort results by Tm descending, then by offset ascending
    @results = sort { $b->{tm} <=> $a->{tm} || $a->{offset} <=> $b->{offset} } @results;
    
    my @filtered_results;
    my @offsets;  # Keep track of offsets in parallel with @filtered_results
    
    foreach my $result (@results) {
        my $offset = $result->{offset};
        # Binary search to find the insertion point
        my ($left, $right) = (0, scalar @offsets);
        while ($left < $right) {
            my $mid = int(($left + $right) / 2);
            if ($offsets[$mid] < $offset) {
                $left = $mid + 1;
            } else {
                $right = $mid;
            }
        }
        # Check for nearby offsets
        my $has_nearby = 0;
        if ($left < @offsets && abs($offsets[$left] - $offset) <= 1) {
            $has_nearby = 1;
            # Replace if the new Tm is higher
            if ($result->{tm} > $filtered_results[$left]->{tm}) {
                $filtered_results[$left] = $result;
                $offsets[$left] = $offset;
            }
        } elsif ($left > 0 && abs($offsets[$left - 1] - $offset) <= 1) {
            $has_nearby = 1;
            if ($result->{tm} > $filtered_results[$left - 1]->{tm}) {
                $filtered_results[$left - 1] = $result;
                $offsets[$left - 1] = $offset;
            }
        }
        # Add to filtered results if no nearby offset
        unless ($has_nearby) {
            splice @offsets, $left, 0, $offset;
            splice @filtered_results, $left, 0, $result;
        }
    }
    
    @results = @filtered_results;

    # Filter out results below the minimum Tm threshold
    @results = grep { $_->{tm} >= $min_tm } @results;

    # Sort filtered results by Tm (descending)
    @results = sort { $b->{tm} <=> $a->{tm} } @results;

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
    print $out "\n# $label\n";
    print $out "# Specific Sequences: $template_name, $primer_name\n\n";

    foreach my $result (@results) {
        # Extract template and primer contexts for Tm calculation
        my $template_context = $result->{template_context};
        my $primer_context = reverse_complement($result->{primer_context});
        my $alignment = $result->{alignment};

        # Retrieve already computed Tm
        my $tm = $result->{tm};

        # Add 5- and -3 to the sequences
        my $template_with_direction = "5-" . $template_context . "-3";
        my $primer_with_direction = "3-" . $result->{primer_context} . "-5";

        # Calculate the padding for alignment
        my $alignment_padding = ' ' x 2;  # Indentation
        my $template_padding = ' ' x 2;   # Indentation for template
        my $primer_padding = ' ' x 2;

        # Adjust the alignment to account for the `5-` marker
        my $adjusted_alignment = (' ' x 2) . (' ' x 2) . $alignment;

        # Format and output the results
        my $formatted_tm = sprintf("%.1f", $tm);  # Round Tm to 1 decimal place
        print $out "Tm: $formatted_tm °C   Offset: $result->{offset}\n";
        print $out "$template_padding$template_with_direction\n";
        print $out "$adjusted_alignment\n";
        print $out "$primer_padding$primer_with_direction\n\n";
    }
}


# After performing analysis, output the annotated dsDNA
output_annotated_dsDNA($template, $template_annotation, $revcom_annotation);

close $out;

### --- Subroutines ---

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
            if ($offset >= 0) {
                $before_context = substr($template_seq, $offset, 0);
            } else {
                $before_context = (' ' x (0 - $offset)) . substr($template_seq, 0, $offset);
            }

            if ($offset + $primer_len <= $template_len) {
                $after_context = substr($template_seq, $offset + $primer_len, 0);
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

    return ($alignment, $complementary_count, $tm);
}

# Check if bases are complementary
sub is_complementary {
    my ($b1, $b2) = @_;
    return ($b1 eq 'A' && $b2 eq 'T') || ($b1 eq 'T' && $b2 eq 'A') ||
           ($b1 eq 'C' && $b2 eq 'G') || ($b1 eq 'G' && $b2 eq 'C');
}

# Function to calculate Tm with SantaLucia1996 parameters
sub calculate_tm {  
    # Nearest-neighbor parameters (enthalpy in kcal/mol, entropy in cal/(mol·K))
    my %nn_params = (
        "AA" => [-8.4, -23.6], "TT" => [-8.4, -23.6],
        "AT" => [-6.5, -18.8], "TA" => [-6.3, -18.5],
        "CA" => [-7.4, -19.3], "TG" => [-7.4, -19.3],
        "GT" => [-8.6, -23.0], "AC" => [-8.6, -23.0],
        "CT" => [-6.1, -16.1], "GA" => [-6.1, -16.1],
        "GA" => [-7.7, -20.3], "CT" => [-7.7, -20.3],
        "CG" => [-10.1, -25.5], "GC" => [-11.1, -28.4],
        "GG" => [-6.7, -15.6], "CC" => [-6.7, -15.6],
    );
    
    # Penalties and corrections for initiation and terminal effects
    my %penalties = (
        "initiation_AT" => [0, -9.0],  # A•T initiation enthalpy and entropy
        "initiation_GC" => [0, -5.9],  # G•C initiation enthalpy and entropy
        "symmetry"      => [0, -1.4],  # Symmetry correction
        "terminal_TA"   => [0.4, 0.0], # Penalty for terminal T•A
    );

    my ($template_seq, $primer_seq, $alignment) = @_;

    # Initialize enthalpy and entropy
    my ($total_dh, $total_ds) = (0, 0);

    # Add initiation penalty based on GC content
    if ($template_seq =~ /G|C/) {
        $total_dh += $penalties{"initiation_GC"}->[0];
        $total_ds += $penalties{"initiation_GC"}->[1];
    } else {
        $total_dh += $penalties{"initiation_AT"}->[0];
        $total_ds += $penalties{"initiation_AT"}->[1];
    }

    # Evaluate each nearest-neighbor pair in the alignment
    for my $i (0 .. length($template_seq) - 2) {
        my $template_pair = substr($template_seq, $i, 2);
        my $primer_pair = substr($primer_seq, $i, 2);
        my $align_char = substr($alignment, $i, 1);

        if ($align_char eq ':') {
            # Perfect match
            if (exists $nn_params{$template_pair}) {
                $total_dh += $nn_params{$template_pair}->[0];
                $total_ds += $nn_params{$template_pair}->[1];
            }
        } elsif ($align_char eq '-') {
            # Bulge penalty
            $total_dh += 0.2;
            $total_ds += -5.7;
        } else {
            # Mismatch penalty (not provided in this dataset, approximate with bulge penalty)
            $total_dh += 0.2;
            $total_ds += -5.7;
        }
    }

    # Add terminal penalties (for example, T•A)
    if ($template_seq =~ /^T|A$/) {
        $total_dh += $penalties{"terminal_TA"}->[0];
        $total_ds += $penalties{"terminal_TA"}->[1];
    }

    # Add symmetry correction for self-complementary sequences
    if ($template_seq eq reverse_complement($template_seq)) {
        $total_dh += $penalties{"symmetry"}->[0];
        $total_ds += $penalties{"symmetry"}->[1];
    }

    # Calculate Tm
    my $R = 1.987; # Gas constant in cal/(mol·K)
    my $strand_conc = 500e-9; # Primer concentration 500 nM
    my $salt_conc = 0.1;      # Salt concentration 100 mM (KCl)
    my $mg_conc = 5e-3;       # Mg concentration 5 mM

    # Salt correction for Tm
    my $effective_salt = $salt_conc + 3.3 * sqrt($mg_conc);  # Adjust for Mg2+
    my $salt_correction = 16.6 * log($effective_salt) / log(10);

    # Calculate Tm (convert ΔH from kcal/mol to cal/mol)
    my $tm = (($total_dh * 1000) / ($total_ds + ($R * log($strand_conc)))) - 273.15 + $salt_correction;

    return $tm;
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
    my $reversed = CORE::reverse($string);  # Explicitly call the built-in reverse function
    return $reversed;
}

# Output the annotated dsDNA
sub output_annotated_dsDNA {
    my ($template, $template_annotation, $revcom_annotation) = @_;

    # Reverse the reverse complement and its annotation for proper alignment
    my $reversed_revcom = CORE::reverse($template->{revcom});  # Explicitly call CORE::reverse
    my $reversed_annotation = reverse_with_padding($revcom_annotation);
    
    # Create the base-pairing line
    my $base_pairing = '';
    for my $i (0 .. length($template->{seq}) - 1) {
        $base_pairing .= ':';
    }

    print $out "\n# Annotated dsDNA\n";
    print $out "     $template_annotation\n";
    print $out "   5-" . $template->{seq} . "-3\n";
    print $out "     $base_pairing\n";
    print $out "   3-" . $reversed_revcom . "-5\n";
    print $out "     $reversed_annotation\n";
}


