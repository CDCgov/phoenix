#!/usr/bin/perl

# create_sample_sheet.pl v0.1.0
# Authors: @hseabolt, @mjcipriano, @sateeshperi

# SYNOPSIS:
# Reads the contents of a given directory (expected to contain sequencing read data files in FASTQ format, optionally gzip-compressed),
# then uses standard file naming conventions to match up samples with their read mate information (including potentially multiple sequencing lanes),
# and finally, formats and outputs a csv sample sheet.

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# GetOpts Variable declarations
my $input = "--";
my $output = "--";
my $prefix;
my $fullpath;
my $sanitize;
my $version;
my $help;

sub usage {
    my $usage = "\ncreate_sample_sheet.pl\nv0.1.0\n
    PURPOSE: Reads the contents of a given directory (expected to contain sequencing read data files in FASTQ format, optionally gzip-compressed),
        then uses standard file naming conventions to match up samples with their read mate information (including potentially multiple sequencing lanes),
        and finally, formats and outputs a csv sample sheet.

    USAGE:    create_sample_sheet.pl -i <input directory> -o <output>

    ARGUMENTS:
    -i | --input        DIR (Required).  Input directory name.
    -o | --output        STR (Optional).  Output sample sheet file name.  If no argument given, prints to STDOUT.
    -p | --prefix        STR (Optional).  Prepend a STR prefix to sequence names (to set directories not on current system).
    -f | --fullpath     (Optional flag). Call realpath and set the fully qualified path in the sequence ID fields.
    -s | --sanitize        (Optional flag). Find and replace spaces in the sequence ID field.
    -v | --version         Print version number and exit.
    -h | --help        Print this help message and exit.

    \n";
    print $usage;
}

GetOptions(    'input|i=s'  => \$input,
            'out|o=s'    => \$output,
            'prefix|p=s' => \$prefix,
            'sanitize|s' => \$sanitize,
            'fullpath|f' => \$fullpath,
            'version|v'  => \$version,
            'help|h'     => \$help,
) or die usage();

# Print the version number or the help message and exit if -v or -h is activated
if ( $version )     { die "create_sample_sheet.pl version 0.1.0\n"; }
if ( $help    )     { die usage();                                            }

################################################################################
# Read the input directory, capturing all files that end in .fastq (optionally .gz) or some flavor of this extension
opendir( my $dh, $input ) or die "create_sample_sheet::ERROR --> Cannot read the given input directory.  $!\n";
my @fastqs;
while ( readdir $dh )    {
    chomp $_;

    # Keep this file if it is a FASTQ file
    # WARNING: A curse upon ye who put .fastq somewhere in the filename other than the extension
    # This *should* work even with full paths that might contain 'fastq' in the path, e.g. path/to/data/fastq/file1
    push @fastqs, $_ if ( $_ =~ /\.fq/ || $_ =~ /\.fastq/ );
}

# For each FASTQ file:
# 1) deconstruct the string to get the sample/file name, lane info, and read information.
# 2) Then, populate a lookup table to store matching information per sample.
my %Samples = ();
foreach my $file ( @fastqs )    {

    # Deconstruct the sample name -- here I'm doing it with regex capture groups
    my $filetmp = $file;
    $filetmp =~ s/\s/_/g;
    my ( $id, $fr_read ) = $filetmp =~ /(.*)_[Rr]?(1|2)/;
    my ( $lane ) = $id =~ s/_L(\d{3})?//;
    $lane = ( $1 eq $id )? "001" : $1;            # Default to 001 if there is no lane information

    # Add the file to the lookup table
    # Hierarchy is SAMPLENAME ==> LANE ==> READ (1|2) = FILENAME
    $Samples{$id}->{$lane}->{$fr_read} = $file;
}

# Set output filehandles
my $fileout = open( OUT, ">", "$output" ) if $output ne "--";
my $fhout;
if ( $fileout )        {    $fhout = *OUT;            }
else                {    $fhout = *STDOUT;        }

# Print out a csv sample sheet with format <name>, <r1>, <r2>, <r1b>, <r2b>, ... <r1n>, <r2n>
foreach my $sample ( keys %Samples )        {

    # Prepend the prefix if one is set by the user.
    # Also do a quick sanity check to remove any double // that might be introduced in various ways.
    my @line_out = ( $sample );
    $line_out[0] =~ s/\/\//\//;            # Well that's an ugly subst regex...
    $line_out[0] =~ s/\s/_/g if ( $sanitize );

    # Loop to get the lanes/read 1 or 2 information
    foreach my $lane ( sort {$a <=> $b} keys %{$Samples{$sample}} )    {
        foreach my $read ( sort {$a <=> $b} keys %{$Samples{$sample}->{$lane}} )    {
            if ( $fullpath )     {
                my $filetmp = $Samples{$sample}->{$lane}->{$read};
                my $realpath = `realpath "$input/$filetmp"`;
                chomp $realpath;
                push @line_out, $realpath;
            }
            else    {
                if ( $prefix )     {
                    push @line_out, "$prefix/$Samples{$sample}->{$lane}->{$read}";        # May need to add $input between prefix/INPUT/name...
                }
                else     {
                    push @line_out, "$input/$Samples{$sample}->{$lane}->{$read}";
                }
            }
        }
    }
    print $fhout join(",", @line_out), "\n";
}
close $fhout if ( $fileout );

exit;
