#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  tangram_filter.pl
#
#        USAGE:  ./tangram_filter.pl  
#
#  DESCRIPTION:  

#         BUGS:  ---
#       AUTHOR:  Jiantao Wu (), 
#  INSTITUTION:  Boston College
#      VERSION:  0.2.1
#      CREATED:  08/29/2012 05:01:04 PM
#===============================================================================

use strict;
use warnings;

use Getopt::Long;
use Scalar::Util qw(looks_like_number);

if (!@ARGV)
{
    ShowHelp();
}

# command line options
my $filterType = "MEI";

# input VCF file
my $vcfFile = "";

# input mask file
my $maskFile = "";

my $outputFile = "";

# threshold for the number of supporting RP fragments
my $rpFragCut = 2;

# threshold for the number of supporting SR fragments
my $srFragCut = 2;

# filter genotype flag
my $filterGenotype = 0;

# show help flag
my $showHelp = 0;

# filter type hash
my %typeHash = (
    "DEL" => 1,
    "DUP" => 1,
    "INV" => 1,
    "MEI" => 1
);

# array of temporary files
my @tmpFiles = ();

# mei family types
my %meiTypes = ();

# mask files
my %maskFiles = ();
my %windowSize = ();

# get the command line arguments
my $result = GetOptions("vcf=s" => \$vcfFile,
                        "msk=s"  => \$maskFile, 
                        "out=s"  => \$outputFile, 
                        "type=s" => \$filterType,
                        "rpf=i"  => \$rpFragCut, 
                        "srf=i"  => \$srFragCut, 
                        "gt"  => \$filterGenotype, 
                        "help" => \$showHelp);

# show the help message if required
if ($showHelp)
{
    ShowHelp();
}

# check the filter type
$filterType = uc($filterType);
if (!defined($typeHash{$filterType}))
{
    DieAndClean("ERROR: $filterType is not a valid filter type.\n");
}

# open the temporary vcf file
my $tmpVcfFileSR = "$vcfFile.tmp.sr.vcf";
my $tmpVcfFileRP = "$vcfFile.tmp.rp.vcf";
open(TMP_SR, ">", $tmpVcfFileSR) or DieAndClean("ERROR: Cannot open the temporary VCF file: $tmpVcfFileSR.\n");
open(TMP_RP, ">", $tmpVcfFileRP) or DieAndClean("ERROR: Cannot open the temporary VCF file: $tmpVcfFileRP.\n");
push(@tmpFiles, $tmpVcfFileSR);
push(@tmpFiles, $tmpVcfFileRP);

# open the vcf file and parse each line of it to remove those low coverage event
my $vcfHeader = "";
if ($vcfFile ne "stdin")
{
    open(IN, "<", $vcfFile) or DieAndClean("ERROR: Cannot open the input VCF file: $vcfFile.\n");
}
else
{
    open(IN, "<-");
}

while (my $line = <IN>)
{
    # save the header
    if ($line =~ /^#/)
    {
        $vcfHeader .= $line;
        next;
    }

    chomp $line;

    my @data = split("\t", $line);

    if ($filterGenotype && $data[9] =~ /^\.:\./)
    {
        next;
    }

    # filter the SV events according to their supporting fragments (reads)
    if ($filterType eq "MEI")
    {
        my $ret = 0;
        if (($ret = SpecialFragFilter($data[7])) > 0)
        {
            my $family = SetSpeicalFamily($data[7]);
            if ($family eq "PO")
            {
                next;
            }

            # events only supported by read pair signal
            # needs futher filtering
            if ($ret == 1)
            {
                print TMP_RP "$line\n";
            }
            else
            {
                print TMP_SR "$line\n";
            }
        }
    }
}

close IN;
close TMP_RP;
close TMP_SR;

# read in the mask file information
open(IN, "<", $maskFile) or DieAndClean("ERROR: Cannot open the input mask list file: $maskFile.\n");
while (my $line = <IN>)
{
    chomp $line;
    next if length($line) == 0;
    next if $line =~ /^#/;

    my @data = split(" ", $line);

    DieAndClean("ERROR: invalid format of the mask list file.\n") if (@data < 3);
    CheckMaskFile(@data);

    $windowSize{uc($data[0])} = $data[1];
    $maskFiles{uc($data[0])} = $data[2];
}
close IN;

# filter the SV events with mask files
if ($filterType eq "MEI")
{
    my $tmpSort = "$vcfFile.tmp.sort.vcf";
    push(@tmpFiles, $tmpSort);

    FilterMEI($tmpSort);
    SortVcf($tmpSort);
}

CleanUp();



#===============================================
# Subroutines
#===============================================
sub SpecialFragFilter
{
    my $infoStr = shift(@_);
    my ($rpFrag5, $rpFrag3, $srFrag5, $srFrag3) = $infoStr =~ /RP5=(\d+);RP3=(\d+);SR5=(\d+);SR3=(\d+)/;

    my $ret = 0;

    if ($srFrag5 >= $srFragCut || $srFrag3 >= $srFragCut)
    {
        $ret = 1;
    }
    elsif ($rpFrag5 >= $rpFragCut && $rpFrag3 >= $rpFragCut)
    {
        $ret = 1;
    }
    elsif (($rpFrag5 >= $rpFragCut || $rpFrag3 >= $rpFragCut) && ($srFrag5 > 0 || $srFrag3 > 0))
    {
        $ret = 1;
    }

    if ($ret == 1 && ($srFrag5 > 0 || $srFrag3 > 0))
    {
        $ret = 2;
    }

    return $ret;
}

sub SetSpeicalFamily
{
    my $alt = shift(@_);
    my $family = 1;
    if ($alt =~ /TYPE=(\w+);/)
    {
        $family = $1
    }

    if ($family ne "PO")
    {
        $meiTypes{$family} = 1;
    }

    return $family;
}

sub CheckMaskFile
{
    my (@data) = @_;
    if (!looks_like_number($data[1])) 
    {
        DieAndClean("ERROR: \"$data[1]\" is not a valid window size in the mask list file.\n");
    }

    DieAndClean("ERROR: \"$data[2]\" is not a valid mask file name.\n") unless (-f $data[2]);
}

sub FilterMEI
{
    my ($tmp) = @_;
    foreach my $type (keys %meiTypes)
    {
        if (defined($maskFiles{$type}))
        {
            my $command = "grep INS:ME:$type $tmpVcfFileRP | windowBed -w $windowSize{$type} -a stdin -b $maskFiles{$type} -v";

            # applied the general masking file
            foreach my $key (keys %maskFiles)
            {
                if (!defined($meiTypes{$key}))
                {
                    $command .= " | windowBed -w $windowSize{$key} -a stdin -b $maskFiles{$key} -v";
                }
            }

            system("$command >> $tmp");
        }
        else
        {
            print STDERR "WARNING: Filter of $type is not provided.\n";
            system("grep INS:ME:$type $tmpVcfFileRP >> $tmpVcfFileSR");
        }
    }

    my $first = 0;
    my $command = "";
    foreach my $key (keys %maskFiles)
    {
        if (!defined($meiTypes{$key}))
        {
            if ($first == 0)
            {
                $command = "windowBed -w $windowSize{$key} -a $tmpVcfFileSR -b $maskFiles{$key} -v";
                $first = 1;
            }
            else
            {
                $command .= " | windowBed -w $windowSize{$key} -a stdin -b $maskFiles{$key} -v";
            }
        }
    }

    if ($command ne "")
    {
        system("$command >> $tmp");
    }
    else
    {
        system("cat $tmpVcfFileSR >> $tmp");
    }

}

sub SortVcf
{
    my ($tmp) = @_;
    my $sortCmd = "perl -lane '\$F[0] =~ s/^chr//g; \$F[0] =~ s/^X/23/; \$F[0] =~ s/^Y/24/; \$F[6] = \"PASS\"; print join(\"\t\", \@F)' $tmp";
    $sortCmd .= " | sort -n -k 1 -n -k 2";
    $sortCmd .= " | perl -ne 's/^23/X/; s/^24/Y/; print \"chr\$_\"'";
    if ($outputFile ne "")
    {
        open(OUT, ">", $outputFile) or DieAndClean("ERROR: Cannot open the output VCF file: $outputFile.\n");
        print OUT "$vcfHeader";
        close OUT;

        $sortCmd .= " >> $outputFile";
    }
    else
    {
        print "$vcfHeader";
    }

    system("$sortCmd") == 0 or DieAndClean ("ERROR: Cannot sort the vcf file.\n");
}

sub CleanUp
{
    foreach my $tmpFile (@tmpFiles)
    {
        system("rm -rf $tmpFile");
    }
}

sub DieAndClean
{
    my ($msg) = @_;
    CleanUp();
    die ($msg);
}

sub ShowHelp
{
    print <<HELP;

    Usage: ./tangram_filter [options] --vcf <input_vcf> --msk <mask_input_list>

    Mandatory Arguments:          --vcf  FILE   input vcf file for filtering
                                  --msk  FILE   input list of mask files with window size information

    Options:                      --type STRING SV event type for filtering, choosing from "DEL", "DUP", "INV" and "MEI" (case insensitive) [MEI]
                                  --rpf  INT    Minimum number of supporting fragments (reads) for read-pair events. For MEI events, this threshold
                                                must be satisfied for reads from both 5' and 3' [2]
                                  --srf  INT    Minimum number of supporting fragments (reads) for split-read events [2]
                                  --out  FILE   Output of filtered and sorted VCF file [stdout]
                                  --help        Show this help message

    Note:
        1. This script require the installation of "bedtools" package and Unix
           sort in the default directory.

        2. Each entry of the list of mask files is a tab delimited file 
           with following format:
        
           "TYPE WINDOW_SIZE FILE_NAME"

           "TYPE" (string) is the type of this mask file. For a referenced MEI 
           mask file, it must match the first two characters of the family name 
           in the VCF file (For example AL: ALU, L1: L1, SV: SVA and HE: HERV).
           This mask file will only be applied to the corresponding type of MEI
           events. For example, AL mask file will only be applied to ALU insertions. 
           The rest of the mask files, such as segmental duplication
           mask and simple repeat mask, their "TYPE" string can be anything and
           it will be applied to all the entries in the VCF file. No space is allowed
           in the type name.

           "WINDOW_SIZE" (integer) is the window size around each entry of
           the mask file.

           "FILE_NAME" (string) is the path to the corresponding mask file

           All the mask files must be in the BED format. For detailed information
           about this format, please check http://genome.ucsc.edu/FAQ/FAQformat.html

HELP

    exit(0);
}
