#!/usr/bin/perl -w
## calculate snp number in specific genome window for each individual

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw{uniq};


my ($in, $out, $popmap, $window, $help);
GetOptions (
    "in=s"      =>\$in,
    "out=s"     =>\$out,
    "popmap=s"  =>\$popmap,
    "window=i"  =>\$window,
    "help:1"    =>\$help
) or die "$!";

my $usage = '
    -i: input vcf file.
    -o: output file.
    -p: popmap file.
    -w: window size.
    -h: help message.
';

die $usage if ($help || !($in && $out && $popmap && $window));

`vcftools --vcf $in --SNPdensity $window`;
` sed '1d' out.snpden |awk '{print \$1"\t"\$2"\t"\$2+'$window'}' > vcf_snpden`;
`rm out.snpden`;
`vcftools --vcf $in --extract-FORMAT-info GT`;

my($in_fh);

my @ind = parse_popmap($popmap);
sub parse_popmap {
    my @ind;
    my $pop = shift;
    open(my $in_fh, $pop) or die "$!";
    while(<$in_fh>) {
        next if /^#|^$/;
        $_ =~ s/\r//g;
        chomp;
        # indv sex
        my @parts = split;
        die if scalar(@parts) ne 2;
        push @ind, $parts[0];
    }
    close $in_fh;
    return @ind;
    
}


my @window_list;
open (my $bed, "vcf_snpden") or die "$!";
while(<$bed>){
    my (@parts,$chr,$start,$end);
    next if /^#|^$/;
    chomp;
    @parts=split;
    ($chr,$start,$end) = ($parts[0], $parts[1], $parts[2]);
    my $region = "$chr:$start-$end";
    push @window_list, $region;
}
close $bed;

open($in_fh, "out.GT.FORMAT") or die "$!";
open(my $out_fh, ">$out") or die "$!";
my $win_i=0;
my @snp_num;
foreach my $ind_i(0..$#ind){
	$snp_num[$ind_i]=0;
}
while (<$in_fh>) {
    my ($chr, $start, $end);
    my ($chrom, $pos, @gt, @parts);
    next if /^#|^$/;
    chomp;
    @parts = split;
    $chrom = $parts[0];
    $pos   = $parts[1];
    if ($chrom eq 'CHROM') {
        print $out_fh join ("\t", "chrom","win_start","win_end",@ind),"\n";
        next;
    }
    my @region_part = split /:|-/, $window_list[$win_i];
    $chr   = $region_part[0];
    $start = $region_part[1];
    $end   = $region_part[2];
    if ($chrom eq $chr && $pos >= $start && $pos < $end) {
		
        foreach my $ind_i (2..$#parts){
            $snp_num[$ind_i-2]++ if ($parts[$ind_i] eq '0/1' or $parts[$ind_i] eq '0|1');
        }                
    } else {
        print $out_fh join ("\t", $chr, $start, $end, @snp_num), "\n";
        $win_i++;
        foreach my $ind_i (2..$#parts){
            if ($parts[$ind_i] eq '0/1' or $parts[$ind_i] eq '0|1'){
                $snp_num[$ind_i-2]=1;
            } else {
                $snp_num[$ind_i-2]=0;                    
            }
            
        }                         
    }
}
close $in_fh;
close $out_fh;

`rm out.GT.FORMAT && rm vcf_snpden`;
