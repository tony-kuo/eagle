#! /usr/bin/perl -w
use strict;
use warnings;
use Fatal;
use IO::File;
use IO::Handle;
use Getopt::Std;
sub usage {
print  <<EOF;
Usage: $0 A_B.txt A_C.txt A_D.txt B_C.txt B_D.txt C_D.txt
  -b  make bedfile by giving the fasta.fai index
EOF
exit;
}
!@ARGV && &usage();

my %opt = (); &getopts( "b:", \%opt );
my $faidx = defined $opt{b} ? $opt{b} : ""; 

my %ab = &read_file($ARGV[0]);
my %ac = &read_file($ARGV[1]);
my %ad = &read_file($ARGV[2]);
my %bc = &read_file($ARGV[3]);
my %bd = &read_file($ARGV[4]);
my %cd = &read_file($ARGV[5]);

my %bed = ();
if (length($faidx) > 0) {
    my $fh = &open_maybe_compressed($faidx);
    while (my $line = $fh->getline) {
        chomp $line;
        my @t = split(/\s+/, $line);
	$bed{$t[0]} = "$t[0]\t0\t" . ($t[1]-1);
    }
}

while (my ($a, $b) = each %ab) {
    next if (!exists $ac{$a} || !exists $bc{$b}); # A in AC && B in BC
    next if (!exists $ad{$a} || !exists $bd{$b}); # A in AD && B in BD
    next if ($ac{$a} ne $bc{$b} || $ad{$a} ne $bd{$b}); 

    my $c = $ac{$a};
    my $d = $ad{$a};
    next if (!exists $cd{$c});
    next if ($cd{$c} ne $d);

    if (length($faidx) > 0) {
       print "$bed{$a}\n$bed{$b}\n$bed{$ad{$a}}\n";
    }
    else {
        print "$a\t$b\t$c\t$d\n";
    }
}

# --------------------------------------------------------------------------------------------------
sub open_maybe_compressed {
    my($fname) = @_;
    my @try = ("pbzip2 -dc", "pigz -dc", "bzip2 -dc", "gzip -dc", "cat");
    my $fh;
    for my $exe (@try) {
        my $io = IO::File->new("$exe \Q$fname\E 2> /dev/null |");
        my $c = $io->getc;
        if (defined $c) {
            $io->ungetc(ord $c);
            #print STDERR "Using $exe for $fname\n";
            return $io;
        }
    }
    die "could not open $fname";
}

sub read_file {
    my($fname) = @_;
    my %data = ();
    my $fh = &open_maybe_compressed($fname);
    while (my $line = $fh->getline) {
        chomp $line;
        my @t = split(/\s+/, $line);
        if (scalar @t == 3) {
            $data{$t[1]} = $t[2];
        }
        else {
            $data{$t[0]} = $t[1];
        }
    }
    return %data;
}
