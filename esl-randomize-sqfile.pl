#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::Easel::SqFile;
use Bio::Easel::Random;

my $do_outfile         = 0;
my $ofh                = "";
my $do_length_fraction = 0;
my $do_split           = 0;
my $nsplit;
my $cur_split_idx      = 1;
my $do_info            = 0;
my $outfile_root       = "";
my $outfile;
my $seed = 1801;
my $do_sfetch_output   = 0;
my $required_tail      = "";

&GetOptions( "O=s" => \$outfile_root, 
             "I"   => \$do_info, 
             "L"   => \$do_length_fraction,
             "N=s" => \$nsplit, 
             "F"   => \$do_sfetch_output,
             "S=s" => \$seed);


# Option processing
if(defined $nsplit) { 
  $do_split = 1; 
  $cur_split_idx = 1;
}
else { 
  $do_split = 0;
  $nsplit = 1;
}
if($do_sfetch_output) { 
  $required_tail = ".sfetch";
}
else { 
  $required_tail = ".fa";
}  
if($outfile_root ne "") { 
  if($outfile_root !~ m/$required_tail$/) { 
    if($do_sfetch_output) { die "ERROR, with -O and -F, output file must end with .sfetch"; }
    else                  { die "ERROR, with -O, output file must end with .fa"; }
  }
}
  
if($outfile_root ne "") { 
  $do_outfile = 1;
  $outfile = $outfile_root;
  $cur_split_idx = 1;
  if($do_split) { 
    $outfile_root =~ s/$required_tail$//;
    $outfile = $outfile_root . "_" . $cur_split_idx . $required_tail; 
  }
  open(OUT, ">" . $outfile) || die "ERROR unable to open $outfile for writing"; 
}

my $usage;
$usage  = "esl-randomize-sqfile.pl [OPTIONS] <seqfile to randomize> <F: fraction of sequences to include>\n";
$usage .= "esl-randomize-sqfile.pl -I        <seqfile to get info on>\n\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t-I    : print info on seq file (nseq, total residues, etc) and exit\n";
$usage .= "\t\t-O <f>: output to file <f>, not stdout\n";
$usage .= "\t\t-L    : include F fraction of residues, (and split based on # residues) not seqs\n";
$usage .= "\t\t-N <n>: split output into <n> files (requires -O)\n";
$usage .= "\t\t-S <n>: seed RNG with seed <n> (default: 1801)\n\n";
$usage .= "\tEXAMPLES:\n";
$usage .= "\t\t'esl-randomize-sqfile.pl input.fa 0.5':\n";
$usage .= "\t\t\toutput half the sequences in input.fa in random order to stdout.\n\n";
$usage .= "\t\t'esl-randomize-sqfile.pl -O new.fa -N 10 input.fa 0.10':\n";
$usage .= "\t\t\toutput 10\% of the sequences in input.fa split into 10 files [new.fa.1 --> new.fa.10]\n";
$usage .= "\t\t\tin random order\n\n";
$usage .= "\t\t'esl-randomize-sqfile.pl -I new.fa\n";
$usage .= "\t\t\toutput information on seqs in new.fa and exit (don't randomize new.fa)\n\n\n";
$usage .= "\tNOTE: it's currently not possible to output seqs in the order\n";
$usage .= "\t      they appear in the input file (due to SSI dependence)\n";
$usage .= "\t      they will always be in random order\n\n";

if((! $do_info && scalar(@ARGV) != 2) || ($do_info && scalar(@ARGV != 1))) { die $usage; }
my ($in_sqfile, $fraction);
if($do_info) { 
  $in_sqfile = $ARGV[0]; 
  $fraction = 1.0;
}
else {
  ($in_sqfile, $fraction) = @ARGV;
}

# option/command-line checks
if($fraction > 1.0 || $fraction <= 0.) { print "ERROR fraction F must fall in range (0..1]\n"; die $usage; }
if($do_split && (! $do_outfile))  { print "ERROR -N requires -O\n"; die $usage; }
if($do_info && ($do_outfile || $do_split || $do_length_fraction))  { print "ERROR -I incompatible with all other options\n"; die $usage; }

# variables
my $seqname; # sequence name
my $L;       # sequence length
my $i;       # counter

# open sequence file
my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $in_sqfile });

# create RNG
my $rng = Bio::Easel::Random->new({ seed => $seed });
# get number of seqs
my $nseq = $sqfile->nseq_ssi();

# if -L enabled, get total length of all seqs
my $totL;
my $maxL = 0;
my $minL = -1;
if($do_length_fraction || $do_info) { 
  for($i = 0; $i < $nseq; $i++) { 
    (undef, $L) = $sqfile->fetch_seq_name_and_length_given_ssi_number($i);
    $totL += $L;
    if($L > $maxL) { $maxL = $L; }
    if($minL == -1 || $L < $minL) { $minL = $L; }
  }
}
if($do_info) { 
  printf("Number of sequences: $nseq\n");
  printf("Total length:        %.6f Mb\n", $totL / 1000000.); 
  printf("Average length:      %.6f Mb\n", $totL / ($nseq * 1000000.)); 
  printf("Maximum length:      %.6f Mb\n", $maxL / 1000000.);
  printf("Minimum length:      %.6f Mb\n", $minL / 1000000.);
  exit(0);
}

# determine stopping condition
my $ndesired_tot; # number of total sequences we'll pick
my $ndesired_cur; # number of sequences we'll pick to current file
my $Ldesired_tot; # total length of sequences we'll pick
my $Ldesired_cur; # length of sequences we'll pick to current file
if($do_length_fraction) { 
  $Ldesired_tot = int($totL * $fraction);
  $Ldesired_cur = int($Ldesired_tot / $nsplit);
  if($Ldesired_tot % $nsplit != 0) { $Ldesired_cur++; }
  #printf("Ldesired_cur: $Ldesired_cur\n");
}
else { 
  $ndesired_tot = int($nseq * $fraction);
  $ndesired_cur = int($ndesired_tot / $nsplit);
  if($ndesired_tot % $nsplit != 0) { $ndesired_cur++; }
}

# We define an array @mapA with an element for each seq.
# Initially $mapA[$i] == $i, but when if we pick seq $i
# we set $mapA[$i] to $mapA[$nremaining-1], then choose
# a random int between 0 and $nremaining-1. This gets us
# a random sample without replacement.
my @mapA = ();
for($i = 0; $i < $nseq; $i++) { $mapA[$i] = $i; }

# randomly pick seqs, and fetch them to stdout as we go
my $nremaining = $nseq;
my $ndesired   = int($nseq * $fraction);
my $npicked_tot = 0; # number seqs output thus far
my $npicked_cur = 0; # number seqs output to current file (will be equal to $npicked_tot unless -S)

my $Lpicked_tot = 0;  # total length of all sequences output thus far
my $Lpicked_cur = 0;  # total length of sequences output to current file 
my $idx         = 0; # index of a sequence in SSI
my $keep_going  = 1; # TRUE to keep picking, set to FALSE differently if -L used

while($keep_going) { 
  $idx = $rng->roll($nremaining - 1);
  if($mapA[$idx] == -1) { die "ERROR algorithm for randomly picking seqs is flawed..." }

  ($seqname, $L) = $sqfile->fetch_seq_name_and_length_given_ssi_number($mapA[$idx]);

  # deal with -S, potentially close current file and open a new one
  if($do_split) { 
    if(((! $do_length_fraction) && $npicked_cur == $ndesired_cur) || # simple case, sequence number based splitting
       ($do_length_fraction &&                                       # more complicated case, residue number based splitting
        (($Lpicked_cur >= $Ldesired_cur) ||                            # we've exceeded our limit OR
         (($Lpicked_cur > (0.5 * $Ldesired_cur)) &&                    # we've exceeded half our limit AND
          $L > (($Ldesired_cur - $Lpicked_cur) * 2))))) {              # next sequence to add is more than twice length adding next sequence of length L would exceed limit by more than L/2
      close OUT;
      printf("Finished writing %-20s [%7d seqs; %10d residues]\n", $outfile, $npicked_cur, $Lpicked_cur);
      if($cur_split_idx < $nsplit) { 
        $cur_split_idx++;
        $outfile = $outfile_root . "_" . $cur_split_idx . $required_tail;
        open(OUT, ">" . $outfile) || die "ERROR unable to open $outfile for writing"; 
      }
      $npicked_cur = 0;
      $Lpicked_cur = 0;
    }
  }

  my $seqstring;
  if($do_sfetch_output) { # print only sequence name
    if(defined $outfile) { print OUT $seqname . "\n"; }
    else                 { print     $seqname . "\n"; }
  }
  else { # print actual sequence
    my $seqstring = $sqfile->fetch_seq_to_fasta_string($seqname, 60);
    if(defined $outfile) { print OUT $seqstring; }
    else                 { print     $seqstring; }
    undef $seqstring;
  }
  
  $Lpicked_tot += $L;
  $Lpicked_cur += $L;
  $npicked_tot++;
  $npicked_cur++;
  
  #printf("Lpicked_tot: $Lpicked_tot\n");
  #printf("Lpicked_cur: $Lpicked_cur\n");
  #printf("Ldesired_cur: $Ldesired_cur\n");

  # update mapA
  if($idx != ($nremaining-1)) { # edge case
    $mapA[$idx] = $mapA[($nremaining-1)];
  }
  $mapA[($nremaining-1)] = -1; # sanity check; if we pick one with a -1 value, we'll know something's wrong
  $nremaining--;

  if(($do_length_fraction     && $Lpicked_tot >= $Ldesired_tot) || 
     ((! $do_length_fraction) && $npicked_tot == $ndesired_tot)) { 
    $keep_going = 0;
  }
}  
if($outfile ne "" && tell(OUT) != -1) { 
  close(OUT); 
  if(-s $outfile) { 
    printf("Finished writing %-20s [%7d seqs; %10d residues]\n", $outfile, $npicked_cur, $Lpicked_cur);
  }
  else { # final file was empty, unlink it
    unlink $outfile;
  }
}

