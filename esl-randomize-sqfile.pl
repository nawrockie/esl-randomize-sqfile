#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::Easel::SqFile;
use Bio::Easel::Random;
use Time::HiRes qw(gettimeofday);

# make *STDOUT file handle 'hot' so it automatically flushes whenever we print to it
# it is printed to
select *STDOUT;
$| = 1;

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
my $do_gzip            = 0;

&GetOptions( "O=s" => \$outfile_root, 
             "I"   => \$do_info, 
             "L"   => \$do_length_fraction,
             "N=s" => \$nsplit, 
             "F"   => \$do_sfetch_output,
             "Z"   => \$do_gzip,
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
$usage  = "esl-randomize-sqfile.pl v0.02\n\n";
$usage  = "Usage:\n";
$usage  = "esl-randomize-sqfile.pl [OPTIONS] <seqfile to randomize> <F: fraction of sequences to include>\n";
$usage .= "esl-randomize-sqfile.pl -I        <seqfile to get info on>\n\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t-I    : print info on seq file (nseq, total residues, etc) and exit\n";
$usage .= "\t\t-O <f>: output to file <f>, not stdout\n";
$usage .= "\t\t-L    : include F fraction of residues, (and split based on # residues) not seqs\n";
$usage .= "\t\t-N <n>: split output into <n> files (requires -O)\n";
$usage .= "\t\t-S <n>: seed RNG with seed <n> (default: 1801)\n";
$usage .= "\t\t-Z    : gzip each file after it is completed\n\n";
$usage .= "\tEXAMPLES:\n";
$usage .= "\t\t'esl-randomize-sqfile.pl input.fa 0.5':\n";
$usage .= "\t\t\toutput half the sequences in input.fa in random order to stdout.\n\n";
$usage .= "\t\t'esl-randomize-sqfile.pl -O new.fa -N 10 input.fa 0.10':\n";
$usage .= "\t\t\toutput 10\% of the sequences in input.fa split into 10 files [new.fa.1 --> new.fa.10]\n";
$usage .= "\t\t\tin random order\n\n";
$usage .= "\t\t'esl-randomize-sqfile.pl -I new.fa\n";
$usage .= "\t\t\toutput information on seqs in new.fa and exit (don't randomize new.fa)\n\n\n";
$usage .= "\tNOTE: See esl-ssplit.pl in Bio-Easel for more efficient randomization of\n";
$usage .= "\t      sequence files, that doesn't use SSI. *This* script is only superior\n";
$usage .= "\t      to esl-ssplit.pl for one use case: when you want a random subset of the\n";
$usage .= "\t      sequences.\n\n";

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
if($do_info && ($do_outfile || $do_split || $do_length_fraction || $do_gzip))  { print "ERROR -I incompatible with all other options\n"; die $usage; }
if($do_gzip && $do_sfetch_output) { print "ERROR -F and -Z are incompatible.\n"; die $usage; }

# variables
my $seqname; # sequence name
my $L;       # sequence length
my $i;       # counter

# open sequence file
my $progress_w = 50;
my $start_secs = output_progress_prior("Opening the sequence file", $progress_w, *STDOUT);
my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $in_sqfile });
output_progress_complete($start_secs, undef, *STDOUT);

# create RNG
my $rng = Bio::Easel::Random->new({ seed => $seed });

# get number of seqs
$start_secs = output_progress_prior("Getting number of seqs and possibly indexing", $progress_w, *STDOUT);
my $nseq = $sqfile->nseq_ssi();
output_progress_complete($start_secs, undef, *STDOUT);

# if -L enabled, get total length of all seqs
my $totL = 0;
if($do_length_fraction || $do_info) { 
  $start_secs = output_progress_prior("Getting total number of nucleotides", $progress_w, *STDOUT);
  $totL = $sqfile->nres_ssi();
  output_progress_complete($start_secs, undef, *STDOUT);
}

if($do_info) { 
  printf("Number of sequences: $nseq\n");
  printf("Total length:        %.6f Mb\n", $totL / 1000000.); 
  printf("Average length:      %.6f Mb\n", $totL / ($nseq * 1000000.)); 
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

my $Lpicked_tot = 0; # total length of all sequences output thus far
my $Lpicked_cur = 0; # total length of sequences output to current file 
my $idx         = 0; # random index in @mapA
my $seqidx      = 0; # index of a sequence in SSI (mapA[$idx])
my $keep_going  = 1; # TRUE to keep picking, set to FALSE differently if -L used

$start_secs = output_progress_prior("Writing file $outfile ", $progress_w, *STDOUT);
while($keep_going) { 
  $idx = $rng->roll($nremaining - 1);
  $seqidx = $mapA[$idx];
  if($seqidx == -1) { die "ERROR algorithm for randomly picking seqs is flawed..." }

  ($seqname, $L) = $sqfile->fetch_seq_name_and_length_given_ssi_number($seqidx);

  # deal with -S, potentially close current file and open a new one
  if($do_split) { 
    if(((! $do_length_fraction) && $npicked_cur == $ndesired_cur) || # simple case, sequence number based splitting
       ($do_length_fraction &&                                       # more complicated case, residue number based splitting
        (($Lpicked_cur >= $Ldesired_cur) ||                            # we've exceeded our limit OR
         (($Lpicked_cur > (0.5 * $Ldesired_cur)) &&                    # we've exceeded half our limit AND
          $L > (($Ldesired_cur - $Lpicked_cur) * 2))))) {              # next sequence to add is more than twice length adding next sequence of length L would exceed limit by more than L/2
      close OUT;
      output_progress_complete($start_secs, sprintf("(%d seqs; %d residues)", $npicked_cur, $Lpicked_cur), *STDOUT);
      if($do_gzip) { 
        $start_secs = output_progress_prior("Gzip'ing $outfile ", $progress_w, *STDOUT);
        run_command("gzip -f $outfile");
        output_progress_complete($start_secs, undef, *STDOUT);
      }
      if($cur_split_idx < $nsplit) { 
        $cur_split_idx++;
        $outfile = $outfile_root . "_" . $cur_split_idx . $required_tail;
        open(OUT, ">" . $outfile) || die "ERROR unable to open $outfile for writing"; 
        $start_secs = output_progress_prior("Writing file $outfile ", $progress_w, *STDOUT);
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
    my $seqstring = $sqfile->fetch_seq_to_fasta_string_given_ssi_number($seqidx, 60);
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
    output_progress_complete($start_secs, sprintf("[%7d seqs; %10d residues]", $npicked_cur, $Lpicked_cur), *STDOUT);
    if($do_gzip) { 
      $start_secs = output_progress_prior("Gzip'ing $outfile ", $progress_w, *STDOUT);
      run_command("gzip -f $outfile");
      output_progress_complete($start_secs, undef, *STDOUT);
    }
  }
  else { # final file was empty, unlink it
    unlink $outfile;
  }
}

#################################################################
# Subroutine:  run_command()
# Incept:      EPN, Mon Dec 19 10:43:45 2016 [ribo_RunCommand()]
#
# Purpose:     Runs a command using system() and exits in error 
#              if the command fails. 
#
# Arguments:
#   $cmd:         command to run, with a "system" command;
#
# Returns:    void
#
# Dies:       if $cmd fails
#################################################################
sub run_command {
  my $sub_name = "run_command()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd) = @_;
  
  system($cmd);

  if($? != 0) { 
    die "ERROR in $sub_name, the following command failed:\n$cmd\n";
  }
}

#################################################################
# Subroutine : output_progress_prior()
# Incept:      EPN, Fri Feb 12 17:22:24 2016 [outputProgressPrior() dnaorg.pm]
#              
# Purpose:      Output to $FH1 (and possibly $FH2) a message indicating
#               that we're about to do 'something' as explained in
#               $outstr.  
#
#               Caller should call *this* function, then do
#               the 'something', then call outputProgressComplete().
#
#               We return the number of seconds since the epoch, which
#               should be passed into the downstream
#               outputProgressComplete() call if caller wants to
#               output running time.
#
# Arguments: 
#   $outstr:     string to print to $FH
#   $progress_w: width of progress messages
#   $FH:         file handle to print to
# 
# Returns:     Number of seconds and microseconds since the epoch.
#
################################################################# 
sub output_progress_prior { 
  my $nargs_expected = 3;
  my $sub_name = "output_progress_prior()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($outstr, $progress_w, $FH) = @_;

  if(defined $FH) { printf $FH ("# %-*s ... ", $progress_w, $outstr); }

  return secondsSinceEpoch();
}

#################################################################
# Subroutine : output_progress_complete()
# Incept:      EPN, Fri Feb 12 17:28:19 2016 [outputProgressComplete() dnaorg.pm]
#
# Purpose:     Output to $FH a 
#              message indicating that we've completed 
#              'something'.
#
#              Caller should call *this* function,
#              after both a call to output_progress_prior()
#              and doing the 'something'.
#
#              If $start_secs is defined, we determine the number
#              of seconds the step took, output it, and 
#              return it.
#
# Arguments: 
#   $start_secs:    number of seconds either the step took
#                   (if $secs_is_total) or since the epoch
#                   (if !$secs_is_total)
#   $extra_desc:    extra description text to put after timing
#   $FH:            file handle to print to
# 
# Returns:     Number of seconds the step took (if $secs is defined,
#              else 0)
#
################################################################# 
sub output_progress_complete { 
  my $nargs_expected = 3;
  my $sub_name = "output_progress_complete()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($start_secs, $extra_desc, $FH) = @_;

  my $total_secs = undef;
  if(defined $start_secs) { 
    $total_secs = secondsSinceEpoch() - $start_secs;
  }

  if(defined $FH) { printf $FH ("done."); }

  if(defined $total_secs || defined $extra_desc) { 
    if(defined $FH) { printf $FH (" ["); }
  }
  if(defined $total_secs) { 
    if(defined $FH) { printf $FH (sprintf("%.1f seconds%s", $total_secs, (defined $extra_desc) ? ", " : "")); }
  }
  if(defined $extra_desc) { 
    if(defined $FH) { printf $FH $extra_desc };
  }
  if(defined $total_secs || defined $extra_desc) { 
    if(defined $FH) { printf $FH ("]"); }
  }

  if(defined $FH) { printf $FH ("\n"); }
  
  return (defined $total_secs) ? $total_secs : 0.;
}

#################################################################
# Subroutine : secondsSinceEpoch()
# Incept:      EPN, Sat Feb 13 06:17:03 2016
#
# Purpose:     Return the seconds and microseconds since the 
#              Unix epoch (Jan 1, 1970) using 
#              Time::HiRes::gettimeofday().
#
# Arguments:   NONE
# 
# Returns:     Number of seconds and microseconds
#              since the epoch.
#
################################################################# 
sub secondsSinceEpoch { 
  my $nargs_expected = 0;
  my $sub_name = "secondsSinceEpoch()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($seconds, $microseconds) = gettimeofday();
  return ($seconds + ($microseconds / 1000000.));
}
