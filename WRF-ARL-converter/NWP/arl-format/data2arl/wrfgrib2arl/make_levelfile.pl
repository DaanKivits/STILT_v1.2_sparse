#!/usr/bin/perl

#  # $Id: make_levelfile.pl,v 1.3 2005/08/25 17:53:45 trn Exp $

use strict;

use Getopt::Long ;
use File::Basename;

my ($verbose, $help, $gribfile, $halffile, $fullfile, $etacode, $wcode);
$verbose=1;
my %options = (
                "verbose!"        => \$verbose, 
                "help|?!"         => \$help,
                "gribfile=s"       => \$gribfile,
                "halffile=s"       => \$halffile,
                "fullfile=s"       => \$fullfile,
                "etacode=i"       => \$etacode,
                "wcode=i"       => \$wcode,
              );
my %helpopt = (
	       "verbose!"        => "printout control", 
	       "help|?!"         => "generate help message",
                "gribfile=s"       => "Required: input grib file",
                "halffile=s"       => "name of half-levels file (output)",
                "fullfile=s"       => "name of full-levels file (output)",
                "etacode=i"       =>"GRIB parameter code (kpds5) for eta",
                "wcode=i"       =>"GRIB parameter code (kpds5) for w
  (specify same as etacode if full and half levels are in GRIB file for kpds5=etacode)"
	      );
$help=1 unless GetOptions (%options);
$help=1 unless $gribfile; #required
$halffile="half_levels.txt" unless $halffile; #default
$fullfile="full_levels.txt" unless $fullfile; #default
$etacode=149 unless $etacode; #default
$wcode=$etacode unless $wcode; #default
$verbose = 1 if $help;

my @call = caller(0);
if (@call > 6) {
    $0 = $call[6] if ($call[3] eq "(eval)");
}
if ($verbose) {
  print "At ",scalar(localtime)," running $0 with:\n";
  foreach my $key (keys(%options)) {
    if (substr($key,length($key)-1) eq "@") {
      print "  $key = @{$options{$key}}\n";
    } else {
      print "  $key = ${$options{$key}}\n";
    }
  }
}
&help(\%options, \%helpopt) if $help;

my ($cmd0,$pwd0,$cmde)=fileparse($0,'\.pl');
my $tmpfile="$cmd0" . "." . "$$";
my $cmd="wgrib $gribfile | grep kpds5=$etacode | wgrib -i -text $gribfile -o $tmpfile";
print "Running $cmd\n" if $verbose;
if (system("$cmd")) {
  die "Error running: $cmd\n";
}
open(TMP,"<$tmpfile") or die "Could not open tempfile $tmpfile\n";
my @rawlevs;
my $ilin=0;
foreach my $line (<TMP>) {
  $ilin++;
  next unless (2*int($ilin/2) == $ilin); #skip odd-numbered lines
  push @rawlevs,$line;
}
close(TMP) or warn "Could not close tempfile $tmpfile\n";
unlink $tmpfile or warn "Could not unlink tempfile $tmpfile\n";

# sort numerically descending: bottom up
my @sorted = sort {$b <=> $a} @rawlevs;
my (@half_levs, @full_levs);
my $ioff=0;
if ($etacode == $wcode) {
  for (my $i=0;$i < @sorted; $i=$i+2) {
    push @full_levs,$sorted[$i];
  }
  $ioff=1;
} else {
  $cmd="wgrib -o /dev/null $gribfile | grep kpds5=$wcode | sed -e 's/[^:]*:[^:]*:[^:]*:[^:]*:[^:]*:[^:]*:kpds7=//' -e 's/:.*//' > $tmpfile";
  print "Running $cmd\n" if $verbose;
  if (system("$cmd")) {
    die "Error running: $cmd\n";
  }
  open(TMP,"<$tmpfile") or die "Could not open tempfile $tmpfile\n";
  @full_levs = <TMP>;
  chomp, @full_levs;
  close(TMP) or warn "Could not close tempfile $tmpfile\n";
  unlink $tmpfile or warn "Could not unlink tempfile $tmpfile\n";
  @full_levs = sort {$b <=> $a} @full_levs
}
for (my $i=$ioff;$i < @sorted; $i=$i+1+$ioff) {
  push @half_levs,$sorted[$i];
}
foreach my $ltype ('half','full') {
  my $fname=$ltype eq 'half' ? $halffile : $fullfile;
  my @levs=$ltype eq 'half' ? @half_levs : @full_levs;
  open(OUTF,">$fname") or die "Could not open $ltype level output file $fname\n";
  print OUTF scalar(@levs)."\n";
  print OUTF "Level values (from kpds5=etacode=$etacode) for $ltype levels\n" if ($ltype eq 'half');
  print OUTF "Level values (from kpds5=wcode=$wcode) for $ltype levels\n" if ($ltype eq 'full');
  print OUTF "Created by ".'$Id: make_levelfile.pl,v 1.3 2005/08/25 17:53:45 trn Exp $'."\n";
  print OUTF "        from input grib file $gribfile\n";
  print OUTF "ENDHEADER\n";
  for (my $i=0; $i < @levs; $i++) {
    print OUTF $levs[$i];
  }
  close(OUTF) or warn "Could not close output file $fname\n";
}
1;

sub help {
  my %options=%{$_[0]};
  my %helpopt=%{$_[1]};
  print "
Help for $0 
 this script generates half- and full-level text files of level values
 suitable for input to wrfgrib2arl: -H half_name -F full_name
 Options:
";
foreach my $key (keys(%options)) {
  print "  $key : $helpopt{$key}\n";
}
die 
"Usage: $0 [options]\n";
}
