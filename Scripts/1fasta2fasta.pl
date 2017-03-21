#!/usr/bin/perl -w

# split and put onto two lines
# e.g. for o/p from uniptoy_process_TMs.pl, for seq processing
# this version separates (I hope) on _TM_ to get back original seq IDs

open (IN, "<process.in") or die;
open (OUT, ">process.out") or die;
while ($line 		= <IN>) {
  chomp $line;
  @words 		= split (" ",$line);
  if (exists $words[0]) {
#    printf OUT "$words[0]\n";
    @wordsB	= split ("_TM_",$words[0]);
    printf OUT "$wordsB[0]\n";
    printf OUT "$words[1]\n";
  }
}
close (IN);
close (OUT);

exit;

