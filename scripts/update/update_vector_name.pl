#! /usr/bin/perl
$|=1;
BEGIN{

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use Data::Dumper;	
use Bio::Unitrap::Fetch;
my $fetch = Bio::Unitrap::Fetch->new;

open (IN, shift @ARGV) || die $!;
while (<IN>){
	chomp $_;
	$_ =~ s/"//g;
	my @arr = split /,/, $_;
	my $name = $arr[1];
	my $vector = $arr[6];
	
	my $sql = qq{UPDATE esclone SET vector_name = \"$vector\" WHERE clone = \"$name\"};
	
	$fetch->update($sql);

}
