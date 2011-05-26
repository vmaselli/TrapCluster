#! /usr/bin/perl

BEGIN{

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use Data::Dumper;
use Bio::Unitrap::Utils::DBSQL::DB;
my $db = Bio::Unitrap::Utils::DBSQL::DB->new;

my $file = shift @ARGV;
open(IN, $file);

while (my $row = <IN>){
	chomp $row;
	my @fields = split /\t/,$row;
	my $trap_name = $fields[5];

	my $sql = qq {UPDATE trap, trapcheck SET trapcheck.tested = 1 WHERE trap.trap_id=trapcheck.trap_id AND trap.trap_name = \'$trap_name\';};	
	$db->update_set($sql);

}
