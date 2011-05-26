#! /usr/bin/perl

BEGIN{

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use Bio::Unitrap::Utils::DBSQL::DB;
my $db = Bio::Unitrap::Utils::DBSQL::DB->new;

my $file = shift @ARGV;
open(IN, $file);
my $fileout2 = "notmapped.txt";
open(OUT2, ">>$fileout2");
my $fileout1 = "mapped_but_not_annotated.txt";
open(OUT1, ">>$fileout1");

while (my $row = <IN>){
	chomp $row;
	my @fields = split /\t/,$row;
	my $trap_name = $fields[5];
	next if $fields[6] ne "NULL";
	
	my $sql1 = qq{select trapcheck.checkid from trapcheck, trap  where trapcheck.trap_id = trap.trap_id and trap.trap_name = \'$trap_name\' and trapcheck.mapped = 1 and trapcheck.annotated = 0};	
	my $check1 = $db->select_many_from_table($sql1);
	print "$sql1\n";	
	if (scalar @{$check1}){
		print OUT1 "$trap_name\tmapped but not annotated\n";
	}
	else{
		my $sql2 = qq{select trapcheck.checkid from trapcheck, trap  where trapcheck.trap_id = trap.trap_id and trap.trap_name = \'$trap_name\' and trapcheck.mapped = 0};
        	my $check2 = $db->select_many_from_table($sql2);

		print "$sql2\n";	
        	if (scalar @{$check2}){
                	print OUT2 "$trap_name\tnotmapped\n";
        	}	
	}

}
