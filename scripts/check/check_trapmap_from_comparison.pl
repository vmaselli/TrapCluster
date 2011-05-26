#! /usr/bin/perl

BEGIN{

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use Bio::Unitrap::Utils::DBSQL::DB;
use Bio::Unitrap::Fetch;
use Bio::Unitrap::LoadTrap;
my $fetch = Bio::Unitrap::Fetch->new;
my $load = Bio::Unitrap::LoadTrap->new;
my $db = Bio::Unitrap::Utils::DBSQL::DB->new;

my $file = shift @ARGV;
open(IN, $file);

while (my $row = <IN>){
	chomp $row;
	my ($trap_name) = split /\t/, $row;
	my $sql = qq{select trapmap.* from trapmap, trap where trap.trap_id=trapmap.trap_id and trap.trap_name = \"$trap_name\"};

	foreach my $item (@{$fetch->select_many_from_table($sql)}){
		my $hash;
		my $id =  $item->{'trapmap_id'};
		my $chr = $item->{'hit_id'};
		my $start = $item->{'start'};
		my $end = $item->{'end'};
		$hash->{'trapmap_region'}{'trapmap_id'} = $id;
		$hash->{'trapmap_region'}{'annotation_ambiguous'} = 0;	
		$hash->{'trapmap_region'}{'number_trapblocks'}=$item->{'num_hsps'};
		$hash->{'trapmap_region'}{'number_annotated_trapblocks'} = $item->{'num_hsps'};
		$hash->{'trapmap_region'}{'total_number_exons'} = 0;
		$hash->{'trapmap_region'}{'overlap'}=100;
		$hash->{'trapmap_region'}{'type'}='genomic';
		my $date = `date`;
		my $sel = qq{select region.* from region  where region.region_name = \"$chr\" and region.description = 'genomic' and region.region_start = $start and region.region_end = $end};
		print STDOUT "$sql\n$sel\n";
		my $resref = $fetch->select_many_from_table($sel); 
		unless (scalar @{$resref}){
			my $rhash;
			$rhash->{'region'}{'name'} = $chr;
        		$rhash->{'region'}{'seq_id'} = $chr;
        		$rhash->{'region'}{'strand'} = "1";
        		$rhash->{'region'}{'start'} = $start;
        		$rhash->{'region'}{'end'} = $end;
        		$rhash->{'region'}{'refseq'} = "NULL";
        		$rhash->{'region'}{'parent_id'} = 0;
        		$rhash->{'region'}{'rank'} = 1;
        		$rhash->{'region'}{'description'} = 'genomic';
			my $region_id = $load->load_region($rhash);	
			$hash->{'trapmap_region'}{'region_id'} = $region_id;
			$load->load_trapmap_region($hash);
		}
		
		foreach my $res (@{$resref}){
		
			$hash->{'trapmap_region'}{'region_id'} = $res->{'region_id'};
			$load->load_trapmap_region($hash);
		}
		
		$date = `date`;
		print STDOUT "done!\n$date\n\n";
	}	
}
1;
