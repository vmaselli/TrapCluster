#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script load the database for trap

=head2 Usage

            ./anntrap.pl [-f filename] [-r] [-id] [-idcount] [-h help]

           
=head2 CODE BEGIN

=cut

BEGIN{
  print "\nReading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use Data::Dumper;
use Getopt::Long;
use Bio::Unitrap::AnnotationTrap;
use Bio::Unitrap::LoadTrap;
use Bio::Unitrap::Fetch;
use Bio::Unitrap::BuildUniTrap;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Collection;
use Bio::Unitrap::Utils::DBSQL::DB;

my $load = Bio::Unitrap::LoadTrap->new;
my $fetch = Bio::Unitrap::Fetch->new;
my $db = Bio::Unitrap::Utils::DBSQL::DB->new;
my $build = Bio::Unitrap::BuildUniTrap->new;
my %conf =  %::conf;

my ($file, $r, , $id, $a,$idcount, $help);
my $opt = &GetOptions( 
	    "file|f=s" =>    	\$file,
	    "run|r" => \$r,
	    "annotate|a" => \$a,
	    "id=s" => \$id,
	    "idcount=s" => \$idcount,
	    "help|h" => \$help
	    );
	    
	
my %seen;
my %group;
my $rank = 0;
my $debug = $conf{'global'}{'debug'};
print "START query\n";


my $dbh = $fetch->db_connection;

my $insertion = qq{select i.* from insertion i order by i.gene_id, i.putative_insertion_start};
my $isth = $dbh->prepare($insertion);

my $trapmap = qq{select tm.* from trapmap tm, trapmap_region tmr where tm.trapmap_id = tmr.trapmap_id and tm.chosen = 1 and tmr.annotation_ambiguous = 0 and tm.trapmap_id = ? and tmr.region_id = ?};
my $tmsth = $dbh->prepare($trapmap);

my $region = qq{select seq_id from region where region_id = ?};
my $rsth = $dbh->prepare($region);

$isth->execute;

while (my $res = $isth->fetchrow_hashref){
	
	my $unitrap_id;
	
	my $trapblock_id = $res->{'trapblock_id'};
	my $trap_id = $res->{'trap_id'};
	my $trapmap_id = $res->{'trapmap_id'};
	my $start = $res->{'putative_insertion_start'};
	my $end = $res->{'putative_insertion_end'};
	my $gene_id = $res->{'gene_id'};
	$tmsth->execute($trapmap_id, $gene_id);
	my ($tm_res) = $tmsth->fetchrow_hashref;
	next unless defined $tm_res->{'trapmap_id'};
	$rsth->execute($gene_id);
	my $seq_id = $rsth->fetchrow;
	print "START $start END $end\n";
	
	unless ($seen{$start}{$end}){
		$debug && print STDOUT "new unitrap ";
		my $toinsert;
		if ($group{$seq_id}){$rank ++}
		else{$rank = 1}
		$toinsert->{'chr'} = $tm_res->{'hit_id'};
		$toinsert->{'hit_db'} = $tm_res->{'hit_db'};
		$toinsert->{'gene_id'} = $res->{'gene_id'};
		$toinsert->{'start'} = $start;
		$toinsert->{'end'} = $end;
		$toinsert->{'seq_id'} = $seq_id;
		$toinsert->{'rank'} = $rank;
		$unitrap_id = $build->build_unitrap($toinsert);
		$seen{$start}{$end} = $unitrap_id;
		print "rank\n";
	}
	$unitrap_id=$seen{$start}{$end};
	#print "seq_id = $seq_id, unitrap_id = $unitrap_id, rank = $rank\n";
	$group{$seq_id} ++;	
	my $test = qq{SELECT trap_unitrap_id FROM trap_unitrap WHERE trap_id = $trap_id AND unitrap_id = $unitrap_id};
	my $test_id = $db->select_from_table($test);
	next if $test_id->{'trap_unitrap_id'};
	my $insert = qq{INSERT INTO trap_unitrap SET trap_id = $trap_id, trapblock_id = $trapblock_id, unitrap_id = $unitrap_id};
	print "$insert\n";
	$db->insert_set($insert);
	
	#$build->calculate_mutated_protein_foreach_transcript ($toinsert);
	#$build->calculate_primers_for_vector_validation ($toinsert);		

}

print "END loop\n";
