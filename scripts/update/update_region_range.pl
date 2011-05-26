#! /usr/bin/perl -w
=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            

=head2 Usage

           
=head2 CODE BEGIN

=cut

BEGIN{

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use Data::Dumper;

use Bio::Unitrap::Utils::DBSQL::DB;
my $db = Bio::Unitrap::Utils::DBSQL::DB->new;

my $sql = qq{select region_id, seq_id, parent_id from region where region_end < region_start};

my $regions = $db->select_many_from_table($sql);

foreach my $intron (@{$regions}){
	my $region_id = $intron->{'region_id'};
	my ($prename,$postname) = split /_/,$intron->{'seq_id'};
	
	my $presql = qq{SELECT region_end FROM region WHERE seq_id = \"$prename\"};
	my $preregion = $db->select_from_table($presql);
	my $start = $preregion->{'region_end'} + 1;
	
	my $postsql = qq{SELECT region_start FROM region WHERE seq_id = \"$postname\"};
	my $postregion = $db->select_from_table($postsql);
	my $end = $postregion->{'region_start'} - 1;
	
	my $update = qq{UPDATE region SET region_start = $start, region_end = $end WHERE region_id = $region_id};
	
	$db->insert_set($update);
	
	


}