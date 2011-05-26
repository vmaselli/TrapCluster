#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by

             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
=head2 Usage


=head2 Options

            
           
=head2 CODE BEGIN

=cut

$|=1;
BEGIN{

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use Data::Dumper;
use Bio::Unitrap::Fetch;
use Bio::Unitrap::LoadTrap;
use Bio::Unitrap::AnnotationTrap;
my $fetch = Bio::Unitrap::Fetch->new;
my $load = Bio::Unitrap::LoadTrap->new;
my $ann = Bio::Unitrap::AnnotationTrap->new;
my $hit_id = shift @ARGV;
my $sql = qq{select * from trapmap where strand = "." and hit_id like \"$hit_id%\"};

foreach my $item (@{$fetch->select_many_from_table($sql)}){
	my $hash;
	my $region_hash;
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
	print STDOUT "$date\nCHR $chr ";
	my $sel = qq{select region.* from region  where region.region_name = \"$chr\" and region.description = 'genomic' and region.region_start = $start and region.region_end = $end};
	my $resarray = $fetch->select_many_from_table($sel);
	my $slice_core_adaptor = $ann->slicecoreadpt;
	my $region;
	if ($chr =~ /NT/) {$region = "supercontig";} 
	else {$region = "chromosome";}
	my $slice = $slice_core_adaptor->fetch_by_region($region,$chr, $start, $end);
	my $strand = $slice->strand;
	my $rew_strand = $strand;
	my $fw_strand = "1"; 
	$fw_strand = "-1" if $strand eq "1"; 
	
	my $updateF = qq{update trap, trapmap, trapblock, region, trapmap_region  set trapmap.strand = trapblock.strand where trap.trap_id = trapmap.trap_id and trapmap.trapmap_id = trapblock.trapmap_id and trapmap_region.region_id = region.region_id and trapmap_region.trapmap_id = trapmap.trapmap_id and region.region_strand = \"$strand\" and trapblock.strand = \"$fw_strand\" and trap.sequencing_direction = '5' and trapmap.trapmap_id = $id};
	$fetch->update($updateF);
	my $updateR = qq{update trap, trapmap, trapblock, region, trapmap_region  set trapmap.strand = trapblock.strand where trap.trap_id = trapmap.trap_id and trapmap.trapmap_id = trapblock.trapmap_id and trapmap_region.region_id = region.region_id and trapmap_region.trapmap_id = trapmap.trapmap_id and region.region_strand = \"$strand\" and trapblock.strand = \"$rew_strand\" and trap.sequencing_direction = '3' and trapmap.trapmap_id = $id};
	$fetch->update($updateR);
		
	unless (scalar @{$resarray}){
		
		$region_hash->{'region'}{'name'} = $slice->seq_region_name;
		$region_hash->{'region'}{'seq_id'} = $slice->seq_region_name;
		$region_hash->{'region'}{'strand'} = $slice->strand;
		$region_hash->{'region'}{'start'} = $slice->start;
		$region_hash->{'region'}{'end'} = $slice->end;
		$region_hash->{'region'}{'refseq'} = "na";
		$region_hash->{'region'}{'parent_id'} = 0;
		$region_hash->{'region'}{'rank'} = 1;
		$region_hash->{'region'}{'description'} = "genomic";
		$region_hash->{'region'}{'biotype'} = "genomic";
		my $region_id = $load->load_region($region_hash);
		$hash->{'trapmap_region'}{'region_id'} = $region_id;
		$load->load_trapmap_region($hash);
		
		
	}
	
	foreach my $res (@{$resarray}){
		my $region_id = $res->{'region_id'};
		$hash->{'trapmap_region'}{'region_id'} = $region_id;
		$load->load_trapmap_region($hash);
		my $update = qq{update region set region_strand = $strand where region_id = $region_id};
		$fetch->update($update);
	}
	$date = `date`;
	print STDOUT "done!\n$date\n\n";
}
