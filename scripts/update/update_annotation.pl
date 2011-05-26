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
use Getopt::Long;

use Bio::Unitrap::AnnotationTrap;
use Bio::Unitrap::Fetch;
use Bio::Unitrap::LoadTrap;
my $load = Bio::Unitrap::LoadTrap->new;
my $ann = Bio::Unitrap::AnnotationTrap->new;
my $fetch = Bio::Unitrap::Fetch->new;


my $sql = qq{select i.insertion_id, i.trapblock_id, t.sequencing_direction, tba.region_id, tba.flanking_exon_id, tba.intronic from insertion i, trapblock_annotation tba, trap t, region r, region r2 where i.trap_id = t.trap_id AND r.region_id =i.gene_id and r2.region_id = tba.region_id and r2.parent_id=i.gene_id and i.trapblock_id = tba.trapblock_id and tba.intronic = 2 and tba.region_id != tba.flanking_exon_id and tba.display_name not like "ENSMUSE%_ENSMUSE%" order by i.trap_id, i.insertion_id ASC};


foreach my $res (@{$fetch->select_many_from_table($sql)}){
	
		my $hash_ref = $ann->find_insertion_site($res);
		my $insID = $res->{'insertion_id'};
		my $tbID = $res->{'trapblock_id'};
		
		my $pis = $hash_ref->{'insertion'}{'putative_insertion_start'};
		my $pie = $hash_ref->{'insertion'}{'putative_insertion_end'};
		my $mycase = $hash_ref->{'insertion'}{'insertion_case'};
		my $ic = $load->load_logic_name($mycase);
		my $name = $hash_ref->{'annotation'}{'display_name'};
		
		my $update1 = qq{update insertion set putative_insertion_start = $pis, putative_insertion_end = $pie, insertion_case = $ic where insertion_id = $insID};
		print STDOUT "$update1\n";		
		$fetch->update($update1);
		my $update2 =qq{update trapblock_annotation set display_name = \"$name\" where trapblock_id = $tbID};
		print STDOUT "$update2\n";
		$fetch->update($update2);


}

