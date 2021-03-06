#!/usr/bin/perl -w

=head1 Bio::TrapCluster::BuildUniTrap

=head2 Authors

=head3 Created by

             Vincenza Maselli
             v.maselli@ucl.ac.uk

=head3 Original script by

              Guglielmo Roma
              guglielmoroma@libero.it

=head2 Description
        
             This module run blast or read information from blastoutput or from database and load the tables `trapmap_region` and `unitrap`	
             
=head2 Usage

	    my $obj = Bio::TrapCluster::BuildUniTrap->new;
            
=cut

package Bio::TrapCluster::BuildUniTrap;

use strict;
use DBI;
use Carp;
use Data::Dumper;
use vars qw(@ISA);
use File::Spec;
use Bio::TrapCluster::Utils::Argument qw(rearrange);
use Bio::TrapCluster::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
require "$ENV{'TrapCluster'}/trapcluster_conf.pl";
@ISA = qw(Bio::Root::Root Bio::TrapCluster::AnnotationTrap);

use Bio::TrapCluster::Fetch;
use Bio::TrapCluster::LoadTrap;
my $fetch = Bio::TrapCluster::Fetch->new;
my $load = Bio::TrapCluster::LoadTrap->new;

sub new{
  	my $caller = shift;

 	my $class = ref($caller) || $caller;
  	my $self = $class->SUPER::new(@_);
	
  	#my ()= rearrange( [],@_);

}


sub build_unitrap{
	my ($self, $toinsert) = @_;
	
	my $last_id = $fetch->get_last_unitrap_id + 1;
	unless ($last_id){$last_id = 1}
	my $uni_accession = $toinsert->{'seq_id'}."-UNI";
	my $num_zero = 9 - (length $last_id);
	for (my $i = 1; $i <$num_zero; $i ++){
		$uni_accession .= "0";
	}
	
	
	$uni_accession.=$last_id."#".$toinsert->{'rank'};
	
	$toinsert->{'accession'} = $uni_accession;

	my $unitrap_id = $load->load_unitrap($toinsert);

	return $unitrap_id;
}



1;