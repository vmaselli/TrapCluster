 #! /usr/bin/perl -w 

=head2 Authors

=head3 Created by

	    Guglielmo Roma
	    guglielm.roma@gmail.com

=head3 Modified by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script run the code to annotate trap

=head2 Usage

           

           
=head2 CODE BEGIN

=cut

BEGIN{

  print "Reading settings from $ENV{'TrapCluster'}/trapcluster_conf.pl\n";
  require "$ENV{'TrapCluster'}/trapcluster_conf.pl";
  
}

use strict;
use vars;
use Data::Dumper;
use Getopt::Long;

use Bio::Root::IO;

use Bio::TrapCluster::ClusterTrap;
use Bio::TrapCluster::Fetch;
my %conf =  %::conf;
my $swarm = $conf{'global'}{'swarm'};
my $debug = $conf{'global'}{'debug'};

my $fetch = Bio::TrapCluster::Fetch->new;

my ($mol_type,$help);
my $opt = &GetOptions( 
	    "mol_type|m=s" =>    	\$mol_type,
	    "help|h" => \$help
	    );

my $usage = "perl clustring_traps.pl -m mRNA";

if ($help){

	print "USAGE: $usage\n";
	exit;
}

unless ($mol_type){
	print "Choose a mol_type between mRNA and Genomic. [mRNA]\n";
	$mol_type = <STDIN>;
	if ($mol_type eq "\n"){$mol_type = 'mRNA';}
	print "Choosed <$mol_type>\n";
}

my $sql = qq{select distinct tm.hit_id from trapmap tm, trap t where t.trap_id = tm.trap_id and tm.chosen = 1 and t.mol_type = 'mRNA' order by tm.hit_id limit 1};

foreach my $res (@{$fetch->select_many_from_table($sql)}){

	my $chr = $res->{'hit_id'};
	my $region = "chromosome"; 
	my $version = "NCBIM37";
		
	my $trapmap_fwd = qq{select tm.* from trapmap tm, trap t, trapblock tb where t.mol_type = '$mol_type' and tm.hit_id = '$chr' and tm.trap_id = t.trap_id and tm.trapmap_id = tb.trapmap_id and tm.chosen = '1' and tm.strand = "1" order by tm.start limit 10};
	my $trapmap_rev = qq{select tm.start, tm.end from trapmap tm, trap t, trapblock tb where t.mol_type = '$mol_type' and tm.hit_id = '$chr' and tm.trap_id = t.trap_id and tm.trapmap_id = tb.trapmap_id and tm.chosen = '1' and tm.strand = "-1" order by tm.start};
	print $trapmap_fwd,"\n";
	
	my $fwd = $fetch->select_many_from_table($trapmap_fwd);
	#my $rev = $fetch->select_many_from_table($trapmap_rev);
	
	print "RUN ClusterTrap on forward strand\n";
	my $cluster_fwd = Bio::TrapCluster::ClusterTrap->new($fwd);
	$cluster_fwd->run($chr, $version,$region);
	
	#my $cluster_rev = Bio::TrapCluster::ClusterTrap->new($rev);
	#$cluster_rev->run($chr, $version,$region);
}


	



