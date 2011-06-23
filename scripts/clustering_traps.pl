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

use Bio::TrapCluster::Fetch;
use Bio::TrapCluster::ClusterTrap;

my %conf =  %::conf;
my $swarm = $conf{'global'}{'swarm'};
my $debug = $conf{'global'}{'debug'};

my $fetch = Bio::TrapCluster::Fetch->new;

my ($mol_type,$id,$help);
my $opt = &GetOptions( 
	    "mol_type|m=s" =>    	\$mol_type,
	    "id|c=s" => \$id,
	    "help|h" => \$help
	    );

my $usage = "perl clustring_traps.pl -m mRNA -c 1\nChr X c= 20\tChr Y c = 21\t Mt c = 22";

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

if ($id == 20){$id = "X"};
if ($id == 21){$id = "Y"};
if ($id == 22){$id = "M"};

my $chr = "Chr".$id;


print "<$chr>\n";

my $region = "chromosome"; 
my $version = "NCBIM37";
	
my $trapmap_fwd = qq{select tm.* from trapmap tm, trap t, trapblock tb where t.mol_type = '$mol_type' and tm.hit_id = '$chr' and tm.trap_id = t.trap_id and tm.trapmap_id = tb.trapmap_id and tm.chosen = '1' and tm.strand = "1" order by tm.start};
my $trapmap_rev = qq{select tm.start, tm.end from trapmap tm, trap t, trapblock tb where t.mol_type = '$mol_type' and tm.hit_id = '$chr' and tm.trap_id = t.trap_id and tm.trapmap_id = tb.trapmap_id and tm.chosen = '1' and tm.strand = "-1" order by tm.start};

my $fwd = $fetch->select_many_from_table($trapmap_fwd);
my $rev = $fetch->select_many_from_table($trapmap_rev);

print "RUN ClusterTrap on forward strand\n";
my $cluster_fwd = Bio::TrapCluster::ClusterTrap->new($fwd);
$cluster_fwd->run($chr, $version,$region,1);
print "RUN ClusterTrap on reverse strand\n";
my $cluster_rev = Bio::TrapCluster::ClusterTrap->new($rev);
$cluster_rev->run($chr, $version,$region,-1);

