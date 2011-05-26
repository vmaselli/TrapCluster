 #! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script run the code to annotate trap

=head2 Usage

            [-id] [-idcount] || [-file] || [-blastfile] || [-map] [-l] || [-name] || [-h help] 

           
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

use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Collection;

use Bio::Unitrap::AnnotationTrap;
use Bio::Unitrap::LoadTrap;
use Bio::Unitrap::Fetch;
use Bio::Unitrap::MapTrap;

my $map = Bio::Unitrap::MapTrap->new;
my $ann = $map->annotation;
my $load = $ann->load;
my $fetch = $load->fetch;
my %conf =  %::conf;
my $swarm = $conf{'global'}{'swarm'};
my $debug = $conf{'global'}{'debug'};
my ($id, $idcount, $file,$blastfile,$trapname,$m,$l,$help);
my $opt = &GetOptions( 
	    "file|f=s" =>    	\$file,
	    "blastfile|b=s" =>    	\$blastfile,
	    "name|n=s" =>\$trapname,
	    "map|m" =>\$m,
	    "limit|l=s" =>\$l,
	    "id=s" => \$id,
	    "idcount=s" => \$idcount,
	    "help|h" => \$help
	    );


my $USAGE = "./anntrap.pl [-id] [-idcount]  [-file] [-blastfile]  [-map]   [-name] || [-h help] ";

if ($help){print STDERR "$USAGE
                  -id = the starting trap db id
         		  -idcount = the amount of trap to annotate
         		  
         		  -f file : the alignment output complete path
         		  
                  -n name  : the trap to annotate
                  
                  -m map : if you want annotate the mapped trap
                  -h : this help\n";
			exit;
}

unless($file){
	unless($blastfile){
		unless($trapname){
			unless ($map){
				unless ($id || $idcount) {
					print STDERR "Please set options. Type -help for informations\n";
					exit;
				}
			}
		}
	}
}
	    

if ($file){
	open (IN, $file) || die;
	my $timestmp = `date`;
	print STDERR "-- annotate from file $file\n";
	while (my $row = <IN>){
		chomp $row;
		my ($trapname) = $row; 
		my $trap = $fetch->get_trap_by_name($trapname);
		my $trap_id = $trap->{'trap_id'};
		print STDERR "call annotate from file for trap ID $trap_id\n";
		$ann->annotate($trap_id,1); 
	}
}elsif($trapname){
	my $timestmp = `date`;
	print STDERR "-- annotate $trapname\n";
	my $trap = $fetch->get_trap_by_name($trapname);
	my $trap_id = $trap->{'trap_id'};
	print STDERR "call annotate for trap ID $trap_id\n";
	$ann->annotate($trap_id,1); 
}elsif($m){
	my $timestmp = `date`;
	#print STDOUT "113 $timestmp";
	print STDOUT "-- START annotate mapped in anntrap.pl\n";
	my ($trapmap_id, $lastid);
	unless($swarm){
		$trapmap_id = $id;
		$lastid = $trapmap_id + $idcount ;
	}
	else{
		$trapmap_id = ($id * $idcount) - $idcount + 1;
		$lastid = $id * $idcount;
	}
	print STDOUT "-- START WITH TRAPMAP ID $trapmap_id LASTID $lastid in anntrap.pl\n";
	for(my $idx = $trapmap_id; $idx <= $lastid; $idx ++){
		unless ($fetch->test_trapmap_by_id($idx)){
			$debug && print STDOUT "NO trapmap for TRAPMAP ID $idx\n";
			next;
		}
		
		my $test_id = $fetch->test_ambiguous_by_id($idx);
		if ($test_id && $debug){
			print STDOUT "TRAPMAP ID $idx ALREADY ANNOTATED\n";
			#next;
		}
		else{
			my $test = $fetch->get_trapmap_by_id($idx);
			if ($test->{'trapmap_id'} && $test->{'annotation_ambiguous'} == 0 && $test->{'type'} ne 'genomic' && $debug){
				print STDOUT "TRAPMAP ID $idx ALREADY ANNOTATED\n";
				#next;
			}
		}
		
		my $trap_id = $fetch->get_trap_id_by_trapmap_id($idx);
		
		print STDOUT "-- call annotate from mapped for trap ID $trap_id\n";
		my $ok = 0;
		eval{$ok = $ann->annotate($trap_id,1)};
		unless ($ok){
			$timestmp = `date`;
			print STDOUT "-- NOT DONE annotate for trapmap $idx for trap ID $trap_id\n";
			next;
		}
		$timestmp = `date`;
		print STDOUT "-- FINISHED annotate from mapped $idx for trap ID $trap_id\n";
		
	}
	$timestmp = `date`;
	print STDOUT "-- DONE annotate mapped\n";
}elsif($blastfile){
	my $searchio = new Bio::SearchIO->new (-file => $blastfile, -format => 'blast') || die $!;
	my $timestmp = `date`;
	print STDERR "\n-- annotate from blast\n";
	while (my $result = $searchio->next_result) {
		my $trap = $map->_get_trap($result->query_name);
		my $trap_id = $trap->{'trap_id'};
		$timestmp = `date`;
		print STDERR "call annotate from blastfile for trap ID $trap_id\n";
		$ann->annotate($trap_id,1);
	}
}else{
	my $timestmp = `date`;
	print STDERR "-- annotate from $id to ".($id + $idcount -1 )."\n";
	$ann->annotate($id, $idcount);
}



1;

