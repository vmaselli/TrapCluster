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
my $load = Bio::Unitrap::LoadTrap->new;
my $fetch = Bio::Unitrap::Fetch->new;
my %conf =  %::conf;

my ($file, $r, , $id, $idcount, $help);
my $opt = &GetOptions( 
	    "file|f=s" =>    	\$file,
	    "run|r" => \$r,
	    "id=s" => \$id,
	    "idcount=s" => \$idcount,
	    "help|h" => \$help
	    );


my $USAGE = "./anntrap.pl [-id] [-idcount] [-h help] ";

if ($help){print STDERR "$USAGE
                  -id = the starting trap db id
         		  -idcount = the amount of trap to annotate
                  -h : this help\n";
			exit;
}


unless ($id || $idcount) {
	print STDERR "Please set options. Type -help for informations\n";
	exit;
}

insertion_id  trap_id 	trapmap_id 	trapblock_id 	transcript_id 	putative_insertion_start 	putative_insertion_end 	insertion_case 	insertion_ambiguous 	region_start 	region_end


	my $trap = $fetch->get_trap_by_id($tid);
	
	my @trapmap = @{$fetch->get_trapmap_by_trap_id($trap_id)};

	
	my $trapblock = $fetch->get_trapblock_by_trapmap_id($trapmap_id);
		
	$load->update_insertion($insertionhash);
			


