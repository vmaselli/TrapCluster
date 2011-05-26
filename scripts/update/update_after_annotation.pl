#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script update the trap status

=head2 Usage

            

           
=head2 CODE BEGIN

=cut

BEGIN{

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}
use DBI;
use Getopt::Long;
use Data::Dumper;
use File::Spec;

use Bio::Unitrap::AnnotationTrap;
use Bio::Unitrap::LoadTrap;
use Bio::Unitrap::Fetch;

my $load = Bio::Unitrap::LoadTrap->new;
my $fetch = Bio::Unitrap::Fetch->new;

my %conf =  %::conf;

# Update splinkerette pair

my $select = qq{SELECT e.project_id,t.trap_id, t.trap_name FROM trap t, trapcheck tc, esclone e WHERE e.esclone_id = t.esclone_id AND tc.trap_id=t.trap_id AND tc.splk = 1 ORDER BY t.trap_name ASC};

my $update = qq{UPDATE trap SET paired_tag_id = $paired_id WHERE trap_id = $trap_id};


my $common;
if ($project_id == 4){$common = substr($trap_name, 2);}
elsif ($project_id == 13){($common) = split /\./,$trap_name;}
elsif ($project_id == 12){$common = substr($trap_name,0,-2)}
else{$common = $trap_name}

