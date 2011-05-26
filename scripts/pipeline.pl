#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            Script to manage and run the whole pipeline

=head2 Usage

            ./pipeline.pl

           
=head2 CODE BEGIN

=cut

BEGIN{

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use Data::Dumper;

`$ENV{'Unitrap'}/script/create_databse.pl`;
`$ENV{'Unitrap'}/script/get_genome.pl`;
`$ENV{'Unitrap'}/script/loadtrap.pl`;
`$ENV{'Unitrap'}/script/update_trap_pair.pl`;
`$ENV{'Unitrap'}/script/maptrap.pl`;
`$ENV{'Unitrap'}/script/get_output.pl`;
`$ENV{'Unitrap'}/script/buildunitrap.pl`;
`$ENV{'Unitrap'}/script/unitrap_history.pl`;