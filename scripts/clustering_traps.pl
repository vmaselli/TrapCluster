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

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use Data::Dumper;
use Getopt::Long;
