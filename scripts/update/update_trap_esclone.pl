#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script load the database for trap

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
use Bio::Unitrap::Utils::DBSQL::DB;

my $db = Bio::Unitrap::Utils::DBSQL::DB->new;

my $file = shift @ARGV;
open(IN,$file);
while (my $row = <IN>){
	chomp $row;
	my ($trap_name,$esclone) = split /;/,$row;
	my $sql = qq{UPDATE trap SET esclone_id = $esclone WHERE trap_name = \"$trap_name\"};
	my $sth = $db->prepare_stmt($sql);
	$sth->execute;
}