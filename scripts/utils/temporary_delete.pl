#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script report some statistics

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
use Bio::Unitrap::Fetch;
use Bio::Unitrap::Utils::File;

my $fetch = Bio::Unitrap::Fetch->new;
my $dbh = $fetch->db_connection;



my $trapmap = qq{select trapmap_id, trap_id from trapmap};
my $trapmap_sth = $dbh->prepare($trapmap);

my $trapmap_rm = qq{delete tm.* from trapmap tm where tm.trapmap_id = ?};
my $trapmap_rm_sth = $dbh->prepare($trapmap_rm);


$trapmap_sth->execute;


$trapmap_sth->execute();
while (my ($trapmap_id, $trap_id) = $trapmap_sth->fetchrow_array){
	my $test_id = $fetch->test_trap_by_id($trap_id);
	next if $test_id;
	$trapmap_rm_sth->execute($trapmap_id) || die $!;
	print STDOUT "trapmap removed $trapmap_id\n";
}


