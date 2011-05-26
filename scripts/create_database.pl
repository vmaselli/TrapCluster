#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script create the database if not exists and load the sql schema into the database

=head2 Usage

            ./create_database.pl

           
=head2 CODE BEGIN

=cut

BEGIN{

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

my %conf =  %::conf;

use strict;
use vars;
use Data::Dumper;
use Bio::Unitrap::Utils::DBSQL::DB;

my $db_obj = Bio::Unitrap::Utils::DBSQL::DB->new(-dbname=>'mysql');
$db_obj->dbname($conf{'dbaccess'}{'dbname'});
if (defined $db_obj->create_db){print "Created Database\n ";}
else{print "Could not create database";}

my $sql_dir = "$ENV{'Unitrap'}/sql/";
my $schema = $sql_dir."unitrap.sql";
if ($db_obj->exec_import($schema)){ print "Imported Schema\n";}
else{
	print "Could not import schema";
}

