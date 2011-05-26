# !/usr/bin/perl

use strict;
use vars;
use Data::Dumper;

use Bio::Unitrap::LoadTrap;
use Bio::Unitrap::RetrieveTrap;

BEGIN{

  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  print "Test begins\n";
  
}

my %conf =  %::conf;
my $tmpdir = $conf{'default'}{'tmp_dir'};
my $dir = "t/data/";

my $ok = 0;
my $totsub = `grep -c 'sub' $ENV{'Unitrap'}/modules/Bio/Unitrap/LoadTrap.pm`;
my $totcomm = `grep -c '#sub' $ENV{'Unitrap'}/modules/Bio/Unitrap/LoadTrap.pm`;
my $totprivate = `grep -c "sub _" $ENV{'Unitrap'}/modules/Bio/Unitrap/LoadTrap.pm`;
#my $tot = $totsub - $totcomm - $totprivate;
my $tot=3;
my $db_obj = Bio::Unitrap::Utils::DBSQL::DB->new(-dbname=>'test');
$db_obj->dbname("unitrap");
if (defined $db_obj->create_db){$ok ++; print "ok create_db ";}
else{print "not ok create_db ";}

my $sql_dir = "$ENV{'Unitrap'}/sql/";
my $schema = $sql_dir."unitrap.sql";
if ($db_obj->exec_import($schema)){$ok ++; print "ok exec_import ";}
else{
	print "not ok exec_import ";
}

my $obj = Bio::Unitrap::RetrieveTrap->new(-DIR=>$dir);
my $traps = $obj->get_from_file();
my ($test) = @{$traps};

my $loadtrap = Bio::Unitrap::LoadTrap->new(-fetch => $obj->fetch);


if ($loadtrap->load_db($traps)){$ok ++; print "ok load_db ";}
else{print "not ok load_db ";}
print "test $ok/$tot\n";





if ($ok == $tot){print "All test succeed\n";}
else{print $tot - $ok, " test failed\n";}


END {
	$db_obj->drop_db;
}