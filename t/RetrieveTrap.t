# !/usr/bin/perl

use strict;
use vars;
use Data::Dumper;
use Bio::Unitrap::RetrieveTrap;
use Bio::Unitrap::Utils::DBSQL::DB;

BEGIN{

  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  print "Test begins\n";
  
}


my %conf =  %::conf;
my $tmpdir = $conf{'default'}{'tmp_dir'};
my $dir = "t/data/";
my $ok = 0;
my $totsub = `grep -c 'sub' $ENV{'Unitrap'}/modules/Bio/Unitrap/RetrieveTrap.pm`;
my $totcomm = `grep -c '#sub' $ENV{'Unitrap'}/modules/Bio/Unitrap/RetrieveTrap.pm`;
my $totprivate = `grep -c "sub _" $ENV{'Unitrap'}/modules/Bio/Unitrap/RetrieveTrap.pm`;
#my $tot = $totsub - $totcomm - $totprivate;
my $tot =17;
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

if (defined $obj){$ok ++; print "ok new ";}
else{print "not ok new ";}
print "test $ok/$tot\n";  

if(defined $obj->mindate){$ok ++; print "ok mindate ";}
else{print "not ok mindate ";}
print "test $ok/$tot\n";

if(defined $obj->update){$ok ++; print "ok update ";}
else{print "not ok update ";}
print "test $ok/$tot\n";

if(defined $obj->dir){$ok ++; print "ok dir ";}
else{print "not ok dir ";}
print "test $ok/$tot\n";

if(defined $obj->formatfile){$ok ++; print "ok format ";}
else{print "not ok format ";}
print "test $ok/$tot\n";

if(defined $obj->swarm){$ok ++; print "ok swarm ";}
else{print "not ok swarm ";}
print "test $ok/$tot\n";

if(defined $obj->fromfile){$ok ++; print "ok fromfile ";}
else{print "not ok fromfile ";}
print "test $ok/$tot\n";

if(defined $obj->fromquery){$ok ++; print "ok fromquery ";}
else{print "not ok fromquery ";}
print "test $ok/$tot\n";

if(defined $obj->fromdb){$ok ++; print "ok fromdb ";}
else{print "not ok fromdb ";}
print "test $ok/$tot\n";


if(defined $obj->todb){$ok ++; print "ok todb ";}
else{print "not ok fromdb ";}
print "test $ok/$tot
";

if(defined $obj->togff){$ok ++; print "ok togff ";}
else{print "not ok fromdb ";}
print "test $ok/$tot
";

if(defined $obj->tofasta){$ok ++; print "ok tofasta ";}
else{print "not ok fromdb ";}
print "test $ok/$tot
";

if(defined $obj->dir){$ok ++; print "ok dir ";}
else{print "not ok fromdb ";}
print "test $ok/$tot
";

# my $traps = $obj->get_from_db($conf{'dbaccess'}{'testdbname'}, $conf{'dbaccess'}{'testhost'}, $conf{'dbaccess'}{'testuser'}, $conf{'dbaccess'}{'testport'}, $conf{'dbaccess'}{'testpasswd'});
# my ($test) = @{$traps};
# if (defined $test->{'name'} && defined $test->{'project_name'} ){$ok ++; print "ok get_from_db ";}
# else{print "not ok fromdb ";}
# print "test $ok/$tot";
# 
 my $traps2 = $obj->get_from_query;
 my ($test2) = @{$traps2};
 if (defined $test2->{'name'} && defined $test2->{'project_name'} ){$ok ++; print "ok get_from_query ";}
 else{print "not ok fromquery ";}
 print "test $ok/$tot
# ";

my $traps3 = $obj->get_from_file();
my ($test3) = @{$traps3};
if (defined $test3->{'name'} && defined $test3->{'project_name'} ){$ok ++; print "ok get_from_file ";}
else{print "not ok get_from_file ";}
print "test $ok/$tot\n";

if ($ok == $tot){print "All test succeed\n";}
else{print $tot - $ok, " test failed\n";}

END {
	$db_obj->drop_db;
}