# !/usr/bin/perl

use strict;
use vars;
use Data::Dumper;
use Bio::Unitrap::Utils::File;
use Bio::Unitrap::Utils::DBSQL::DB;

BEGIN{

  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  print "Test begins\n";
  
}

my %conf =  %::conf;
my $tmpdir = $conf{'default'}{'tmp_dir'};
my $ok = 0;
my $totsub = `grep -c 'sub' $ENV{'Unitrap'}/modules/Bio/Unitrap/Utils/DBSQL/DB.pm`;
my $totcomm = `grep -c '#sub' $ENV{'Unitrap'}/modules/Bio/Unitrap/Utils/DBSQL/DB.pm`;
my $totprivate = `grep -c "sub _" $ENV{'Unitrap'}/modules/Bio/Unitrap/Utils/DBSQL/DB.pm`;
print "ALL SUB $totsub - COMMENTED $totcomm - PRIVATE $totprivate\n";
#my $tot = $totsub - $totcomm - $totprivate;
my $tot = 18;
my $db_obj = Bio::Unitrap::Utils::DBSQL::DB->new(-dbname=>'test');

if (defined $db_obj->db_connection){$ok ++; print "ok db_connection ";}
else{print "not ok db_connection ";}
my $sql_dir = "$ENV{'Unitrap'}/sql/";
my $schema = $sql_dir."unitrap.sql";

#Bio::Unitrap::Utils::File->purge_sql($sql_dir);


my $table_name = "logic_name";
my $table_name_id = $table_name."_id";
my $column_name = "name";
my $value = "testword";
my $primary = $table_name."_id";

my $first_insert = qq{INSERT INTO $table_name SET $column_name  = \"$value\"};

my $sec_value = "second insert";
my $par = qq{INSERT INTO $table_name SET $column_name = \"$sec_value\"};

if (defined $db_obj){$ok ++; print "defined new ";}
else{print "not defined new ";}
print "test $ok/$tot\n";

if (defined $db_obj->db_connection){$ok ++; print "defined db_connection ";}
else{print "not defined db_connection ";}
print "test $ok/$tot\n";

if (defined $db_obj->host){$ok ++; print "defined host ";}
else{print "not defined host ";}
print "test $ok/$tot\n";

if (defined $db_obj->dbname()){$ok ++; print "defined dbname ";}
else{print "not defined dbname ";}
print "test $ok/$tot\n";

if (defined $db_obj->user){$ok ++; print "defined user ";}
else{print "not defined user ";}
print "test $ok/$tot\n";

if (defined $db_obj->port){$ok ++; print "defined port ";}
else{print "not defined port ";}
print "test $ok/$tot\n";

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	      -host => $conf{'dbaccess'}{'enshost'},
	      -user => $conf{'dbaccess'}{'ensuser'},
	      -pass => $conf{'dbaccess'}{'enspass'},
	      -database => $conf{'dbaccess'}{'ensdbname'},
	      -port => $conf{'dbaccess'}{'ensport'}
);
$db_obj->ensembl_connection($registry);
if (defined $registry->get_adaptor( 'Mouse', 'Core', 'Slice' )){$ok ++; print "ok ensembl_connection ";}
else{print "not ok ensembl_connection ";}
print "test $ok/$tot\n";


$db_obj->dbname("unitrap");
if (defined $db_obj->create_db){$ok ++; print "ok create_db ";}
else{print "not ok create_db ";}
print "test $ok/$tot\n";
if ($db_obj->exec_import($schema)){$ok ++; print "ok exec_import ";}
else{
	print "not ok exec_import ";
	}

print "test $ok/$tot\n";
my $dbID;
if($dbID = $db_obj->insert_set($first_insert)){$ok ++;print "ok first insert got dbID = $dbID ";}
else{print "not ok first insert ";}
print "test $ok/$tot\n";

if (defined $db_obj->sth){$ok ++; print "ok sth ";}
else{print "not ok sth ";}
print "test $ok/$tot\n";

my $condition = "$primary=$dbID";

my $uppar = qq{UPDATE $table_name SET $column_name = \"$sec_value\" WHERE $condition};
if ($db_obj->update_set($uppar)){$ok ++; print "ok update_set ";}
else{print "not ok update_set ";}
print "test $ok/$tot\n";

my $condition = "$primary=$dbID";                                                                              
my $selpar = qq{SELECT  $column_name FROM $table_name WHERE $condition};
my $res = $db_obj->select_from_table($selpar);
if ($res->{$column_name}){$ok ++; print "ok select_from_table ";}
else {print "not ok select_from_table $selpar";}
print "test $ok/$tot\n";

my $selpar = qq{SELECT  $column_name FROM $table_name };
my $res = $db_obj->select_many_from_table($selpar);
if (scalar @{$res}){$ok ++; print "ok select_many_from_table ";}
else{print "not ok select_many_from_table ";}
print "test $ok/$tot\n";

if($db_obj->check_return($sec_value, $table_name, $column_name)){$ok ++; print "ok check_return ";}
else{print "not ok check_return ";}
print "test $ok/$tot\n";

my $file = $db_obj->exec_dump(1);
my $read_content = Bio::Unitrap::Utils::File->readFileData($file,$tmpdir,1);

if (length $read_content){$ok ++; print "ok exec_dump " }
else{print "not ok exec_dump "}
print "test $ok/$tot\n";

my $delpar = qq{DELETE FROM $table_name WHERE $condition};
my $deldbID = $db_obj->delete_set($delpar);
if ($deldbID){$ok ++; print "ok delete_from_table ";}
else{print "not ok delete_from_table ";}
print "test $ok/$tot\n";

if ($ok == $tot){print "All test succeed\n";}
else{print $tot - $ok, " test failed\n";}

END {
	$db_obj->drop_db;
}
	

