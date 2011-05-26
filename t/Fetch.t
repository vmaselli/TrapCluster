# !/usr/bin/perl

use strict;
use vars;
use Data::Dumper;


BEGIN{

  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  print "Test begins\n";
  
}

use Bio::Unitrap::Fetch;
my %conf =  %::conf;
my $tmpdir = $conf{'default'}{'tmp_dir'};
my $ok = 0;
my $totsub = `grep -c 'sub' $ENV{'HOME'}/src/UniTrap/modules/Bio/Unitrap/Fetch.pm`;
my $totcomm = `grep -c '#sub' $ENV{'HOME'}/src/UniTrap/modules/Bio/Unitrap/Fetch.pm`;
my $totprivate = `grep -c "sub _" $ENV{'HOME'}/src/UniTrap/modules/Bio/Unitrap/Fetch.pm`;
my $tot = $totsub - $totcomm - $totprivate;

my $obj = Bio::Unitrap::Fetch->new;

if (defined $obj){$ok ++; print "ok new ";}
else{print "not ok new ";}
print " test $ok/$tot\n";  	

my $newdbID;
my $pr_dbID = $obj->get_project_id_by_name('cmhd');
if (defined $pr_dbID){$ok ++; print "ok get_project_id_by_name $pr_dbID";}
else{
	my $pr_insert = qq{INSERT INTO project SET project_name = "cmhd", url="http://www.cmhd.ca/", wordkey="cmhd"};
	$newdbID =$obj->store($pr_insert);
	$pr_dbID = $obj->get_project_id_by_name('cmhd');
	if (defined $pr_dbID){$ok ++; print "ok get_project_id_by_name $pr_dbID";}
	else{print "not ok get_project_id_by_name ";}
}
print " test $ok/$tot\n";
if (defined $newdbID){$ok ++; print "ok store $newdbID";}
else{print "not ok store ";}
print " test $ok/$tot\n";

my $ln_dbID = $obj->get_logic_name_id_by_name('test');
if (defined $ln_dbID){$ok ++; print "ok get_logic_name_id_by_name $ln_dbID";}
else{
	my $ln_insert = qq{INSERT INTO logic_name SET name = "test"};
	$obj->store($ln_insert);
	$ln_dbID = $obj->get_logic_name_id_by_name('test');
	if (defined $ln_dbID){$ok ++; print "ok get_logic_name_id_by_name $ln_dbID";}
	else{print "not ok get_logic_name_id_by_name ";}
}
print " test $ok/$tot\n";

my $esc_dbID = $obj->get_esclone_id_by_clone('test');
if (defined $esc_dbID){$ok ++; print "ok get_esclone_id_by_clone $esc_dbID";}
else{
	my $esc_insert = qq{INSERT INTO ESClone SET clone = "test", project_id = $pr_dbID};
	$obj->store($esc_insert);
	$esc_dbID = $obj->get_esclone_id_by_clone('test');
	if (defined $esc_dbID){$ok ++; print "ok get_esclone_id_by_clone $esc_dbID";}
	else{print "not ok get_esclone_id_by_clone ";}
}
print " test $ok/$tot\n";

if ($ok == $tot){print "All test succeed\n";}
else{print $tot - $ok, " test failed\n";}