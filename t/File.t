# !/usr/bin/perl

use strict;
use vars;
use Data::Dumper;

BEGIN{

  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}
my %conf =  %::conf;
my $dir = $conf{'default'}{'tmp_dir'};
use Bio::Unitrap::Utils::File;

my $ok = 0;
my $totsub = `grep -c 'sub' $ENV{'Unitrap'}/modules/Bio/Unitrap/Utils/File.pm`;
my $totcomm = `grep -c '#sub' $ENV{'Unitrap'}/modules/Bio/Unitrap/Utils/File.pm`;
my $totprivate = `grep -c "sub _" $ENV{'Unitrap'}/modules/Bio/Unitrap/Utils/File.pm`;
#my $tot = $totsub - $totcomm - $totprivate;
my $tot = 5;

my $data = ">Seq1\nATGGCGATGAGCTAGACTAGTCGATCGATTAGCTAGCTAGCTAGGGATAGCTGGTAAATGC\n";
my $add_data = ">Seq2\nATGGCGATGAGCTAGACTAGTCGATCGATTAGCTAGCTAGCTAGGGATAGCTGGTAAATGC\n";
my $file = "file.t.out";

print "FILE $file DIR $dir\n";

my $file = Bio::Unitrap::Utils::File->writeDataToFile($file,$dir, $data, 1);
if (open(IN,$dir.$file)){$ok ++; print "ok writeDataToFile ";}
else{print "not ok writeDataToFile $file ";}
print "test $ok/$tot\n";
my $read_content = Bio::Unitrap::Utils::File->readFileData($file,$dir,1);
if ($read_content eq $data){$ok ++; print "ok readFileData " }
else{print "not ok readFileData "}
print "test $ok/$tot\n";
my $file = Bio::Unitrap::Utils::File->appendDataToFile($file,$dir, $add_data, 1);
if (Bio::Unitrap::Utils::File->readFileData($file,$dir,1) eq $data.$add_data){$ok ++; print "ok appendFileData "}
else{print "not ok readFileData "}
print "test $ok/$tot\n";
my %hash = Bio::Unitrap::Utils::File->fromFastaFileToHash($file,$dir,1);
if (scalar keys %hash){$ok ++; print "ok fromFastaFileToHash " }
else{print "not ok fromFastaFileToHash "}
print "test $ok/$tot\n";
my $text = Bio::Unitrap::Utils::File->readFileData($file,$dir,1);
my %hash = Bio::Unitrap::Utils::File->fromFastaFormatToHash($text,1);
if (scalar keys %hash){$ok ++; print "ok fromFastaFormatToHash " }
else{print "not ok fromFastaFormatToHash "}
print "test $ok/$tot\n";

if ($ok == $tot){print "All test succeed\n";}
else{print $tot - $ok, " test failed\n";}
