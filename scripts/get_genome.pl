#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script load and format the mouse genome 

=head2 Usage

            ./get_genome.pl -i [index] -d [dir] -rm -all
           
=head2 CODE BEGIN

=cut

BEGIN{

  print STDERR "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use File::Spec;
use Bio::Unitrap::Utils::File;
use Getopt::Long;
my %conf =  %::conf;

my ($index, $dir, $all,$rm,$help);
my $opt = &GetOptions( 
	    "index|i=s" =>    	\$index,
	    "dir|d=s" => \$dir,
	    "all|a" => \$all,
	    "rm|r=s" => \$rm,
	    "help|h" => \$help
	    );


my $USAGE = "./get_genome.pl -i [index] -d [dir] -rm -all\n";

if ($help){print "$USAGE
                  -i index : the format for the genome (blast, blat)
                  -a all:  1 only chromosome, 2 for NT only, 3 complete
                  -rm: 1 repeat masked, 0 not masked
                  -d dir : the db directory if different from the setting
                  -h : this help\n";
			exit;
}

if (!$index || !$all){print "$USAGE. -help for details\n";exit;}
if ($index =~ /blat/){print "Sorry not implemented yet for blat\n";exit;}

my $wget_base = "ftp://ftp.ensembl.org/pub/release-".$conf{'default'}{'ensembl_version'}."/fasta/mus_musculus/dna";
my $base_name = "Mus_musculus.".$conf{'default'}{'hit_db'}.".".$conf{'default'}{'ensembl_version'};

my $fasta_name = lc($base_name);

if ($rm == 1){$base_name .= ".dna_rm";}
if ($rm == 0){$base_name .= ".dna";}
print "you choose RM $rm ALL $all INDEX $index DIR $dir\n";
if ($all == 1 || $all == 3){
	for (my $i = 1; $i <= 19; $i++){
		my $name = $base_name.".chromosome.".$i.".fa.gz";
		my $wget_cmd = $wget_base."/".$name;
		`wget  $wget_cmd `;
	}
	for my $c ("MT", "X", "Y"){
		my $name = $base_name.".chromosome.".$c.".fa.gz";
		my $wget_cmd = $wget_base."/".$name;
		`wget  $wget_cmd `;
	}
}
if ($all == 2 || $all == 32){
my $name = $base_name.".nonchromosomal.fa.gz";
my $wget_cmd = $wget_base."/".$name;
`wget  $wget_cmd `;
}

$fasta_name .= ".chromosomal.fa" if $all == 1;
$fasta_name .= ".nonchromosomal.fa" if $all == 2;
$fasta_name .= ".all.fa" if $all == 3;

my $indexdir;
$indexdir = $ENV{'Blastdir'} if $ENV{'Blastdir'} && $index =~ /blast/;
$indexdir = $ENV{'Blatdir'} if $ENV{'Blatdir'} && $index =~ /blat/;
$indexdir = $dir if defined $dir;

unless (defined $indexdir){
	print "please enter a destination directory for the genome\n";
	$indexdir = <STDIN>;
}

opendir(DIR, $indexdir) || die "$! ";

while (my $file = readdir(DIR)){
	next if $file =~ /^\./;
	next unless $file =~ /gz$/;
	my $file_gz = File::Spec->catfile("./",$file);
	print "gunzip $file_gz\n";
	`gunzip $file_gz`;
	$file =~ s/\.gz//;
	print "add $file\n";
	my $data = Bio::Unitrap::Utils::File->readFileData($file,"./");
	Bio::Unitrap::Utils::File->appendDataToFile ($fasta_name,$indexdir,$data);
	print "rm $file $file_gz\n";
	`rm -rf $file $file_gz`;
}

my $fasta = File::Spec->catfile($indexdir,$fasta_name);
open (IN, $fasta) || die $!;
close (IN);
if ($index =~ /blast/){`formatdb -i $fasta -pF -oT `;}

`rm -rf $fasta`;
