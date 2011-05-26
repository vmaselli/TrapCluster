#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script retrieves the data from the source indicated into the config file and loads the database for trap or writes the file as indicated into the config file

=head2 Usage

            ./loadtrap.pl [-id] [-idcount]  [-file] [-type] || [-h help] 

           
=head2 CODE BEGIN

=cut

BEGIN{

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use Data::Dumper;

use Bio::Unitrap::LoadTrap;
use Bio::Unitrap::RetrieveTrap;
use Getopt::Long;

my %conf =  %::conf;
my $tmpdir = "/tmp";
my $swarm = $conf{'global'}{'swarm'};
my $debug = $conf{'global'}{'debug'};
my $dir = $conf{'default'}{'tmp_dir'};
my ($id, $idcount, $file,$type,$help);
my $opt = &GetOptions( 
	    "file|f=s" =>    	\$file,
	    "type|t=s" =>\$type,
	    "id=s" => \$id,
	    "idcount=s" => \$idcount,
	    "help|h" => \$help
	    );


my $USAGE = "./loadtrap.pl [-id] [-idcount]  [-file] [-type] || [-h help] ";

if ($help){print STDERR "$USAGE
                  -id = the starting trap db id
         		  -idcount = the amount of trap to annotate
         		  
         		  -f file : the alignment output complete path
         		  
                  -t type  : file format 
                  
                  -h : this help\n";
			exit;
}


=pod

	1. RetrieveTrap module. This module reads information from file or web or local database and load the tables: `project`, `esclone`, and `trap` or writes the info into files (not implemented yet).

	If the config file is correctly read the script will procede as indicated in it.
	If from file the dir containing one or more files and the file format are specified into the confing file,as well as the db name of from another database or the query to submit to genbank.
	If the $check variable is defined it means that only an update is necessary, if it is not defined, it means that it is the first database loading.

=cut

$debug && print STDOUT "loadtrap.pl line 78\n";
my $retrieve = Bio::Unitrap::RetrieveTrap->new(-dir =>"/tmp");
$debug && print STDOUT "loadtrap.pl line 80\n";
my $load = Bio::Unitrap::LoadTrap->new(-fetch => $retrieve->fetch);
$debug && print STDOUT "loadtrap.pl line 82\n";
my $check = $retrieve->update;
$debug && print STDOUT "loadtrap.pl line 84\n";

=pod

	load to database or file (not implemented yet)
	
=cut

if ($retrieve->todb){
	$debug && print STDOUT "loadtrap.pl line 93\n";
	my $traps;
	if ($retrieve->fromfile){
		$debug && print STDOUT "loadtrap.pl line 96\n";
		my $version;
		unless ($type){
			print STDOUT "Please provide a file type, choose from the following GFF (default), FASTA, GenBank\n";
			$type = <STDIN>;
			$type = "gff" unless $type;
		}
		if ($type eq "gff"){
			unless ($file){
				print STDOUT "Please provide a GFF file name in this format name.gffX where X is the gff format. The gff will split automatically in several files in order to make the query on GenBank\n";
				$file = <STDIN>;
				die unless $file;
			}
			$version = substr($file,-1);
			if ($version ne 1 && $version ne 2 && $version ne 3){
				print STDOUT "Please provide a GFF file name in this format name.gffX where X is the gff format. The gff will split automatically in several files in order to make the query on GenBank\n";
				$file = <STDIN>;
				die unless $file;
			}
		}
		my $wc = `wc -l $file`;
		my ($nro) = split / /, $wc;
		my $count = int($nro/100) + 1;
		if ($swarm){
			unless ($id || $idcount){
				print STDOUT "Please provide from which row you want to start\n";
				$id = <STDIN>;
				print STDOUT "Please, now provide how manu row do you want to read (default $count)\n";
				$idcount = <STDIN>;
				$idcount = $count unless $idcount;
				die unless $id;
			}
		}
		$debug && print STDOUT "loadtrap.pl line 129\n";
		my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => $version);
	
		# loop over the input stream
		my $countid = 1;
		my $countidcount = 1;
		my $newgff;
		$debug && print STDOUT "ID $id, idcount $idcount\n";
		
		while( my $feature = $gffio->next_feature()) {
			my @tags = $feature->get_all_tags;
			my ($trap_name) = $feature->get_tag_values('Sequence_tag_ID');
			my $test_id = $retrieve->fetch->get_trap_id_by_name($trap_name);
			next if $test_id;
			if ($countid < $id){$countid ++;next;}
			else{
				if ($countidcount <= $idcount){$countidcount ++}
				else{last}
			}
			
# 			my $newgff = File::Spec->catfile($dir, $trap_name.".gff".$version);
#  			my $gwriter = new Bio::Tools::GFF(  -gff_version => $version,
#  												-file        => ">$newgff");
# 			$debug && print STDOUT "===> WORK ON $dir\t$countid\t$countidcount\t$newgff\n";	
# 			$gwriter->write_feature($feature);
			#$traps = $retrieve->get_from_file();
			$traps = $retrieve->get_from_gff_feature($feature);
#			`rm $newgff`;
		}
		$debug && print STDOUT "loadtrap.pl line 151\n";
	}
	
	
	elsif ($retrieve->fromdb){$traps = $retrieve->get_from_db($conf{'retrieve'}{'dbname'})}
	elsif ($retrieve->fromquery){$traps = $retrieve->get_from_query()}

	$load->load_db($traps);
	
}

elsif ($retrieve->tofile){
	print STDOUT "Sorry, writing into file is not implemented yet.\n";
	exit;
}


