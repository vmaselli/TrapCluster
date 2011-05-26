#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script load the database for trap

=head2 Usage

            ./maptrap.pl [-f file] [-r] [-p] [-dir] [-t] [-id] [-idcount] 

			 -t = blast      BLAST (WUBLAST, NCBIBLAST,bl2seq)   
 				  blasttable BLAST -m9 or -m8 output (both NCBI and WUBLAST tabular)	
  				  blastxml   NCBI BLAST XML
  		     -p blast, blat, exonerate
           
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
use Bio::Unitrap::MapTrap;
use Bio::Unitrap::Utils::TextUtility;
my %conf =  %::conf;
my $swarm = $conf{'global'}{'swarm'};
my $debug = $conf{'global'}{'debug'};
my $tmpdir = $conf{'default'}{'tmp_dir'};

my $map = Bio::Unitrap::MapTrap->new;

my ($file, $id,$idcount,$t,$r, $dir,$program,$help);
my $opt = &GetOptions( 
	    "file|f=s" =>    	\$file,
	    "dir|d=s" => \$dir,
	    "prg|p=s" => \$program,
	    "run|r" => \$r,
	     "id=s" => \$id,
	    "idcount=s" => \$idcount,
	    "type|t=s" => \$t,
	    "help|h" => \$help
	    );


my $USAGE = "./maptrap.pl [-f file] [-r] [-dir] [-p] [-t] [-id] [-idcount]";

if ($help){print "$USAGE
                  -f file : the alignment output complete path or fasta file if run is required
                  -t format: blast format output
                  -p program: the program to use to make the alignment
                  -d dir : the alignment output directory complete path 
                  -r : type if you want to run comparison analysis
                  -id = the index for the job array 
         	  -idcount = the amount of trap to analyze
                  -h : this help\n";
			exit;
}

unless ($program){
	print "Please specifiy the program you want to use to make the alignment. Choose between Blast, Blat or Exonerate [Blast]\n";
	$program = <STDIN>;
	$program = "Blat" unless $program;
}

unless($r){
	unless($file){
		unless($dir){
			print "Please set options. Type -help for informations\n";
			exit;
		}
	}
}else{
	unless ($id){
		print "Please specifiy the trap_id to start with or type enter to use the default [1]\n";
		$id = <STDIN>;
		$id = 1 unless $id;
	}
	unless ($idcount){
		print "Please specifiy how many trap have to be analyzed or type enter to use the default [1]\n";
		$idcount = <STDIN>;
		$idcount = 1 unless $idcount;
	}
}

if ($program =~ /blast/){
	unless ($t){
		print "Please set blast output format\n";
		$t = <STDIN>;
	}
}

if ($r){
	my $timestmp = `date`; 
	unless ($file){
		#write idcount fasta seq in a multifasta
		print "START AT $timestmp get_trap_by_id_range_fasta\n";
		my $trap_id = ($idcount * ($id - 1)) + 1;
		my $lastid = $trap_id + $idcount ;
		print STDOUT "TRAP ID $trap_id LASTID $lastid\n";
		for(my $idx = $trap_id; $idx < $lastid; $idx ++){
		#	my $trap = $traps[$idx];
		#	my $trap_id = $trap->{'trap_id'};
			$file = Bio::Unitrap::Utils::TextUtility->get_trap_by_id_range_fasta($idx, $idcount);
		#	$file is complete path
		}
	}
	else{
	    	my $fasta_id = $id;
		$file = $id."_".$file;
	}
	if (lc($program) =~ /blast/){		
		#run a single blast for idcount trap
		$timestmp = `date`;
		print "START AT $timestmp run_blast on file $file\n";

		my $blastfile = $map->run_blast($file,'fasta');
		$timestmp = `date`;
        print "$timestmp -- finished blast, starting split\n";
        #system("perl $ENV{'Unitrap'}/scripts/splitme.pl $blastfile");
        Bio::Unitrap::Utils::TextUtility->splitme($blastfile);
        $timestmp = `date`;
        print "$timestmp -- finished split, starting parse\n";
        for (my $split_index = 1; $split_index <= 51; $split_index++){
            my $split_file_name = "${blastfile}.${split_index}";
          	my $grep = `grep "Query=" $split_file_name`;
            print "counting hsps in file $split_file_name for query $grep\n";
            my $hsps = count_hsps($split_file_name);
            print "counted $hsps hsps in file $split_file_name\n";
            if($hsps >  50){
              print "skipping\n";
              next;
            }
            $timestmp = `date`;
            print "$timestmp -- running pars_blast on $split_file_name\n";
            my $res = $map->pars_blast_results($split_file_name,$t);
            $timestmp = `date`;
            print "$timestmp -- finished pars_blast on $split_file_name\n" if $res;
            print "$timestmp -- no results on $split_file_name\n" unless $res;
        }

		#`rm $blastfile`;
	}
	elsif (lc($program) =~ /blat/){
		print "Sorry $program run is not still implemented\n";
		exit;
	}
	elsif (lc($program) =~ /exonerate/){
		$timestmp = `date`;
		print "START AT $timestmp run_exonerate\n";
		my $result = $map->run_exonerate($file,'fasta');
		$timestmp = `date`;
		print "START AT $timestmp pars_results\n";
		$map->pars_exonerate_results($result,"exonerate");
	}
	else{
		print "Sorry $program run is not still implemented\n";
		exit;
	}
}

elsif($dir){

	opendir(DIR,$dir) || die "$! $dir";
	while (my $file = readdir(DIR)){
		next if $file =~ /^\./;
		
		my $blastfile = File::Spec->catfile($dir,$file);
		print STDERR "BLASTFILE $blastfile\n";
		$map->pars_results($blastfile,$t);
	}
}
elsif ($file && !$r){
	print STDOUT "FILE $file FORMAT $t\n";
	if ($t =~ /gff/){
		# loop over the input stream
		my $countid = 1;
		my $countidcount = 1;
		my $newgff;
		$debug && print STDOUT "ID $id, idcount $idcount\n";
		my ($row_id, $lastid);
		unless($swarm){
			$row_id = $id;
			$lastid = $row_id + $idcount ;
		}
		else{
			$row_id = ($id * $idcount) - $idcount + 1;
			$lastid = $id * $idcount;
		}
		my $version = substr($t,-1);
		# specify input via -fh or -file
		my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => $version);
		while( my $feature = $gffio->next_feature()) {
			my @tags = $feature->get_all_tags;
			my $location = $feature->location;
			my $start = $location->start;
			my $end = $location->end;
			my $strand = $location->strand;
			my $chr = $feature->seq_id;
	
			my ($trap_name) = $feature->get_tag_values('Sequence_tag_ID');
			my $trap = $map->_get_trap($trap_name);
			
			unless ($trap){`grep $trap_name $file | cat >> missed.gff `; next;}
			unless ($trap->{'mol_type'}){print STDERR "Molecular type error for $trap_name, next\n"; `grep $trap_name $file | cat >> missed.gff `; next;}
			my $trap_id = $trap->{'trap_id'};
			my $trapmap_test_array = $map->load->fetch->get_trapmap_by_trap_id($trap_id);
			my $idtest = 0;
			foreach my $trapmap_test (@{$trapmap_test_array}){
			if ($trapmap_test->{'start'} == $start && $trapmap_test->{'end'} == $end  && $trapmap_test->{'strand'} == $strand && $trapmap_test->{'hit_id'} eq $chr) 
			{
				$idtest = $trapmap_test->{'trapmap_id'};}	
			}
			$debug && print STDOUT "TEST for TRAP ID $trap_id TRAPMAP $idtest in maptrap.pl\n";
			next if $idtest;
			if ($countid < $row_id){$countid ++;next;}
			else{
				if ($countidcount <= $idcount){$countidcount ++}
				else{last}
			}
			$map->parse_gff_feat($feature,$file);
			$countid ++;
			$debug && print STDOUT "ROW $row_id, COUNTID $countid, countidcount $countidcount\n";
		}
	}
	else{
		$map->parse_results($file,$t);
	}
}
	    
### SUBROUTINES

sub count_hsps {
  my $file_name = shift;
  my $count = 0;
  open (BLAST_FILE, "<$file_name") or die "cant open file $file_name\n";
  while(<BLAST_FILE>){
    if($_ =~ /Score =/){
      $count++;
    }
    last if $count >= 100;
  }
  return $count;
}
