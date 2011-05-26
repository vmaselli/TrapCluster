#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script return file in different format

=head2 Usage

            ./get_output.pl

           
=head2 CODE BEGIN

=cut

BEGIN{

  print STDERR "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}

use strict;
use vars;
use Data::Dumper;
use Getopt::Long;
use Bio::Unitrap::Utils::TextUtility;
use Bio::Unitrap::Utils::File;
use Bio::Unitrap::Fetch;
my $fetch = Bio::Unitrap::Fetch->new;

my ($tag,$file,$help);
my $opt = &GetOptions( 
	    "tag|t=s" =>    	\$tag,
	    "help|h" => \$help
	    );


my $USAGE = "./get_output.pl [-t] || [-h help] ";

if ($help){print STDERR "$USAGE
                  -t tag : the type of output
                  			gff
                  			bed
                  			fasta (for all fasta)
                  			fasta2 (for the trap indicated into a file)
                  -f file   input file with trap
                  -h : this help\n";
			exit;
}


if ($tag eq "gff"){Bio::Unitrap::Utils::TextUtility->get_unitrap_gff}
elsif ($tag eq "fasta"){Bio::Unitrap::Utils::TextUtility->get_all_trap_fasta}
elsif ($tag eq "bed"){Bio::Unitrap::Utils::TextUtility->get_unitrap_bed}
elsif ($tag eq"fasta2"){
	open (IN, $file) || die $!;
	my $count = 1;
	my %seen;
	my $i = 1;
	while (<IN>){
		chomp $_;
		my ($trap_name) = split /\t/, $_;
		next if $seen{$trap_name};
		print "======== $count $trap_name ==========\n";
		$seen{$trap_name} ++;
		my $trap = $fetch->get_trap_by_name($trap_name);
		my $trap_id = $trap->{'trap_id'};
		my $file = Bio::Unitrap::Utils::TextUtility->get_trap_by_id_fasta($trap_id);
		$file =~ s/\/tmp//;
		my $data = Bio::Unitrap::Utils::File->readFileData($file,"/tmp/");
		if ($count == 50){$i ++, $count = 0;}
		if ($count <= 50){
			my $fastaname = $i."_notmapped.fa";
        		Bio::Unitrap::Utils::File->appendDataToFile ($fastaname,"./",$data);
			print "written file $fastaname\n";
		}
		$count ++;
	}




}
else{Bio::Unitrap::Utils::TextUtility->get_all_trap_by_project_fasta($tag)}
