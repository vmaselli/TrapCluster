#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            
            This script report some statistics

=head2 Usage

           
           
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
use Bio::Unitrap::Fetch;
use Bio::Unitrap::Utils::File;

my $fetch = Bio::Unitrap::Fetch->new;
my %conf =  %::conf;
my %stat;

my $tot = qq{select count(t.trap_id) as tot, u.accession, p.project_name from trap t, esclone e, project p, unitrap u, trap_unitrap tu where t.esclone_id = e.esclone_id and e.project_id = p.project_id and t.trap_id = tu.trap_id and u.unitrap_id = tu.unitrap_id group by u.accession, p.project_name};

open(TOT,">>complex.tab");

foreach my $totres  (@{$fetch->select_many_from_table($tot)}){
	print TOT "$totres->{'tot'}\t$totres->{'accession'}\t$totres->{'project_name'}\n";
}


my $unigene = qq{select count(unitrap_id) as tot, substr(accession,1,18) as gene from unitrap group by substr(accession,1,18)};

open(UG, ">>unitrap_x_gene.tab");

foreach my $ugres  (@{$fetch->select_many_from_table($unigene)}){
	print UG "$ugres->{'tot'}\t$ugres->{'gene'}\n";
}


my $sql = qq{select p.project_name, t.trap_name, u.accession from trap t, esclone e, project p, unitrap u, trap_unitrap tu where t.esclone_id = e.esclone_id and e.project_id = p.project_id and t.trap_id = tu.trap_id and u.unitrap_id = tu.unitrap_id};

my %seen_unitrap;
my %seen_trap;
my %seen_pair;
my %seen_trap_unitrap;

foreach my $res (@{$fetch->select_many_from_table($sql)}){
	
	$seen_unitrap{$res->{'accession'}} ++;
	$seen_trap{$res->{'trap_name'}} ++;
	unless ($seen_trap_unitrap{$res->{'accession'}}{$res->{'trap_name'}}){
		$seen_trap_unitrap{$res->{'accession'}}{$res->{'trap_name'}} ++;
		push (@{$stat{'unitrap'}{$res->{'accession'}}{'trap'}}, $res->{'accession'});
	}
	unless ($seen_pair{$res->{'project_name'}}{$res->{'accession'}}){
		$seen_pair{$res->{'project_name'}}{$res->{'accession'}} ++;
		push (@{$stat{'project'}{$res->{'project_name'}}{'unitrap'}}, $res->{'accession'});
	}
	unless ($seen_pair{$res->{'project_name'}}{$res->{'trap_name'}}){
		$seen_pair{$res->{'project_name'}}{$res->{'trap_name'}} ++;		
		push (@{$stat{'project'}{$res->{'project_name'}}{'trap'}}, $res->{'trap_name'});
		push (@{$stat{'project'}{$res->{'project_name'}}{'trap_unitrap'}{$res->{'accession'}}}, $res->{'trap_name'});
	}
}

open (TR, ">>trap.tab");

open (OUT,">>general.tab");

print OUT "TOT UNITRAP ", scalar keys %seen_unitrap, "\n";
print OUT "TOT TRAP ", scalar keys %seen_trap, "\n";


my %unitrap_trap_freq;
print TR "Traps in unitrap\n";
foreach my $unitrap (keys %{$stat{'unitrap'}}){
	my $numbtraps = scalar @{$stat{'unitrap'}{$unitrap}{'trap'}};
	$unitrap_trap_freq{$numbtraps} ++;
	print TR "$numbtraps\t$unitrap\n";
}

print OUT "TOT TRAPS / UNITRAPS PER PROJECT\n";
foreach my $project (keys %{$stat{'project'}}){
	my $numbtraps = scalar @{$stat{'project'}{$project}{'trap'}};
	my $numbunitraps = scalar @{$stat{'project'}{$project}{'unitrap'}};
	print OUT "$project\t$numbtraps\t$numbunitraps\n";
	open (PR, ">>$project.tab");
	print PR "Traps in unitrap per project\n";
	foreach my $unitrap (keys %{$stat{'project'}{$project}{'trap_unitrap'}}){
		my $numbtraps = scalar @{$stat{'project'}{$project}{'trap_unitrap'}{$unitrap}};
		$unitrap_trap_freq{$numbtraps} ++;
		print PR"$project\t$numbtraps\t$unitrap\n";
	}
	
}


