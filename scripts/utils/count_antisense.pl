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
my %sense;
my %antisense;

open (OUT,">>to_load.psl");

#total


my $count_antisense_trap = qq{select tm.trap_id, t.trap_name, tmr.region_id, r.external_name, r.seq_id from insertion i, trap t, region r,trapmap tm, trapmap_region tmr where i.trap_id = tm.trap_id and i.trapmap_id = tm.trapmap_id and i.gene_id = tmr.region_id and tm.trap_id = t.trap_id and tmr.region_id = r.region_id and tm.trapmap_id = tmr.trapmap_id and tmr.annotation_ambiguous = 1};
my $count_sense_trap = qq{select tm.trap_id, t.trap_name, tmr.region_id, r.external_name, r.seq_id from insertion i, trap t, region r,trapmap tm, trapmap_region tmr where i.trap_id = tm.trap_id and i.trapmap_id = tm.trapmap_id and i.gene_id = tmr.region_id and tm.trap_id = t.trap_id and tmr.region_id = r.region_id and tm.trapmap_id = tmr.trapmap_id and tmr.annotation_ambiguous = 0};

foreach my $astrap (@{$fetch->select_many_from_table($count_antisense_trap)}){
	$antisense{$astrap->{'trap_name'}}{$astrap->{'external_name'}} ++;
}

foreach my $strap (@{$fetch->select_many_from_table($count_sense_trap)}){
	$sense{$strap->{'trap_name'}}{$strap->{'external_name'}} ++;
}


my $count_antisense = 0;
my $count_sense = 0;
my $count_both = 0;

foreach my $trap (keys %sense){
	if ($antisense{$trap}){$count_both ++}
	else{$count_sense ++;}
	
}

foreach my $trap (keys %antisense){
	next if $sense{$trap};
	foreach my $gene(keys %{$antisense{$trap}}){
		print STDERR "Check: $trap on $gene\n";	
	}
	$count_antisense ++;
}

print STDOUT "TOT SENSE FROM QUERY ", scalar keys %sense,"\n";
print STDOUT"TOT ANTISENSE FROM QUERY ", scalar keys %antisense,"\n";

print STDOUT"TOT ONLY SENSE $count_sense\n";
print STDOUT"TOT BOTH $count_both\n";
print STDOUT"TOT SENSE ", ($count_sense + $count_both),"\n";
print STDOUT"TOT ANTISENSE $count_antisense\n";

printf STDOUT"CONTROL:\ninitial tot = %d\n final tot = %d\n", ((scalar keys %sense) + (scalar keys %antisense)), ($count_sense + $count_both + $count_antisense);