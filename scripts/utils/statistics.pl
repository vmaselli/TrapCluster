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

#total
my $count_trap = qq{select count(trap_id) as tot from trap};
my $count_trap_per_project = qq{select count(t.trap_id) as tot, p.project_name from trap t, project p, esclone e where t.esclone_id = e.esclone_id and e.project_id = p.project_id group by p.project_id};

my $count_aligned_trap = qq{select count(distinct(tm.trap_id)) as tot from trap t, trapmap tm where tm.trap_id = t.trap_id };
my $count_annotated_trap = qq{select count(distinct(tm.trap_id)) as tot from trap t, trapmap tm, trapmap_region tmr, insertion i where tm.trapmap_id = tmr.trapmap_id and tm.trap_id = t.trap_id and i.trap_id = t.trap_id and i.trapmap_id = tm.trapmap_id and i.gene_id = tmr.region_id };
my $count_antisense_trap = qq{select count(distinct(tm.trap_id)) as tot from trap t, trapmap tm, trapmap_region tmr, insertion i where tm.trapmap_id = tmr.trapmap_id and tm.trap_id = t.trap_id and i.trap_id = t.trap_id and i.trapmap_id = tm.trapmap_id and i.gene_id = tmr.region_id and tmr.annotation_ambiguous = 1};
my $count_sense_trap = qq{select count(distinct(tm.trap_id)) as tot from trap t, trapmap tm, trapmap_region tmr, insertion i where tm.trapmap_id = tmr.trapmap_id and tm.trap_id = t.trap_id and i.trap_id = t.trap_id and i.trapmap_id = tm.trapmap_id and i.gene_id = tmr.region_id and tmr.annotation_ambiguous = 0};

#x project

my $count_aligned_trap_per_project = qq{select count(distinct(tm.trap_id)) as tot, p.project_name from trap t, project p, esclone e, trapmap tm where tm.trap_id = t.trap_id and t.esclone_id = e.esclone_id and e.project_id = p.project_id group by p.project_id};

#total annotated
my $count_annotated_trap_per_project = qq{select count(distinct(tm.trap_id)) as tot, p.project_name from trap t, project p, esclone e, trapmap tm, trapmap_region tmr, insertion i where tm.trapmap_id = tmr.trapmap_id and tm.trap_id = t.trap_id and t.esclone_id = e.esclone_id and e.project_id = p.project_id and i.trap_id = t.trap_id and i.trapmap_id = tm.trapmap_id and i.gene_id = tmr.region_id group by p.project_id};

#sense / antisense
my $count_sense_antisense_trap_per_project = qq{select count(distinct(tm.trap_id)) as tot, p.project_name, tmr.annotation_ambiguous as amb from trap t, project p, esclone e, trapmap tm, trapmap_region tmr, insertion i where tm.trapmap_id = tmr.trapmap_id and tm.trap_id = t.trap_id and t.esclone_id = e.esclone_id and e.project_id = p.project_id and i.trap_id = t.trap_id and i.trapmap_id = tm.trapmap_id and i.gene_id = tmr.region_id group by p.project_id, tmr.annotation_ambiguous};


#misc
my $gene_calls = qq{select count(distinct(external_name)) as tot from region};
my $sense_gene_calls = qq{select count(distinct(r.external_name)) as tot from region r, trapmap_region tmr where tmr.region_id = r.region_id and tmr.annotation_ambiguous = 0};
my $count_unitrap = qq{select count(unitrap_id) as tot from unitrap};

$stat{'total'}{'misc'}{'gene_calls'}{'total'} = $fetch->select_from_table($gene_calls)->{'tot'};
$stat{'total'}{'misc'}{'gene_calls'}{'sense'} = $fetch->select_from_table($sense_gene_calls)->{'tot'};
$stat{'total'}{'misc'}{'unitrap'}{'total'} = $fetch->select_from_table($count_unitrap)->{'tot'};

$stat{'total'}{'TotalTraps'} = $fetch->select_from_table($count_trap)->{'tot'};
foreach my $res (@{$fetch->select_many_from_table($count_trap_per_project)}){
	$stat{'total'}{'project'}{$res->{'project_name'}} = $res->{'tot'};
}

$stat{'total'}{'TotalAligned'} = $fetch->select_from_table($count_aligned_trap)->{'tot'};
foreach my $res (@{$fetch->select_many_from_table($count_aligned_trap_per_project)}){
	$stat{'aligned'}{$res->{'project_name'}} = $res->{'tot'};
}

$stat{'total'}{'TotalAnnotated'} = $fetch->select_from_table($count_annotated_trap)->{'tot'};
foreach my $res (@{$fetch->select_many_from_table($count_annotated_trap_per_project)}){
	$stat{'annotated'}{'total'}{$res->{'project_name'}} = $res->{'tot'};
}

$stat{'total'}{'TotalAntisense'} = $fetch->select_from_table($count_antisense_trap)->{'tot'};
$stat{'total'}{'TotalSense'} = $fetch->select_from_table($count_sense_trap)->{'tot'};
foreach my $res (@{$fetch->select_many_from_table($count_sense_antisense_trap_per_project)}){
	my $amb;
	$amb = "antisense" if $res->{'amb'} == 1;
	$amb = "sense" if $res->{'amb'} == 0;
	$stat{'annotated'}{$amb}{$res->{'project_name'}}= $res->{'tot'};
}


print "General\n";
print "Total Tags\t$stat{'total'}{'TotalTraps'}\n";
printf "Aligned Tags\t$stat{'total'}{'TotalAligned'}\t%.2f%s  of Total Tags\n", (($stat{'total'}{'TotalAligned'}/$stat{'total'}{'TotalTraps'})*100),"%";
printf "Annotated Tags\t$stat{'total'}{'TotalAnnotated'}\t%.2f%s  of Total Aligned\n", (($stat{'total'}{'TotalAnnotated'}/$stat{'total'}{'TotalAligned'})*100),"%";
printf "Antisense\t$stat{'total'}{'TotalAntisense'}\t%.2f%s of Total Annotated\n", (($stat{'total'}{'TotalAntisense'}/$stat{'total'}{'TotalAnnotated'})*100),"%";
printf "Sense\t$stat{'total'}{'TotalSense'}\t%.2f%s  of Total Annotated\n", (($stat{'total'}{'TotalSense'}/$stat{'total'}{'TotalAnnotated'})*100),"%";
print "Total Gene Calls\t$stat{'total'}{'misc'}{'gene_calls'}{'total'}\n";
printf "Sense Gene Calls\t$stat{'total'}{'misc'}{'gene_calls'}{'sense'}\t%.2f%s  of Total Annotated\n",(($stat{'total'}{'misc'}{'gene_calls'}{'sense'}/$stat{'total'}{'misc'}{'gene_calls'}{'total'})*100),"%";
print "No Unitraps\t$stat{'total'}{'misc'}{'unitrap'}{'total'}\n";
print "\n";

print "Total\n";
foreach my $key (sort keys %{$stat{'total'}{'project'}}){
	printf "$key\t$stat{'total'}{'project'}{$key}\t%.2f%s\t\t\n", (($stat{'total'}{'project'}{$key}/$stat{'total'}{'TotalAligned'})*100),"%";
}
print "\nAligned\n";
foreach my $key (sort keys %{$stat{'aligned'}}){
	printf "$key\t$stat{'aligned'}{$key}\t%.2f%s\t\t\n", (($stat{'aligned'}{$key}/$stat{'total'}{'project'}{$key})*100),"%";
}
print "\nAnnotated\n";
foreach my $key (sort keys %{$stat{'annotated'}{'total'}}){
	printf "$key\t$stat{'annotated'}{'total'}{$key}\t%.2f%s\t\t\n", (($stat{'annotated'}{'total'}{$key}/$stat{'total'}{'project'}{$key})*100),"%";
}
print "\n";
print "Annotated Sense\n";
foreach my $key (sort keys %{$stat{'annotated'}{'sense'}}){
	printf "$key\t$stat{'annotated'}{'sense'}{$key}\t%.2f%s\t\t\n",(($stat{'annotated'}{'sense'}{$key}/$stat{'annotated'}{'total'}{$key})*100),"%";
}
print "\n";
print "Annotated AntiSense\n";
foreach my $key (sort keys %{$stat{'annotated'}{'antisense'}}){
	printf "$key\t$stat{'annotated'}{'antisense'}{$key}\t%.2f%s\t\t\n", (($stat{'annotated'}{'antisense'}{$key}/$stat{'annotated'}{'total'}{$key})*100),"%";
}

