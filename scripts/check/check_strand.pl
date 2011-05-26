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
use Bio::SearchIO;
use Bio::Unitrap::Fetch;
use Bio::Unitrap::Utils::File;

my $fetch = Bio::Unitrap::Fetch->new;
my %conf =  %::conf;
my $swarm = $conf{'global'}{'swarm'};
my $debug = $conf{'global'}{'debug'};
open (OUT,">>to_load.psl");

#total

my ($id, $idcount, $file,$help);
my $opt = &GetOptions( 
	    "file|f=s" =>    	\$file,
	    "id=s" => \$id,
	    "idcount=s" => \$idcount,
	    "help|h" => \$help
	    );


my $USAGE = "./check_strand.pl [-id] [-idcount]  [-file] || [-h help] ";

if ($help){print STDERR "$USAGE
                  -id = the starting trap db id
         		  -idcount = the amount of trap to annotate
         		  -f file : the alignment complete path
                  -h : this help\n";
			exit;
}

my %new;
my $searchio;
eval{$searchio = Bio::SearchIO->new (-file => $file, -format => 'psl') || die $!;};	
if ($@) {
	$debug && print STDERR "Could not create Bio::SearchIO object: $@\n";
}

while (my $result = $searchio->next_result) {
	my $name = $result->query_name;
	my @f = split /|/, $name;
	my $trapname;
	if ($f[0] eq 'ti'){$id = $f[1]; $trapname = $f[7]}
	else{$trapname = $name}
	my $res = `grep "|$trapname|" /nfs/users/nfs_v/vm6/Data/Trap/traps.psl`;
	foreach my $row (split /\n/,$res){
		my @fields = split /\t/,$row;
		my $start = $fields[15];
		print "<$trapname> <$start>\n";
		$new{$trapname}{$start} ++;
	}
}
				

my ($trapmap_id, $lastid);
unless($swarm){
	$trapmap_id = $id;
	$lastid = $trapmap_id + $idcount ;
}
else{
	$trapmap_id = ($id * $idcount) - $idcount + 1;
	$lastid = $id * $idcount;
}

for(my $idx = $trapmap_id; $idx <= $lastid; $idx ++){
	my $query = qq{select t.trap_id, t.trap_name, tm.start, tm.strand, tmr.region_id, r.region_strand from region r, trap t, trapmap tm, trapmap_region tmr where r.region_id = tmr.region_id and t.trap_id = tm.trap_id and tm.trapmap_id = tmr.trapmap_id };
	my $ref = $fetch->select_from_table($query);
	my $trap_id = $ref->{'trap_id'};
	my $trapmap_id = $idx;
	my $region_id = $ref->{'region_id'};
	my $tm_start = $ref->{'start'};
	my $tm_strand = $ref->{'strand'};
	my $r_strand = $ref->{'region_strand'};
	my $trap = $ref->{'trap_name'};
	my $res = `grep "|$trap|" /nfs/users/nfs_v/vm6/Data/Trap/traps.psl`;
	foreach my $row (split /\n/,$res){
		my @fields = split /\t/,$row;
		my $start = $fields[15];
		if($new{$trap}{$start} && $start == ($tm_start + 1)){$tm_start = $start}
		my $strand = $fields[8];
		if ($strand eq '-'){$strand = '-1'}
		if ($strand eq '+'){$strand = '1'}
		print "<$strand> <$tm_strand> <$start> <$tm_start>\n";
		if ($strand eq $tm_strand && $start == $tm_start){next;}
		my $amb = 0;
		if ($strand ne $r_strand){$amb = 1}
		my $update = qq{update trapmap tm, trapmap_region tmr set tm.strand = "$strand" and tm.start = $start and tmr.annotation_ambiguous = $amb where tm.trapmap_id = tmr.trapmap_id and tm.trapmap_id = $trapmap_id and tm.trap_id = $trap_id and tmr.region_id = $region_id and tm.start = $start};
		print STDERR "$update\n";
		$fetch->update($update);
	}
}	
	