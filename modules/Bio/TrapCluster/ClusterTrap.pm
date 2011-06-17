#!/usr/bin/perl -w

=head1 Bio::TrapCluster::ClusterTrap

=head2 Authors

=head3 Created by

             Vincenza Maselli
             v.maselli@ucl.ac.uk

=head3 Original script by

              Guglielmo Roma
              guglielmoroma@libero.it

=head2 Description
        
             This module run blast or read information from blastoutput or from database and load the tables `trapblock_Cluster` and `region`	
             
=head2 Usage

	    my $obj = Bio::TrapCluster::ClusterTrap->new;
            
=cut

package Bio::TrapCluster::ClusterTrap;

use strict;
use DBI;
use Carp;
use Data::Dumper;
use vars qw(@ISA);
use File::Spec;
use Bio::TrapCluster::Utils::Argument qw(rearrange);
use Bio::TrapCluster::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

#setting global variables

require "$ENV{'TrapCluster'}/trapcluster_conf.pl";

my %conf =  %::conf;
my $debug = $conf{'global'}{'debug'};

use Bio::SeqIO;
use Bio::SearchIO;
use Bio::MCE::Range;
use Bio::TrapCluster::AnnotationTrap;
@ISA = qw(Bio::Root::Root Bio::TrapCluster::AnnotationTrap);

=head2 new

  Arg [..]: Take a set of named argumnts from a config file
  Example: my $retrieve = Bio::TrapCluster::RetrieveTrap->new
  Description:
  Returntype:
  Exceptions: source (fromfile, fromdb, fromquery) not defined;
  Caller:
  Status: 

=cut

sub new{
  	my $caller = shift;

 	my $class = ref($caller) || $caller;
  	my $self = $class->SUPER::new(@_);
	
  	my (
   		$trapmap_arrayref
    )
    = $self->_rearrange( [
      	'TRAPMAP_ARRAYREF'
    ],
    @_
    );
	
	$trapmap_arrayref && $self->trapmap_arrayref($trapmap_arrayref);
	
  	
	
	my $annotation = Bio::TrapCluster::AnnotationTrap->new;
  	$self->annotation($annotation);
  	
  	return $self;
}

# sub subname{
# 
#   my ($self, $value) = @_;
#   $self->{'subname'} = $value if defined $value;
#   
#   return $self->{'subname'};
# }

sub trapmap_arrayref{

  my ($self, $value) = @_;
  $self->{'trapmap_arrayref'} = $value if defined $value;
  
  return $self->{'trapmap_arrayref'};
}


sub load{
  my ($self) = @_;
  return $self->annotation->load;
}

sub annotation{
  my ($self, $value) = @_;
  $self->{'annotation'} = $value if defined $value;
  
  return $self->{'annotation'};
}


sub disc_ranges () {
	my ($self,$str) = @_;

	my @ranges;
	
	my $min_start = '1000000000000000000';
	my $max_end = '0';
	
	foreach my $coords (@{$self->trapmap_arrayref}) {
		my $start = $coords->{'start'};
		my $end = $coords->{'end'};
		
		if ($start < $min_start) {
			$min_start = $start;
		}
		
		if ($end > $max_end) {
			$max_end = $end;
		}
		
		#$debug && print STDOUT "Trying to range: START $start - END $end\n";
		
		if ($start > $end) {
			print STDERR "ERROR $start - $end... iinverting the coords!!!\n";
			exit;
		}
		
		if ($start && $end) {
			my $range = new Bio::Range (-start=>$start, -end=>$end, -strand=>$str);
			push @ranges, $range;
		} else {
			print STDERR "\n########\nError. There isn't start or end, or both\n\n";
		}
	}
	
	$debug && print STDOUT "Disconnecting ranges\n";
	my @sorted_ranges = sort {$a->start <=> $b->start} @ranges;
	my $disc_ranges = Bio::MCE::Range->disconnected_ranges($self->load->fetch->db_connection,\@sorted_ranges,1);
	if ($disc_ranges) {
		return $disc_ranges, $min_start, $max_end;
	} else {
		return undef;
	}
}



sub get_next_accession () {
    my ($self,$table,$prefix) = @_;
	my $accession;
	my $table_id = $table."_id";
	my $sql_acc = qq{select accession from $table order by $table_id desc limit 1};
	my $acc = $self->load->fetch->select_from_table($sql_acc);
	
	
	my $last = $acc->{'accession'};
	
	if (!$last) {
		$accession = $prefix."1";
	} else {
	    my @numbers = split (/[A-Za-z]+/, $last);
		
		my $number = $numbers[1];
		my $next_number = $number + 1;
		
		$accession = $prefix.$next_number;
	}
	
	$debug && print STDOUT "LAST $table accession: $last \n" if $last;
	$debug && print STDOUT "NEW $table accession: $accession\n";
	
	return $accession;
}

sub get_seq_from_id () {

	my ($self,$table, $field,$id, $strand) = @_;
	my ($sql, $cluster_seq);
	
	if ($strand == 1) {
		$sql = "select sequence from $table where $field = $id order by start;";
	} elsif ($strand == -1) {
		$sql = "select sequence from $table where $field = $id order by start desc;";
	}
	
	;
	
	my $i = 0;
	
	foreach my $seqs (@{$self->load->fetch->select_many_from_table($sql)}) {
		my $seq;
		
		if ($i%2 == 0) {
			$seq = lc($seqs->{'sequence'});
		} else {
			$seq = uc($seqs->{'sequence'});
		}
		
		$cluster_seq .= $seq;
		$i++;
	}
	
	return $cluster_seq;
}

sub run{
	my ($self, $chr, $version, $region,$str) = @_;
	$chr =~ s/chr//;
	my ($disc_ranges, $min_start, $max_end) = $self->disc_ranges;
	my $fetch = $self->load->fetch;
	my $hit_db = $conf{'default'}{'hit_db'};
	my $maxicluster_accession = $self->get_next_accession ('maxicluster','MCL');
	my (%maxiclustermap_toinsert,%maxicluster_toinsert);
	$maxicluster_toinsert{'accession'} = $maxicluster_accession;
	$maxicluster_toinsert{'mol_type'} = 'mRNA';
	
	my $maxicluster_id = $self->load->load_maxicluster(\%maxicluster_toinsert);
	foreach my $r (@{$disc_ranges}) {
		$debug && print "Got here 237\n";
		
		$maxiclustermap_toinsert{'maxicluster_id'} = $maxicluster_id;
		$maxiclustermap_toinsert{'hit_id'} = $chr;
		$maxiclustermap_toinsert{'hit_db'} = $version;
		$maxiclustermap_toinsert{'strand'} = $str;
		$maxiclustermap_toinsert{'start'} = $r->start;
		$maxiclustermap_toinsert{'end'} = $r->end;
		
		print "MAX ",$r->start," ",$r->end,"\n";
		
		my $maxiclustermap_id = $self->load->load_maxiclustermap(\%maxiclustermap_toinsert);
		
		### Retrieve all the trapmaps overlapping the current maxicluster
		my $sql_map = "select distinct tm.trap_id, tm.trapmap_id from trap t, trapmap tm where tm.start <= '".$r->end."' and tm.end >= '".$r->start."' order by tm.start,tm.end;";
		
		my (@ranges,%trap_maxicluster_info);
	
	
		#$debug && print STDERR "\t\t\tline 259\n$sql_map";	
		### Foreach trapmap located within the current maxicluster,
		### push the trapblocks in an array, named ranges
		foreach my $maps (@{$self->load->fetch->select_many_from_table($sql_map)}) {
			my $trap_id = $maps->{'trap_id'};
			my $trapmap_id = $maps->{'trapmap_id'};
			#$debug && print STDERR "\t\tMAP $trap_id $trapmap_id\n";
			$trap_maxicluster_info{$trap_id}{'trapmap'} = $trapmap_id;
			
			### Retrieve all the trapblocks of the current trapmap
			my $sql_block = "select distinct tb.start, tb.end from trapblock tb where tb.trapmap_id = '$trapmap_id' order by tb.start,tb.end;";
			
			foreach my $tbs (@{$self->load->fetch->select_many_from_table($sql_block)}) {
				my $tbs_start = $tbs->{'start'};
				my $tbs_end = $tbs->{'end'};
				#$debug && print STDERR "\t\tBLOCK $tbs_start $tbs_end\n";
				### Push the trapblocks in the array
				if ($tbs_start && $tbs_end) {
					my $range = new Bio::Range (-start=>$tbs_start, -end=>$tbs_end, -strand=>$str);
					push @ranges, $range;
				} else {
					print STDERR "\n########\nError. There isn't start or end, or both\n\n";
				}
			}
		}
		$debug && print STDERR "\t\t\tline 284\n";
		### Make a disconnected ranges of all the trapblocks present within the current maxicluster
		### in order to create maxiclusterblocks
		my @sorted_ranges = sort {$a->start <=> $b->start} @ranges;
	my @tb_disc_ranges = @{Bio::MCE::Range->disconnected_ranges($self->load->fetch->db_connection,\@sorted_ranges,1)};
		
		### These are maxiclusterblocks
		foreach my $r (@tb_disc_ranges) {
			### Filling the table maxiclusterblock
			my %toinsert_maxiclusterblocks;
			$toinsert_maxiclusterblocks{'maxiclustermap_id'} = $maxiclustermap_id;
			$toinsert_maxiclusterblocks{'start'} = $r->start;
			$toinsert_maxiclusterblocks{'end'} = $r->end;
			$toinsert_maxiclusterblocks{'strand'} = $str;					
			my $slice = $self->annotation->slicecoreadpt->fetch_by_region ($region,$chr,$r->start,$r->end,$str);
			
			$toinsert_maxiclusterblocks{'sequence'} = $slice->seq();
			
			my $maxiclusterblocks_id = $self->load->load_maxiclusterblock(\%toinsert_maxiclusterblocks);
		}
		
		### Filling the table trap_maxicluster
		my %toinsert_trap_maxicluster;
		my @inserted_traps;
		### Retrieving the full sequence of the maxicluster simply joining the sequence of each maxiclusterblock
		my $maxicluster_seq = $self->get_seq_from_id ('maxiclusterblock',  'maxiclustermap_id',$maxiclustermap_id, $str);
		my $sth_seq = qq{update maxicluster set sequence = \"$maxicluster_seq\" where maxicluster_id = $maxicluster_id};
		
		$self->load->fetch->update($sth_seq);
		
		foreach my $trap (keys %trap_maxicluster_info) {
			push @inserted_traps, $trap;
			my $trapmap = $trap_maxicluster_info{$trap}{'trapmap'};
			
			$toinsert_trap_maxicluster{'maxicluster_id'} = $maxicluster_id;
			$toinsert_trap_maxicluster{'trap_id'} = $trap;
			$toinsert_trap_maxicluster{'trapmap_id'} = $trapmap;
			
			my $trap_maxicluster_id = $self->load->load_trap_maxicluster(\%toinsert_trap_maxicluster);
		}
		$debug && print STDERR "\t== MAXICLUSTER DONE == CREATE TRAPCLUSTER STARTING\n";
		my $feat_hash = $self->create_trapcluster($region,$maxicluster_id, $maxiclustermap_id, $hit_db, $chr, $str, $debug);
		$debug && print "\n\n\t== TRAPCLUSTER CREATED START WITH ANNOTATION ==\n\n";
		$self->annotate($region,$chr,$feat_hash) if $conf{'annotation'}{'do'};
		
	}
}

sub create_trapcluster{                                                                                                                                                                                                         
	my ($self,$region,$maxicluster_id, $maxiclustermap_id, $hit_db, $chr, $strand, $debug) = @_;
	my %hash;
	my $sql = "select maxiclusterblock_id from maxiclusterblock where maxiclustermap_id = '$maxiclustermap_id' order by start,end;";
	
	$debug && print STDERR "create_trapcluster line 334\n";
	
	foreach my $blocks (@{$self->load->fetch->select_many_from_table($sql)}) {
		my $maxiclusterblock_id = $blocks->{'maxiclusterblock_id'};
		my $sql_block = "select start, end from maxiclusterblock where maxiclusterblock_id = '$maxiclusterblock_id' and checked = '0'";
		#$debug && print STDERR "$sql_block\n";
		
		my $block_info = $self->load->fetch->select_from_table($sql_block);
		my $start = $block_info->{'start'};
		my $end = $block_info->{'end'};
		
		if ($start && $end) {
			my %blocks;
			$blocks{$maxiclusterblock_id}{'start'} = $start;
			$blocks{$maxiclusterblock_id}{'end'} = $end;
			
			my ($features, $trapcluster_id) = $self->find_overlapping_traps($region,$maxicluster_id, $maxiclustermap_id, \%blocks, $hit_db, $chr, $strand, $debug);
			$hash{$trapcluster_id} = $features;
		}
	}
	return \%hash;
}

sub find_overlapping_traps{
	
	my ($self,$region,$maxicluster_id, $maxiclustermap_id, $blocks, $hit_db, $chr, $strand, $debug) = @_;
	my %blocks = %$blocks;
	my ($where,$where2);
	my $i = 0;
	
	my $min_start = 100000000000;
	my $max_end = 0;
	
	$debug && print STDERR "find_overlapping_traps line 366\n";
	
	foreach my $block_id (keys %blocks) {
		if ($i == 0) {
			$where .= " ( ";
		} else {
			$where .= " or ";
		}
		
		my $bstart = $blocks->{$block_id}->{'start'};
		my $bend = $blocks->{$block_id}->{'end'};
		
		$where .= " (tb.start <= '$bend' and tb.end >= '$bstart') ";
		
		if ($bstart < $min_start) {
			$min_start = $bstart;
		}
		
		if ($bend > $max_end) {
			$max_end = $bend;
		}
		
		$i++;
	}
	
	$where .= " ) ";
	
	$debug && print STDERR "find_overlapping_traps line 393\n";
	
	#### retrieving all the blocks of the traps overlapping the input coordinates
	my $sql = "select distinct tx.trap_id, tb.trapmap_id from trap_maxicluster tx, trapblock tb where tb.trapmap_id = tx.trapmap_id and tx.maxicluster_id = '$maxicluster_id' and $where;";
	
	my %trap_trapcluster_info;
	my $j = 0;
	
	foreach my $traps (@{$self->load->fetch->select_many_from_table($sql)}) {
		my $trap_id = $traps->{'trap_id'};
		my $trapmap_id = $traps->{'trapmap_id'};
		$trap_trapcluster_info{$trap_id}{'trapmap'} = $trapmap_id;
		
		my $sql_tm = "select distinct tb.start, tb.end from trapblock tb where trapmap_id = '$trapmap_id';";
		foreach my $trap_blocks (@{$self->load->fetch->select_many_from_table($sql_tm)}) {		
			my $trap_start = $trap_blocks->{'start'};
			my $trap_end = $trap_blocks->{'end'};
			
			if ($j == 0) {
				$where2 .= " ( ";
			} else {
				$where2 .= " or ";
			}
			
			$where2 .= " (mb.start <= '$trap_end' and mb.end >= '$trap_start') ";
			
			$j++;
		}
	}
	
	$where2 .= " ) ";
	
	$debug && print STDERR "find_overlapping_traps line 425\n";
	#print "######### J $j == O $overlapping_trapblocks\n";
	#if ($j == $overlapping_trapblocks) {
	### filling trapcluster
	my %trapcluster_toinsert;
	
	my $trapcluster_accession = $self->get_next_accession ('trapcluster', 'TCL');
	$trapcluster_toinsert{'accession'} = $trapcluster_accession;
	$trapcluster_toinsert{'maxicluster_id'} = $maxicluster_id;
	
	### storing links to ensembl and ucsc genome browsers (adding a padding of 1000)
	my $link_start = $min_start - 1000;
	my $link_end = $max_end + 1000;
	$trapcluster_toinsert{'link_to_ensembl'} = $conf{'default'}{'ensembl_base_url'}."/contigview?chr=$chr&&vc_start=$link_start&vc_end=$link_end&x=0&y=0";
	$trapcluster_toinsert{'link_to_ucsc'} = "http://genome.ucsc.edu/cgi-bin/hgTracks?org=Mouse&position=chr$chr:$link_start-$link_end&db=".$conf{'default'}{'db_name'};
	
	my $trapcluster_id = $self->load->load_trapcluster(\%trapcluster_toinsert);
	my $trapclustermap_id;
	
	if ($trapcluster_id) {
		### filling trapclustermap
		my %trapclustermap_toinsert;
		
		$trapclustermap_toinsert{'trapcluster_id'} = $trapcluster_id;
		$trapclustermap_toinsert{'hit_id'} = $chr;
		$trapclustermap_toinsert{'hit_db'} = $hit_db;
		$trapclustermap_toinsert{'strand'} = $strand;
		$trapclustermap_toinsert{'start'} = $min_start;
		$trapclustermap_toinsert{'end'} = $max_end;
		
		$trapclustermap_id = $self->load->load_trapclustermap(\%trapclustermap_toinsert);	
	}
	
	my @features;
	foreach my $block_id (keys %blocks) {
		my $sql_tclblock = "select * from maxiclusterblock where maxiclusterblock_id = '$block_id';";
		my $tclblocks = $self->load->fetch->select_from_table($sql_tclblock);
		
		### filling trapclusterblocks
		my %toinsert_trapclusterblocks;
		$toinsert_trapclusterblocks{'trapclustermap_id'} = $trapclustermap_id;
		$toinsert_trapclusterblocks{'start'} = $tclblocks->{'start'};
		$toinsert_trapclusterblocks{'end'} = $tclblocks->{'end'};
		$toinsert_trapclusterblocks{'strand'} = $strand;
		$toinsert_trapclusterblocks{'sequence'} = $tclblocks->{'sequence'};
		
		my $trapclusterblock_id = $self->load->load_trapclusterblock(\%toinsert_trapclusterblocks);
		
		### set checked for blocks in the maxiclusterblock table!!!
		my $sql_block = "update maxiclusterblock set checked = '1' where maxiclusterblock_id = '$block_id';";
		$self->load->fetch->update($sql_block);
		
		### Create a feature for each trapclusterblock
		my $feature = Bio::SeqFeature::Generic->new(	-display_name => $trapclusterblock_id,
								-start => $tclblocks->{'start'},
								-end => $tclblocks->{'end'},
								-strand => $strand	);
		push @features,$feature;
	}
	
	### filling trap_trapcluster
	my %toinsert_trap_trapcluster;
	foreach my $trap (keys %trap_trapcluster_info) {
		my $trapmap = $trap_trapcluster_info{$trap}{'trapmap'};
		
		$toinsert_trap_trapcluster{'trapcluster_id'} = $trapcluster_id;
		$toinsert_trap_trapcluster{'trap_id'} = $trap;
		$toinsert_trap_trapcluster{'trapmap_id'} = $trapmap;
		
		my $trap_trapcluster_id = $self->load->load_trap_trapcluster(\%toinsert_trap_trapcluster);
	}
	
	### Once the sequence of each tclblock is available, we update the sequence of the whole TCL, as well as the number of tclblocks and the n of traps
	my $trapcluster_seq = $self->get_seq_from_id ('trapclusterblock',  'trapclustermap_id',$trapclustermap_id, $strand, $debug);
	my $sth_seq = "update trapcluster set sequence = '$trapcluster_seq', trapclusterblocks = '".scalar(keys %blocks)."', traps ='".scalar(keys %trap_trapcluster_info)."' where trapcluster_id = '$trapcluster_id';";
	
	$self->load->fetch->update($sth_seq);
	
	return (\@features, $trapcluster_id);
	#}


}

sub annotate{
	my ($self,$region,$chr,$hash) = @_;
	
	foreach my $trapcluster_id (keys %{$hash}){
		my $feats = $hash->{$trapcluster_id};
		my $trapcluster = $self->load->fetch->get_trapcluster_by_id($trapcluster_id);
		
		foreach my $feat (@{$feats}){
			;
			my $slice_core_adaptor = $self->annotation->slicecoreadpt;
			my $slice_core = $slice_core_adaptor->fetch_by_region($region,$chr, $feat->{'start'}, $feat->{'end'});
			my $found = 0;
			print STDERR "---- FOUND $found try to do ENSGENE ".$self->annotation->do_ensgene."\n";
			if ($self->annotation->do_ensgene) {
				
				if (defined $slice_core){
					$found = $self->annotate_with_ensembl_gene($trapcluster,$feat, $slice_core);
				}
			}
			if ($found) {
				my $timestmp = `date`;
				print STDERR "---- DONE annotate gene with annotate_with_ensembl_gene\n";
				#return $found;
			}
			print STDERR "---- FOUND $found try to do EST ".$self->annotation->do_ensestgene."\n";
			if ($self->annotation->do_ensestgene){
				my $slice_est_adaptor = $self->annotation->sliceestadpt;
				my $slice_est = $slice_est_adaptor->fetch_by_region($region,$chr, $feat->{'start'}, $feat->{'end'});
				print STDERR "---- got here 838\n";
				if (defined $slice_est){
					print STDERR "---- SLICE EST defined\n";
					$found = $self->annotate_with_ensembl_estgene($trapcluster,$feat, $slice_est);
				}
				else{print STDERR "---- SLICE EST not defined\n"; }
				print STDERR "---- got here 844\n";
			}
			else{print STDERR "---- got here 845\n";}
			if ($found) {
				my $timestmp = `date`;
				print STDERR "---- DONE annotate ensembl EST gene with annotate_with_ensembl_estgene\n";
				#return $found;
			}
			
			print STDERR "---- FOUND $found try to do GENESCAN ".$self->annotation->do_genescan."\n";
	
			if ($self->annotation->do_genescan) {
				print STDERR "annotate_with_ensembl_genescan\n";
				$found = $self->annotate_with_ensembl_genescan($trapcluster,$feat, $slice_core);
			}
			
			print STDERR "---- FOUND $found try to do UNIGENE ".$self->annotation->do_unigene."\n";
			if ($self->annotation->do_unigene) {
				print STDERR "annotate_with_unigene\n";
				$found = $self->annotate_with_unigene($trapcluster,$feat, $slice_core);
				
			}
			print STDERR "---- FOUND $found try to do CDNA ".$self->annotation->do_mouse_cdna."\n";
			if ($self->annotation->do_mouse_cdna) {
				print STDERR "annotate_with_cdna\n";
				$found = $self->annotate_with_mouse_cdna($trapcluster,$feat, $slice_core);
			}
			print STDERR "---- FOUND $found try to do ENSEST ".$self->annotation->do_ensest."\n";
			if ($self->annotation->do_ensest) {
				print STDERR "annotate_ensest\n";
				$found = $self->annotate_with_ensest($trapcluster,$feat, $slice_core);
			}
			print STDERR "---- FOUND $found try to do REPEATS ".$self->annotation->do_ensrepeats."\n";
			if ($self->annotation->do_ensrepeats) {
				print STDERR "annotate_repeats\n";
				$found = $self->annotate_with_ensrepeats($trapcluster,$feat, $slice_core);
			}
			print STDERR "---- FOUND $found. END\n";

		}
	}
}


sub annotate_with_ensembl_genescan {
	my ($self, $trapcluster,$feat, $slice) = @_;
	my $trapcluster_id = $trapcluster->{'trapcluster_id'};
	my $type = "genescan prediction";
	my $moltype = $trapcluster->{'mol_type'};
	my $seqdir = $trapcluster->{'sequencing_direction'};
	my $accession = $trapcluster->{'accession'};
	my $hash;
	my $yes_gene = 0;
	my $trapclustermap_id = $feat->display_name;
	my $trapblock = $self->load->fetch->get_trapclusterblock_by_trapclustermap_id($trapclustermap_id);
	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		return 0;
	}
	
	#WHICH IS CALLED GENE HERE IS ACTUALLY A PREDICT TRANSCRIPT
	$hash->{'trapclustermap_region'}{'trapclustermap_id'} = $trapclustermap_id;
	print STDERR "------ got here 706\n";
	foreach my $gene (@{$slice->get_all_PredictionTranscripts}) {
		my $gstart = $gene->seq_region_start;
		my $gend = $gene->seq_region_end;
		my $gstrand = $gene->strand;

	    #test if the feature (trapclustermap) inteserct a gene
		my ($ti) = $self->annotation->trap_utility->intersection($gene,$gstart,$gend,$feat);
		print STDERR "------ got here 714\n";
		if ($ti->{'i'}) {
			$yes_gene ++;
			
			my $gene_id = $self->annotation->region->gene($gene); #load gene in region table
			print STDERR "------ got here 719\n";
			my $amb = $self->annotation->trap_utility->check_strand($trapcluster,$feat,$gene); #used to check if the alignment is on the correct strand of the gene
			#if ($debug && $amb == 1){next;}
			
			#if the strand is not correct it isn't possible to infer the insertion
			$debug && print STDERR "FOUND ".$gene->biotype.": ".$gene->stable_id." and AMB = $amb\n";			

			#set generic info
			$hash->{'trapclustermap_region'}{'region_id'} = $gene_id;
			$hash->{'trapclustermap_region'}{'overlap'} = $ti->{'coverage'};
			$hash->{'trapclustermap_region'}{'type'} = $gene->biotype;
			$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = $amb;				
			$hash->{'trapclustermap_region'}{'number_trapblocks'} = scalar @{$trapblock};
			my $number_annotate_trapblocks = 0;

			my @regions = @{$self->annotation->region->rearrange_predict_exons($gene,$gene_id)};
			my $first = $regions[0]; 
			my $last = $regions[$#regions];
			my $first_exon_id = $first->{'region_id'};
			my $last_exon_id = $last->{'region_id'};
			my $first_start = $first->{'region_start'}; 
			my $last_end = $last->{'region_end'};
			$hash->{'trapclustermap_region'}{'total_number_exons'} = 1+(int(scalar @regions)/2);
			my ($nti) = $self->annotation->trap_utility->intersection($regions[0],$first_start,$last_end,$feat);			
			unless ($nti->{'i'}) {
				print STDERR "no coding block found for gene: $gene_id\n"; 
				$hash->{'trapclustermap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
				$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
				next;
			}

			my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);

			#load hash

			$hash->{'additional'}{'processing'} = "annotate with genescan prediction";
			$hash->{'additional'}{'comment'} = "";
			$hash->{'additional'}{'label'} = $self->annotation->region->get_label($gene);
			$hash->{'additional'}{'trapcluster_id'} = $trapcluster_id;
			$hash->{'additional'}{'user'} = "mysql_dev";
			$hash->{'additional'}{'note'} = 'GENESCAN';
			$self->load->load_trapclusteradditional($hash->{'additional'});

			#$hash->{'trapclustermap_region'}{'trapblock_id'} = $trapblock_id;
			$hash->{'trapclustermap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
			$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});

		}#end intersect gene
		print STDERR "------ got here 448\n";
	}#end gene loop

	if ($yes_gene) {
		my $timestmp = `date`;
		print STDERR "------ DONE annotate gene with annotate_with_ensembl_genescan\n";
		return 1;
	}else{
		print STDERR "------ no genescan prediction found\n";
		#print STDERR Dumper $trapcluster;
		my $region_id = $self->annotation->region->slice($slice);

		$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapclustermap_region'}{'type'} = "genomic";
		$hash->{'trapclustermap_region'}{'overlap'} = 100;
		$hash->{'trapclustermap_region'}{'region_id'} = $region_id;

		$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
		my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_ensembl_genescan\n";
		return 0;
	}
	
	
	
}

sub annotate_with_unigene {
	my ($self, $trapcluster,$feat, $slice) = @_;
	my $unigene_coverage = $conf{'annotation'}{'unigene_coverage'};
	my $unigene_perc_id = $conf{'annotation'}{'unigene_perc_id'};
	my $trapcluster_id = $trapcluster->{'trapcluster_id'};
	my $type = "unigene";
	my $moltype = $trapcluster->{'mol_type'};
	my $seqdir = $trapcluster->{'sequencing_direction'};
	my $accession = $trapcluster->{'accession'};
	my $hash;
	my $yes_gene = 0;
	my $trapclustermap_id = $feat->display_name;
	my $trapblock = $self->load->fetch->get_trapclusterblock_by_trapclustermap_id($trapclustermap_id);	
	
	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		return 0;
	}
	
	$hash->{'trapclustermap_region'}{'trapclustermap_id'} = $trapclustermap_id;
	
	my @dna_dna_align_feats = @{$slice->get_all_DnaAlignFeatures ('Unigene')};
	my @dnadnafeat;
	my $ori_start = $feat->start;
	my $ori_end = $feat->end;
	my $s = 1;
	my $e = ($feat->end - $feat->start) + 1;
	$feat->start($s);
	$feat->end($e);
	
	foreach my $dnadna (@dna_dna_align_feats) {
		if ($dnadna->percent_id >= $unigene_perc_id && $dnadna->hseqname =~ /Mm/) {
			
			my $feature = Bio::SeqFeature::Generic->new (	-display_name => $dnadna->hseqname,
									-start => $dnadna->hstart,
									-end => $dnadna->hend,
									-strand => $dnadna->hstrand);
			push @dnadnafeat ,$feature;
			
			my $col = new Bio::SeqFeature::Collection();
			my $totaladded = $col->add_features(\@dnadnafeat);
			my $amb = 0;
			if ($dnadna->hstrand ne $feat->strand){$amb = 1}
			$debug && print "FOUND ".$dnadna->hseqname."\n";
			my @subset = $col->features_in_range (-range=>$feat, -contain=> 0);	
			my $coverage = $unigene_coverage;
			foreach my $s (@subset) {	
				my $inter = $feat->intersection($s);	
				my $cov = ($inter->length/$s->length)*100;
				unless (defined $coverage){$coverage = $cov};
				if ($cov >= $coverage) {
					$debug && print "COVERAGE IS RIGHT\n"; 
					my $ghash;
					$ghash->{'seq_id'} = $s->display_name;
					my $start = $ori_start - $s->start;
					if ($feat->start < $s->start){$start = $ori_start + $s->start}
					my $end = $ori_start - $s->end;
					if ($feat->end < $s->end){$end = $ori_end + $s->end}
					$ghash->{'start'} = $start;
					$ghash->{'end'} = $end;
					$ghash->{'name'} = $feat->display_name;
					$ghash->{'strand'} = $dnadna->hstrand;
					$ghash->{'parent_id'} = 0;
					$ghash->{'rank'} = 0;
					$ghash->{'description'} = 'Unigene';
					
					my $gene_id = $self->load->load_region($ghash);
					
					$hash->{'trapclustermap_region'}{'region_id'} = $gene_id;
					$hash->{'trapclustermap_region'}{'overlap'} = $coverage;
					$hash->{'trapclustermap_region'}{'type'} = 'Unigene';
					$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = $amb;				
					$hash->{'trapclustermap_region'}{'number_trapblocks'} = scalar @{$trapblock};

					my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);
				
					#load hash
				
					$hash->{'additional'}{'processing'} = "annotate with Unigene DnaAlignFeatures";
					$hash->{'additional'}{'comment'} = "";
					$hash->{'additional'}{'label'} = "";
					$hash->{'additional'}{'trapcluster_id'} = $trapcluster_id;
					$hash->{'additional'}{'user'} = "mysql_dev";
					$hash->{'additional'}{'note'} = 'UNIGENE';
					$self->load->load_trapclusteradditional($hash->{'additional'});
				
					#$hash->{'trapclustermap_region'}{'trapblock_id'} = $trapblock_id;
					$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
					$yes_gene ++;
				}
				else{$debug && print "COVERAGE BELOW THRESHOLD\n";}
			}
		}
	}
	
	print "YES $yes_gene\n";	
	
	if ($yes_gene) {
		my $timestmp = `date`;
		print STDERR "------ DONE annotate gene with annotate_with_unigene\n";
		return 1;
	}else{
		print STDERR "------ no unigene found\n";
		#print STDERR Dumper $trapcluster;
		my $region_id = $self->annotation->region->slice($slice);

		$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapclustermap_region'}{'type'} = "genomic";
		$hash->{'trapclustermap_region'}{'overlap'} = 100;
		$hash->{'trapclustermap_region'}{'region_id'} = $region_id;

		$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
		my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_unigene\n";
		return 0;
	}
}

sub annotate_with_ensest {
	my ($self, $trapcluster,$feat, $slice) = @_;
	my $est_coverage = $conf{'annotation'}{'est_coverage'};
	my $est_perc_id = $conf{'annotation'}{'est_perc_id'};
	my $trapcluster_id = $trapcluster->{'trapcluster_id'};
	my $type = "est";
	my $moltype = $trapcluster->{'mol_type'};
	my $seqdir = $trapcluster->{'sequencing_direction'};
	my $accession = $trapcluster->{'accession'};
	my $hash;
	my $yes_gene = 0;
	my $trapclustermap_id = $feat->display_name;
	my $trapblock = $self->load->fetch->get_trapclusterblock_by_trapclustermap_id($trapclustermap_id);	
	
	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		return 0;
	}
	
	$hash->{'trapclustermap_region'}{'trapclustermap_id'} = $trapclustermap_id;
	
	my @dna_dna_align_feats = @{$slice->get_all_DnaAlignFeatures ('mouse_est')};
	my @dnadnafeat;
	my $ori_start = $feat->start;
	my $ori_end = $feat->end;
	my $s = 1;
	my $e = ($feat->end - $feat->start) + 1;
	$feat->start($s);
	$feat->end($e);
	
	foreach my $dnadna (@dna_dna_align_feats) {
		if ($dnadna->percent_id >= $est_perc_id ) {
			
			my $feature = Bio::SeqFeature::Generic->new (	-display_name => $dnadna->hseqname,
									-start => $dnadna->hstart,
									-end => $dnadna->hend,
									-strand => $dnadna->hstrand);
			push @dnadnafeat ,$feature;
			
			my $col = new Bio::SeqFeature::Collection();
			my $totaladded = $col->add_features(\@dnadnafeat);
			my $amb = 0;
			if ($dnadna->hstrand ne $feat->strand){$amb = 1}
			
			my @subset = $col->features_in_range (-range=>$feat, -contain=> 0);	
			my $coverage = $est_coverage;
			foreach my $s (@subset) {	
				my $inter = $feat->intersection($s);	
				my $cov = ($inter->length/$s->length)*100;
				unless (defined $coverage){$coverage = $cov};
				if ($cov >= $coverage) {
					my $ghash;
					$ghash->{'seq_id'} = $s->display_name;
					my $start = $ori_start - $s->start;
					if ($feat->start < $s->start){$start = $ori_start + $s->start}
					my $end = $ori_start - $s->end;
					if ($feat->end < $s->end){$end = $ori_end + $s->end}
					$ghash->{'start'} = $start;
					$ghash->{'end'} = $end;
					$ghash->{'name'} = $feat->display_name;
					$ghash->{'strand'} = $dnadna->hstrand;
					$ghash->{'parent_id'} = 0;
					$ghash->{'rank'} = 0;
					$ghash->{'description'} = 'EST';
					
					my $gene_id = $self->load->load_region($ghash);
					
					$hash->{'trapclustermap_region'}{'region_id'} = $gene_id;
					$hash->{'trapclustermap_region'}{'overlap'} = $coverage;
					$hash->{'trapclustermap_region'}{'type'} = 'EST';
					$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = $amb;				
					$hash->{'trapclustermap_region'}{'number_trapblocks'} = scalar @{$trapblock};

					my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);
					#load hash
				
					$hash->{'additional'}{'processing'} = "annotate with Mouse Est DnaAlignFeatures";
					$hash->{'additional'}{'comment'} = "";
					$hash->{'additional'}{'label'} = "";
					$hash->{'additional'}{'trapcluster_id'} = $trapcluster_id;
					$hash->{'additional'}{'user'} = "mysql_dev";
					$hash->{'additional'}{'note'} = 'EST';
					$self->load->load_trapclusteradditional($hash->{'additional'});
				
					#$hash->{'trapclustermap_region'}{'trapblock_id'} = $trapblock_id;
					$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
					$yes_gene ++;
				}
			}
		}
	}
	
	if ($yes_gene) {
		my $timestmp = `date`;
		print STDERR "------ DONE annotate gene with annotate_with_ensest\n";
		return 1;
	}else{
		print STDERR "------ no EST found\n";
		#print STDERR Dumper $trapcluster;
		my $region_id = $self->annotation->region->slice($slice);

		$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapclustermap_region'}{'type'} = "genomic";
		$hash->{'trapclustermap_region'}{'overlap'} = 100;
		$hash->{'trapclustermap_region'}{'region_id'} = $region_id;

		$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
		my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_ensest\n";
		return 0;
	}
	
}


sub annotate_with_ensrepeats {
	my ($self, $trapcluster,$feat, $slice) = @_;
	
	my $trapcluster_id = $trapcluster->{'trapcluster_id'};
	my $type = "est";
	my $moltype = $trapcluster->{'mol_type'};
	my $seqdir = $trapcluster->{'sequencing_direction'};
	my $accession = $trapcluster->{'accession'};
	my $hash;
	my $yes_gene = 0;
	my $trapclustermap_id = $feat->display_name;
	my $trapblock = $self->load->fetch->get_trapclusterblock_by_trapclustermap_id($trapclustermap_id);	
	my $ori_start = $feat->start;
	my $ori_end = $feat->end;
	my $s = 1;
	my $e = ($feat->end - $feat->start) + 1;
	$feat->start($s);
	$feat->end($e);

	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		return 0;
	}
	
	$hash->{'trapclustermap_region'}{'trapclustermap_id'} = $trapclustermap_id;
	
	my @repeatfeat;

	foreach my $repeat (@{$slice->get_all_RepeatFeatures()}) {
		$debug && print STDOUT "Repeat Name: ".$repeat->seqname." - Start  : ".$repeat->start." - End:".$repeat->end."Strand:".$repeat->strand." - Score:".$repeat->score." - Hit start:".$repeat->hstart." - Hit end:".$repeat->hend." - Consensus:".$repeat->display_id."\n";
		
		### 
		if ($repeat->display_id ne 'Low_complexity' && $repeat->display_id ne 'Simple_repeat') {
			my $feature = Bio::SeqFeature::Generic->new(	-display_name => $repeat->display_id,
									-start => $repeat->hstart,
									-end => $repeat->hend,
									-strand => $repeat->hstrand);
			push @repeatfeat ,$feature;
			my $amb = 0;
			if ($repeat->hstrand ne $feat->strand){$amb = 1}
			my $col = new Bio::SeqFeature::Collection();
			my $totaladded = $col->add_features(\@repeatfeat);
			
			foreach my $f (@repeatfeat) {
				$debug && print STDOUT "Feature: ".$f->display_name."\t".$f->start."\t".$f->end."\n";
				
				my @subset = $col->features_in_range(-range=>$f, -contain=> 0);
				my $coverage;
				foreach my $s (@subset) {
					$debug && print STDOUT "In range: " . $s->display_name . "\t". $s->start . "\t". $s->end . "\n";
					my $inter = $feat->intersection($s);	
					my $cov = ($inter->length/$s->length)*100;
					unless (defined $coverage){$coverage = $cov};
					if ($cov >= $coverage) {
						my $ghash;
						$ghash->{'seq_id'} = $s->display_name;
						my $start = $ori_start - $s->start;
						if ($feat->start < $s->start){$start = $ori_start + $s->start}
						my $end = $ori_start - $s->end;
						if ($feat->end < $s->end){$end = $ori_end + $s->end}
						$ghash->{'start'} = $start;
						$ghash->{'end'} = $end;
						$ghash->{'name'} = $feat->display_name;
						$ghash->{'strand'} = $repeat->hstrand;
						$ghash->{'parent_id'} = 0;
						$ghash->{'rank'} = 0;
						$ghash->{'description'} = 'Repeats';
						
						my $gene_id = $self->load->load_region($ghash);
						
						$hash->{'trapclustermap_region'}{'region_id'} = $gene_id;
						$hash->{'trapclustermap_region'}{'overlap'} = $coverage;
						$hash->{'trapclustermap_region'}{'type'} = 'Repeats';
						$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = $amb;				
						$hash->{'trapclustermap_region'}{'number_trapblocks'} = scalar @{$trapblock};
	
						my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);
					
						#load hash
					
						$hash->{'additional'}{'processing'} = "annotate with Repeat Features";
						$hash->{'additional'}{'comment'} = "";
						$hash->{'additional'}{'label'} = "";
						$hash->{'additional'}{'trapcluster_id'} = $trapcluster_id;
						$hash->{'additional'}{'user'} = "mysql_dev";
						$hash->{'additional'}{'note'} = 'Repeats';
						$self->load->load_trapclusteradditional($hash->{'additional'});
					
						#$hash->{'trapclustermap_region'}{'trapblock_id'} = $trapblock_id;
						$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
						$yes_gene ++;
					}
				}
			}
		}
	}	
	
	if ($yes_gene) {
		my $timestmp = `date`;
		print STDERR "------ DONE annotate gene with annotate_with_ensrepeats\n";
		return 1;
	}else{
		print STDERR "------ no Repeats found\n";
		#print STDERR Dumper $trapcluster;
		my $region_id = $self->annotation->region->slice($slice);

		$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapclustermap_region'}{'type'} = "genomic";
		$hash->{'trapclustermap_region'}{'overlap'} = 100;
		$hash->{'trapclustermap_region'}{'region_id'} = $region_id;

		$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
		my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_ensrepeats\n";
		return 0;
	}
}	
	

sub annotate_with_mouse_cdna {
	my ($self, $trapcluster,$feat, $slice) = @_;
	my $cdna_coverage = $conf{'annotation'}{'cdna_coverage'};
	my $cdna_perc_id = $conf{'annotation'}{'cdna_perc_id'};
	my $trapcluster_id = $trapcluster->{'trapcluster_id'};
	my $type = "cdna";
	my $moltype = $trapcluster->{'mol_type'};
	my $seqdir = $trapcluster->{'sequencing_direction'};
	my $accession = $trapcluster->{'accession'};
	my $hash;
	my $yes_gene = 0;
	my $trapclustermap_id = $feat->display_name;
	my $trapblock = $self->load->fetch->get_trapclusterblock_by_trapclustermap_id($trapclustermap_id);	
	
	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		return 0;
	}
	
	$hash->{'trapclustermap_region'}{'trapclustermap_id'} = $trapclustermap_id;
	
	my @dna_dna_align_feats = @{$slice->get_all_DnaAlignFeatures ('mouse_cdna')};
	my @dnadnafeat;
	my $ori_start = $feat->start;
	my $ori_end = $feat->end;
	my $s = 1;
	my $e = ($feat->end - $feat->start) + 1;
	$feat->start($s);
	$feat->end($e);
	
	foreach my $dnadna (@dna_dna_align_feats) {
		$debug && print "FOUND ".$dnadna->hseqname."\n";
		if ($dnadna->percent_id >= $cdna_perc_id) {
			$debug && print "Percentage IS RIGHT\n"; 
			my $feature = Bio::SeqFeature::Generic->new (	-display_name => $dnadna->hseqname,
									-start => $dnadna->hstart,
									-end => $dnadna->hend,
									-strand => $dnadna->hstrand);
			push @dnadnafeat ,$feature;
			
			my $col = new Bio::SeqFeature::Collection();
			my $totaladded = $col->add_features(\@dnadnafeat);
			my $amb = 0;
			if ($dnadna->hstrand ne $feat->strand){$amb = 1}
			
			my @subset = $col->features_in_range (-range=>$feat, -contain=> 0);	
			my $coverage = $cdna_coverage;
			foreach my $s (@subset) {	
				my $inter = $feat->intersection($s);	
				my $cov = ($inter->length/$s->length)*100;
				unless (defined $coverage){$coverage = $cov};
				if ($cov >= $coverage) {
					$debug && print "COVERAGE IS RIGHT\n"; 
					my $ghash;
					$ghash->{'seq_id'} = $s->display_name;
					my $start = $ori_start - $s->start;
					if ($feat->start < $s->start){$start = $ori_start + $s->start}
					my $end = $ori_start - $s->end;
					if ($feat->end < $s->end){$end = $ori_end + $s->end}
					$ghash->{'start'} = $start;
					$ghash->{'end'} = $end;
					$ghash->{'name'} = $feat->display_name;
					$ghash->{'strand'} = $dnadna->hstrand;
					$ghash->{'parent_id'} = 0;
					$ghash->{'rank'} = 0;
					$ghash->{'description'} = 'CDNA';
					
					my $gene_id = $self->load->load_region($ghash);
					
					$hash->{'trapclustermap_region'}{'region_id'} = $gene_id;
					$hash->{'trapclustermap_region'}{'overlap'} = $coverage;
					$hash->{'trapclustermap_region'}{'type'} = 'CDNA';
					$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = $amb;				
					$hash->{'trapclustermap_region'}{'number_trapblocks'} = scalar @{$trapblock};

					my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);
				
					#load hash
				
					$hash->{'additional'}{'processing'} = "annotate with Mouse CDNA DnaAlignFeatures";
					$hash->{'additional'}{'comment'} = "";
					$hash->{'additional'}{'label'} = "";
					$hash->{'additional'}{'trapcluster_id'} = $trapcluster_id;
					$hash->{'additional'}{'user'} = "mysql_dev";
					$hash->{'additional'}{'note'} = 'CDNA';
					$self->load->load_trapclusteradditional($hash->{'additional'});
				
					#$hash->{'trapclustermap_region'}{'trapblock_id'} = $trapblock_id;
					$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
					$yes_gene ++;
				}
				else{$debug && print "COVERAGE IS LOW\n"; }
			}
		}
		else{$debug && print "PERCENTAGE IS LOW\n"; }
	}
	
	
	if ($yes_gene) {
		my $timestmp = `date`;
		print STDERR "------ DONE annotate gene with annotate_with_mouse_cdna\n";
		return 1;
	}else{
		print STDERR "------ no cdna found\n";
		#print STDERR Dumper $trapcluster;
		my $region_id = $self->annotation->region->slice($slice);

		$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapclustermap_region'}{'type'} = "genomic";
		$hash->{'trapclustermap_region'}{'overlap'} = 100;
		$hash->{'trapclustermap_region'}{'region_id'} = $region_id;

		$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
		my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_mouse_cdna\n";
		return 0;
	}
	
	
}

=head2 annotate_with_ensembl_gene

	Title    : annotate_with_ensembl_gene
	Usage    : $obj->annotate_with_ensembl_gene(Args)
 	Function : annotate trapcluster wth ensembl gene
 	Returns  : 
 	Args     : trapcluster href, hit feat, slice on genome	  

=cut

sub annotate_with_ensembl_gene{
	my ($self, $trapcluster,$feat, $slice) = @_;
	my $trapcluster_id = $trapcluster->{'trapcluster_id'};
	my $type = "ensembl gene";
	my $moltype = $trapcluster->{'mol_type'};
	my $seqdir = $trapcluster->{'sequencing_direction'};
	my $accession = $trapcluster->{'accession'};
	my $hash;
	my $yes_gene = 0;
	my $trapclustermap_id = $feat->display_name;
	my $trapblock = $self->load->fetch->get_trapclusterblock_by_trapclustermap_id($trapclustermap_id);
	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		return 0;
	}
	$hash->{'trapclustermap_region'}{'trapclustermap_id'} = $trapclustermap_id;
	
	foreach my $gene (@{$slice->get_all_Genes}) {
		my $gstart = $gene->seq_region_start;
		my $gend = $gene->seq_region_end;
		my $gstrand = $gene->strand;
	    #test if the feature (trapclustermap) inteserct a gene
		my ($ti) = $self->annotation->trap_utility->intersection($gene,$gstart,$gend,$feat);
		if ($ti->{'i'}) {
			$yes_gene ++;
			my $gene_id = $self->annotation->region->gene($gene); #load gene in region table
			my $amb = $self->annotation->trap_utility->check_strand($trapcluster,$feat,$gene); #used to check if the alignment is on the correct strand of the gene
			
			if ($debug && $amb == 1){next;}
			#if the strand is not correct it isn't possible to infer the insertion
			$debug && print STDERR "------ FOUND ".$gene->biotype.": ".$gene->stable_id." and AMB = $amb\n";			

			#set generic info
			$hash->{'trapclustermap_region'}{'region_id'} = $gene_id;
			$hash->{'trapclustermap_region'}{'overlap'} = $ti->{'coverage'};
			$hash->{'trapclustermap_region'}{'type'} = $gene->biotype;
			$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = $amb;				
			$hash->{'annotation'}{'trapclustermap_id'} = $feat->display_name;
			$hash->{'trapclustermap_region'}{'number_trapblocks'} = scalar @{$trapblock};
			my $number_annotate_trapblocks = 0;

			my @regions = @{$self->annotation->region->rearrange_exons($gene,$gene_id)};
			my $first = $regions[0]; 
			my $last = $regions[$#regions];
			my $first_exon_id = $first->{'region_id'};
			my $last_exon_id = $last->{'region_id'};
			my $first_start = $first->{'region_start'}; 
			my $last_end = $last->{'region_end'};
			$hash->{'trapclustermap_region'}{'total_number_exons'} = 1+(int(scalar @regions)/2);
			my ($nti) = $self->annotation->trap_utility->intersection($regions[0],$first_start,$last_end,$feat);			
			unless ($nti->{'i'}) {
				print STDERR "------ no coding block found for gene: $gene_id\n"; 
				$hash->{'trapclustermap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
				$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
				next;
			}

			#if the gene is a non coding gene and trapcluster is RACE it is treated like a genomic
			if ($gene->biotype !~ /coding/ && $moltype eq 'mRNA' ){
				print STDERR "------ Got here because $moltype and ".$gene->biotype."\n";
				$hash->{'additional'}{'processing'} = "annotate with ensembl gene";
				$hash->{'additional'}{'comment'} = "non coding gene";
				$hash->{'additional'}{'label'} = $gene->biotype;
				$hash->{'additional'}{'trapcluster_id'} = $trapcluster_id;
				$hash->{'additional'}{'note'} = $ti->{'label'}.": ".$ti->{'comment'};
				$hash->{'trapclustermap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
				$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = 1;
				$self->load->load_trapclusteradditional($hash->{'additional'});
				$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});			
				my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);
				next;	
			}#end if non coding gene and race tag

			my $trapblock_id;
			my $prev_id = 0;
			my $post_id = 0;
			my $seen_trapped_region = 0;
			my %trapped;
			my $trapped_region_rank;

			#check trapped region and assing trapblock
			foreach my $rhref_block (@{$trapblock}) {

				my $tbstart = $rhref_block->{'start'};
				my $tbend = $rhref_block->{'end'};
				my $strand = $rhref_block->{'strand'};

				my $tbf = Bio::SeqFeature::Generic->new(	-display_name => $trapblock_id,
												-start => $tbstart,
												-end => $tbend,
												-strand => $strand	);

				my ($tif) = $self->annotation->trap_utility->intersection($gene,$gstart,$gend,$tbf);				
				if ($tif->{'i'}){$number_annotate_trapblocks ++}
				else{next}#the trapclustermap find more than one gene but this trapblock could not intersect all the genes

				my %check_rank;
				foreach my $region (@regions) {
					my $rstart = $region->{'region_start'};
					my $rend = $region->{'region_end'};
					my ($tir) = $self->annotation->trap_utility->intersection($region,$rstart,$rend,$tbf);

					my $region_id = $region->{'region_id'};


					if ($seen_trapped_region && $post_id == 0){
						# the previous region is trapped						
						$post_id = $region_id; #this is the region after the trapped one
					}

					if ($tir->{'i'}) {
						$debug && print STDERR "------ line 533 ";
						$debug && print STDERR ref $self;
						$debug && print STDERR "\n------ found an intersecant region ".$region->{'region_id'}." -> ".$region->{'seq_id'}."\n";
						next if $check_rank{$region->{'rank'}};
						$check_rank{$region->{'rank'}} ++;
						# found first intersecant region
						$seen_trapped_region ++;
						$trapped{$seen_trapped_region}{'region'} = $region;
						$trapped{$seen_trapped_region}{'posrelative'} = $tir;
						$trapped{$seen_trapped_region}{'trapblock_id'} = $rhref_block->{'trapblock_id'};
						$trapped{$seen_trapped_region}{'ttbf'} = $tbf;
						$post_id = 0; #void in order to get the real post one
						$debug && print STDERR "------ REGION $region_id RANK ".$region->{'rank'}." SEEN RANK no ".$check_rank{$region->{'rank'}}." no $seen_trapped_region \n";

					}#end trapped region
					if ($seen_trapped_region == 0){$prev_id = $region_id;}

				}#end region loop

				if ($seen_trapped_region == 0){
					print STDERR warn "------ got a problem with trapclustermap $trapclustermap_id: found no region intersecation block for gene ".$gene->stable_id."\n";
					next;
				}#end test at least one trapped region must be found
			}#end trapblock loop
			next unless $number_annotate_trapblocks;
			next unless $seen_trapped_region;
			my $first_trapped = $trapped{1}{'region'};
			my $last_trapped = $trapped{$seen_trapped_region}{'region'};
			my $first_id = $first_trapped->{'region_id'};
			my $last_id = $last_trapped->{'region_id'};
			if($post_id == 0){$post_id = $last_id} #the trapped region is the last
			if($prev_id == 0){$prev_id = $first_id} #the trapped region is the first
			$debug && print STDERR "------ PRE $prev_id POST $post_id  FIRST $first_id LAST $last_id\n";
			my $pre_region = $self->load->fetch->get_region_by_id($prev_id);
			my $post_region = $self->load->fetch->get_region_by_id($post_id);
			if ($moltype eq 'mRNA'){
				if ($pre_region->{'description'} =~ /intron/){
					my ($pre,$post) = split /_/, $pre_region->{'seq_id'};
					$pre_region = $self->load->fetch->get_region_by_seq_id($pre);
					$prev_id = $pre_region->{'region_id'};
				}
				if ($post_region->{'description'} =~ /intron/){
					my ($pre,$post) = split /_/, $post_region->{'seq_id'};
					$post_region = $self->load->fetch->get_region_by_seq_id($post);
					$post_id = $post_region->{'region_id'};
				}
			}#end test if pre or post are intron for RACE tag
			$debug && print STDERR "------ PRE $prev_id POST $post_id  FIRST $first_id LAST $last_id\n";
			#defined the trapped region according with the moltype and the strand of the gene
			if ($moltype eq 'mRNA'){
				if ($gstrand eq "-1"){							
					if ($seqdir eq "3"){$trapped_region_rank =  $seen_trapped_region}
					elsif ($seqdir eq "5"){$trapped_region_rank = 1}
				}
				if ($gstrand eq "1"){							
					if ($seqdir eq "5"){$trapped_region_rank =  $seen_trapped_region}
					elsif ($seqdir eq "3"){$trapped_region_rank = 1}
				}
			}
			elsif($moltype eq "genomic DNA"){
				if ($gstrand eq "1"){$trapped_region_rank = 1}
				if ($gstrand eq "-1"){$trapped_region_rank =  $seen_trapped_region}
			}
			else{
				print STDERR "ERROR: problem with moltype $moltype. Exit program\n";
				die;
			}

			$hash->{'annotation'}{'multiple'} = $seen_trapped_region - 1;
			#defined the flanking region
			$hash = $self->annotation->trap_utility->define_internal_trapped_region($hash,$gene,$trapped{$trapped_region_rank}{'tbf'},$first_trapped,$last_trapped,$pre_region,$post_region,$seqdir,$moltype);

			$debug && print STDERR "------ no $trapped_region_rank\n";
			$trapblock_id = $trapped{$trapped_region_rank}{'trapblock_id'};
			$debug && print STDERR "------ trapblock_id $trapblock_id\n";

			my $posrelative = $trapped{$trapped_region_rank}{'posrelative'};

			my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
			$self->load->fetch->update($update_trapclustercheck);

			#load hash

			$hash->{'insertion'}{'trapclustermap_id'} = $trapclustermap_id;
			$hash->{'insertion'}{'trapblock_id'} = $trapblock_id;
			$hash->{'insertion'}{'trapcluster_id'} = $trapcluster_id;
			$hash->{'insertion'}{'gene_id'} = $gene_id;
			$hash->{'insertion'}{'tbf'} = $feat;
			$hash->{'insertion'}{'sequencing_direction'} = $seqdir;						
			$hash->{'insertion'}{'moltype'} = $moltype;

			$hash->{'additional'}{'processing'} = "annotate with ensembl gene";
			$hash->{'additional'}{'comment'} = "";
			$hash->{'additional'}{'label'} = $self->annotation->region->get_label($gene);
			$hash->{'additional'}{'trapcluster_id'} = $trapcluster_id;
			$hash->{'additional'}{'user'} = "mysql_dev";
			$hash->{'additional'}{'note'} = $posrelative->{'label'}.": ".$posrelative->{'comment'};
			$self->load->load_trapclusteradditional($hash->{'additional'});

			$hash->{'annotation'}{'trapblock_id'} = $trapblock_id;
			$hash->{'annotation'}{'display_name'} = $trapped{$trapped_region_rank}{'region'}{'seq_id'};
			$hash->{'annotation'}{'comment'} = "" ;
			$hash->{'annotation'}{'coverage'} = $posrelative->{'coverage'};
			$hash->{'annotation'}{'rank'} = $posrelative->{'rank'};
			$hash->{'annotation'}{'dS_s'} = $posrelative->{'dS_s'};
			$hash->{'annotation'}{'dE_e'} = $posrelative->{'dE_e'};
			$hash->{'annotation'}{'label'} = $posrelative->{'label'};
			$hash->{'annotation'}{'logic_name'} = "EnsEMBL";
			$hash->{'annotation'}{'type'} = lc($trapped{$trapped_region_rank}{'region'}{'description'});

			$hash->{'annotation'}{'display_name'} = $self->annotation->trap_utility->find_insertion_site($hash->{'insertion'}) ;
			$hash->{'annotation'}{'flanking_exon_id'} = $hash->{'insertion'}{'flanking_exon_id'};
			$hash->{'annotation'}{'region_id'} = $hash->{'insertion'}{'region_id'};
			$self->load->load_trapclusterannotation($hash->{'annotation'});			

			$hash->{'trapclustermap_region'}{'trapblock_id'} = $trapblock_id;
			$hash->{'trapclustermap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
			$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});

		}#end intersect gene

	}#end gene loop

	if ($yes_gene) {
		my $timestmp = `date`;
		print STDERR "------ DONE annotate gene with annotate_with_ensembl_gene\n";
		return 1;
	}else{
		print STDERR "------ no gene found\n";
		#print STDERR Dumper $trapcluster;
		my $region_id = $self->annotation->region->slice($slice);

		$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapclustermap_region'}{'type'} = "genomic";
		$hash->{'trapclustermap_region'}{'overlap'} = 100;
		$hash->{'trapclustermap_region'}{'region_id'} = $region_id;

		$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
		my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
		$self->load->fetch->update($update_trapclustercheck);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_ensembl_gene\n";
		return 0;
	}
}


=head2 annotate_with_ensembl_estgene

	Title    : annotate_with_ensembl_estgene
	Usage    : $obj->annotate_with_ensembl_gene(Args)
 	Function : annotate trapcluster wth ensembl gene
 	Returns  : 
 	Args     : trapcluster href, hit feat, slice on genome	  

=cut

sub annotate_with_ensembl_estgene{
	my ($self, $trapcluster,$feat, $slice) = @_;
	my $trapcluster_id = $trapcluster->{'trapcluster_id'};
	my $type = "ensembl est gene";
	my $moltype = $trapcluster->{'mol_type'};
	my $seqdir = $trapcluster->{'sequencing_direction'};
	my $accession = $trapcluster->{'accession'};
	my $hash;
	my $yes_gene = 0;
	my $trapclustermap_id = $feat->display_name;
	my $trapblock = $self->load->fetch->get_trapclusterblock_by_trapclustermap_id($trapclustermap_id);
	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		return 0;
	}
	$hash->{'trapclustermap_region'}{'trapclustermap_id'} = $trapclustermap_id;
	print STDERR "------ got here 706\n";
	foreach my $gene (@{$slice->get_all_Genes}) {
		my $gstart = $gene->seq_region_start;
		my $gend = $gene->seq_region_end;
		my $gstrand = $gene->strand;

	    #test if the feature (trapclustermap) inteserct a gene
		my ($ti) = $self->annotation->trap_utility->intersection($gene,$gstart,$gend,$feat);
		print STDERR "------ got here 714\n";
		if ($ti->{'i'}) {
			$yes_gene ++;
			
			my $gene_id = $self->annotation->region->gene($gene); #load gene in region table
			print STDERR "------ got here 719\n";
			my $amb = $self->annotation->trap_utility->check_strand($trapcluster,$feat,$gene); #used to check if the alignment is on the correct strand of the gene
			#if ($debug && $amb == 1){next;}
			
			#if the strand is not correct it isn't possible to infer the insertion
			$debug && print STDERR "FOUND ".$gene->biotype.": ".$gene->stable_id." and AMB = $amb\n";			

			#set generic info
			$hash->{'trapclustermap_region'}{'region_id'} = $gene_id;
			$hash->{'trapclustermap_region'}{'overlap'} = $ti->{'coverage'};
			$hash->{'trapclustermap_region'}{'type'} = $gene->biotype;
			$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = $amb;				
			$hash->{'trapclustermap_region'}{'number_trapblocks'} = scalar @{$trapblock};
			my $number_annotate_trapblocks = 0;

			my @regions = @{$self->annotation->region->rearrange_est_exons($gene,$gene_id)};
			my $first = $regions[0]; 
			my $last = $regions[$#regions];
			my $first_exon_id = $first->{'region_id'};
			my $last_exon_id = $last->{'region_id'};
			my $first_start = $first->{'region_start'}; 
			my $last_end = $last->{'region_end'};
			$hash->{'trapclustermap_region'}{'total_number_exons'} = 1+(int(scalar @regions)/2);
			my ($nti) = $self->annotation->trap_utility->intersection($regions[0],$first_start,$last_end,$feat);			
			unless ($nti->{'i'}) {
				print STDERR "no coding block found for gene: $gene_id\n"; 
				$hash->{'trapclustermap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
				$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
				next;
			}

			my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);

			#load hash

			$hash->{'additional'}{'processing'} = "annotate with ensembl est gene";
			$hash->{'additional'}{'comment'} = "";
			$hash->{'additional'}{'label'} = $self->annotation->region->get_label($gene);
			$hash->{'additional'}{'trapcluster_id'} = $trapcluster_id;
			$hash->{'additional'}{'user'} = "mysql_dev";
			#$hash->{'additional'}{'note'} = $posrelative->{'label'}.": ".$posrelative->{'comment'};
			$self->load->load_trapclusteradditional($hash->{'additional'});

			#$hash->{'trapclustermap_region'}{'trapblock_id'} = $trapblock_id;
			$hash->{'trapclustermap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
			$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});

		}#end intersect gene
		print STDERR "------ got here 764\n";
	}#end gene loop

	if ($yes_gene) {
		my $timestmp = `date`;
		print STDERR "------ DONE annotate gene with annotate_with_ensembl_estgene\n";
		return 1;
	}else{
		print STDERR "------ no estgene found\n";
		#print STDERR Dumper $trapcluster;
		my $region_id = $self->annotation->region->slice($slice);

		$hash->{'trapclustermap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapclustermap_region'}{'type'} = "genomic";
		$hash->{'trapclustermap_region'}{'overlap'} = 100;
		$hash->{'trapclustermap_region'}{'region_id'} = $region_id;

		$self->load->load_trapclustermap_region($hash->{'trapclustermap_region'});
		my $update_trapclustercheck = qq{UPDATE trapclustercheck SET annotated = 1, checked = 1, mapped = 1 WHERE trapcluster_id = $trapcluster_id};
				$self->load->fetch->update($update_trapclustercheck);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_ensembl_estgene\n";
		return 0;
	}
}

sub best_overlap{
	my ($self) = @_;
=pod to implement		
		### Select the best overlapping foreach trapclusterblock and for the total trapcluster as well
		### Options are ranked as 1) refseq, 2) gene (means ensembl-gene), 3) fantom3, 4) unigene.
		$self->choose_best_overlapping_for_tcl (\@features, $trapcluster_id, $strand);
	} else {
		#### retrieving all the blocks of the maxiclusters overlapping all the blocks of the traps retrieved previously!
		my $sql2 = "select mb.maxiclusterblock_id, mb.start, mb.end from maxiclusterblock mb where mb.maxiclustermap_id = '$maxiclustermap_id' and $where2 ";
		my $sth2 = $trapdb->prepare($sql2);
		$debug && print STDOUT "SQL CODE: $sql2\n";
		$sth2->execute;
		
		my %new_blocks;
		
		while (my $block_info = $sth2->fetchrow_hashref()) {
			my $maxiclusterblock_id = $block_info->{'maxiclusterblock_id'};
			my $start = $block_info->{'start'};
			my $end = $block_info->{'end'};
			
			if ($start && $end) {
				$new_blocks{$maxiclusterblock_id}{'start'} = $start;
				$new_blocks{$maxiclusterblock_id}{'end'} = $end;
			}
		}
		
		&find_overlapping_traps ($maxicluster_id, $maxiclustermap_id, \%new_blocks, $j, $hit_db, $chr, $strand, $debug);
=cut

}

sub choose_best_overlapping_for_tcl () {
    my ($self,$fs, $trapcluster_id, $strand) = @_;
    
=pod to be implemented    
    
    my $tcl_total_overlapping;
    my $tcl_total_nested;
    
    my @features = @{$fs};
    foreach my $f (@features) {
	################ refseq and gene overlapping
	my $sth = $trapdb->prepare("select ensmusg, refseq_id from trapclusterblock_ensg_annotation where trapclusterblock_id = '".$f->display_name."' and ensmuse is not null and strand='$strand' order by refseq_id desc");
	$debug && print STDOUT "SQL CODE: select ensmusg, refseq_id from trapclusterblock_ensg_annotation where trapclusterblock_id = '".$f->display_name."' and ensmuse is not null and strand='$strand' order by refseq_id desc\n";
	my $nres = $sth->execute;
	
	if ($nres >= 1) {
		my $rhref = $sth->fetchrow_hashref;
		my $ensmusg = $rhref->{'ensmusg'};
		my $refseq_id = $rhref->{'refseq_id'};
		
		if ($refseq_id ne '') {
			my $sql = "update trapclusterblock set overlapping = 'refseq' where trapclusterblock_id = ".$f->display_name;
			$debug && print STDOUT "SQL CODE: $sql\n";
			my $sth_update_block_overlapping = $trapdb->prepare($sql);
			$sth_update_block_overlapping->execute() || die "insert failed : $DBI::errstr";
			
			$tcl_total_overlapping = 'refseq';
		} else {
			my $sql = "update trapclusterblock set overlapping = 'gene' where trapclusterblock_id = ".$f->display_name;
			$debug && print STDOUT "SQL CODE: $sql\n";
			my $sth_update_block_overlapping = $trapdb->prepare($sql);
			$sth_update_block_overlapping->execute() || die "insert failed : $DBI::errstr";
			
			if ($tcl_total_overlapping ne 'refseq') {
				$tcl_total_overlapping = 'gene';
			}
		}
		
		next;
	} else {
		### if this trapclusterblock does not overlap any exon of a refseq, check if it is intronic
		################ refseq only
		my $sth_intronic = $trapdb->prepare("select ensmusg, refseq_id from trapclusterblock_ensg_annotation where trapclusterblock_id = '".$f->display_name."' and ensmuse is null and refseq_id != '' and strand = '$strand' order by refseq_id desc");
		$debug && print STDOUT "SQL CODE: select ensmusg, refseq_id from trapclusterblock_ensg_annotation where trapclusterblock_id = '".$f->display_name."' and ensmuse is null and refseq_id != '' and strand = '$strand' order by refseq_id desc\n";
		my $nres_intronic = $sth_intronic->execute;
		
		if ($nres_intronic >= 1) {
			my $rhref_intronic = $sth_intronic->fetchrow_hashref;
			my $ensmusg_intronic = $rhref_intronic->{'ensmusg'};
			my $refseq_id_intronic = $rhref_intronic->{'refseq_id'};
			
			my $sql_intronic_update = "update trapclusterblock set intronic = 1 where trapclusterblock_id = ".$f->display_name;
			$debug && print STDOUT "SQL CODE: $sql_intronic_update\n";
			my $sth_intronic_update = $trapdb->prepare($sql_intronic_update);
			$sth_intronic_update->execute() || die "insert failed : $DBI::errstr";
			
			$tcl_total_nested = '1';
		}
	}
	
	################ fantom3
	my $sth2 = $trapdb->prepare("select display_name from trapclusterblock_fantom_annotation where trapclusterblock_id = '".$f->display_name."' and strand='$strand'");
	$debug && print STDOUT "SQL CODE: select display_name from trapclusterblock_fantom_annotation where trapclusterblock_id = '".$f->display_name."' and strand='$strand'\n";
	my $nres2 = $sth2->execute;
	
	if ($nres2 >= 1) {
		my $rhref2 = $sth2->fetchrow_hashref;
		my $fantom_name = $rhref2->{'display_name'};
		
		my $sql = "update trapclusterblock set overlapping = 'fantom3' where trapclusterblock_id = '".$f->display_name."'";
		$debug && print STDOUT "SQL CODE: $sql\n";
		my $sth_update_block_overlapping = $trapdb->prepare($sql);
		$sth_update_block_overlapping->execute() || die "insert failed : $DBI::errstr";
		
		if ($tcl_total_overlapping ne 'refseq' && $tcl_total_overlapping ne 'gene') {
			$tcl_total_overlapping = 'fantom3';
		}
		
		next;
	}
	
	################ unigene
	my $sth3 = $trapdb->prepare("select display_name from trapclusterblock_unigene_annotation where trapclusterblock_id = '".$f->display_name."'  and strand='$strand'");
	$debug && print STDOUT "SQL CODE: select display_name from trapclusterblock_unigene_annotation where trapclusterblock_id = '".$f->display_name."' and strand='$strand'\n";
	my $nres3 = $sth3->execute;
	
	if ($nres3 >= 1) {
		my $rhref3 = $sth3->fetchrow_hashref;
		my $unigene_name = $rhref3->{'display_name'};
		
		my $sql = "update trapclusterblock set overlapping = 'unigene' where trapclusterblock_id = '".$f->display_name."'";
		$debug && print STDOUT "SQL CODE: $sql\n";
		my $sth_update_block_overlapping = $trapdb->prepare($sql);
		$sth_update_block_overlapping->execute() || die "insert failed : $DBI::errstr";
		
		if ($tcl_total_overlapping ne 'refseq' && $tcl_total_overlapping ne 'gene' && $tcl_total_overlapping ne 'fantom3') {
			$tcl_total_overlapping = 'unigene';
		}
		
		next;
	}
    }
    
    my $sql_update_tcl = "update trapcluster set overlapping = '".$tcl_total_overlapping."' where trapcluster_id = $trapcluster_id";
    $debug && print STDOUT "SQL CODE: $sql_update_tcl\n";
    my $sth_update_tcl = $trapdb->prepare($sql_update_tcl);
    $sth_update_tcl->execute() || die "insert failed : $DBI::errstr";
    
    ### If the trapcluster does not overlap at least one exon of a refseq and it has at least one trapclusterblock intronic to a refseq
    if ($tcl_total_overlapping ne 'refseq' && $tcl_total_nested == 1) {
	    my $sql_update_tcl2 = "update trapcluster set intronic='".$tcl_total_nested."' where trapcluster_id = $trapcluster_id";
	    $debug && print STDOUT "SQL CODE: $sql_update_tcl2\n";
	    my $sth_update_tcl2 = $trapdb->prepare($sql_update_tcl2);
	    $sth_update_tcl2->execute() || die "insert failed : $DBI::errstr";
    }
=cut    
    
    
}


1;