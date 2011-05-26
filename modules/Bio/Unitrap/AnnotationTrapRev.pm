#!/usr/bin/perl -w

=head1 Bio::Unitrap::AnnotationTrap

=head2 Authors

=head3 Created by

             Vincenza Maselli
             v.maselli@ucl.ac.uk

=head3 Original script by

              Guglielmo Roma
              guglielmoroma@libero.it

=head2 Description
        
             This module run blast or read information from blastoutput or from database and load the tables `trapblock_annotation` and `region`	
             
=head2 Usage

	    my $obj = Bio::Unitrap::AnnotationTrap->new;
            
=cut

package Bio::Unitrap::AnnotationTrap;

use strict;
use DBI;
use Carp;
use Data::Dumper;
use vars qw(@ISA);
use File::Spec;
use Bio::Unitrap::Utils::Argument qw(rearrange);
use Bio::Unitrap::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Intron;

@ISA = qw(Bio::Root::Root);

#setting global variables

require "$ENV{'Unitrap'}/unitrap_conf.pl";

my %conf =  %::conf;
my $debug = $conf{'global'}{'debug'};

use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Unitrap::Utils::File;
use Bio::Unitrap::LoadTrap;
use Bio::Unitrap::Region;
use Bio::Unitrap::Utils::TrapUtility;

=head2 new

  Arg [..]: Take a set of named argumnts from a config file
  Example: my $retrieve = Bio::Unitrap::RetrieveTrap->new
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
   		$registry,
   		$dba_core,
   		$dba_est,
   		$slice_core_adaptor,
   		$slice_est_adaptor
    )
    = $self->_rearrange( [
      	'REGISTRY',
      	'DBACORE',
      	'DBAEST',
      	'SLICECOREADPT',
      	'SLICEESTADPT'
    ],
    @_
    );
	
	my $enshost = $conf{'dbaccess'}{'enshost'};
	my $ensuser = $conf{'dbaccess'}{'ensuser'};
	my $ensdbname = $conf{'dbaccess'}{'ensdbname'};
	my $ensestdbname = $conf{'dbaccess'}{'ensestdbname'};
	my $ensdbport = $conf{'dbaccess'}{'ensport'};
	
	unless (defined $registry){
		$registry = 'Bio::EnsEMBL::Registry';
		$registry->load_registry_from_db(
	      -host => $enshost,
	      -user => $ensuser,
	      -port => $ensdbport,
	      -species => 'Mus musculus',
	      -verbose => '0'
		);
	}
	unless (defined $slice_core_adaptor){
		$slice_core_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' ) || die $!;
	}
	unless (defined $slice_est_adaptor){
		$slice_est_adaptor = $registry->get_adaptor( 'Mouse', 'Otherfeatures', 'Slice' ) || die $!;
	}
	
	
	$registry && $self->registry($registry);
	$slice_core_adaptor && $self->slicecoreadpt($slice_core_adaptor);
	$slice_est_adaptor && $self->sliceestadpt($slice_est_adaptor);
	my $load = Bio::Unitrap::LoadTrap->new;
	$load && $self->load($load);
	my $utility = Bio::Unitrap::Utils::TrapUtility->new(-LOAD=>$load);
	$utility && $self->trap_utility($utility);
	my $region = Bio::Unitrap::Region->new(
											-LOAD=>$load,
											-REGISTRY =>$registry);
	$region && $self->region($region);
	
	my $do_ensgene = $conf{'annotation'}{'do_ensgene'};
	my $do_ensestgene = $conf{'annotation'}{'do_ensestgene'};
	my $do_genescan = $conf{'annotation'}{'do_genescan'};
	my $do_unigene = $conf{'annotation'}{'do_unigene'};
	my $unigene_coverage = $conf{'annotation'}{'unigene_coverage'};
	my $unigene_perc_id = $conf{'annotation'}{'unigene_perc_id'};
	my $do_mouse_cdna = $conf{'annotation'}{'do_mouse_cdna'};
	my $cdna_coverage = $conf{'annotation'}{'cdna_coverage'};
	my $cdna_perc_id = $conf{'annotation'}{'cdna_perc_id'};
	my $do_ensest = $conf{'annotation'}{'do_ensest'};
	my $est_coverage = $conf{'annotation'}{'est_coverage'};
	my $est_perc_id = $conf{'annotation'}{'est_perc_id'};
	my $do_tclg = $conf{'annotation'}{'do_tclg'};
	my $tclg_coverage = $conf{'annotation'}{'tclg_coverage'};
	my $do_ensrepeats = $conf{'annotation'}{'do_ensrepeats'};
	
	$do_ensgene && $self->do_ensgene($do_ensgene);
		
  	return $self;
}

# sub subname{
# 
#   my ($self, $value) = @_;
#   $self->{'subname'} = $value if defined $value;
#   
#   return $self->{'subname'};
# }

sub registry{

  my ($self, $value) = @_;
  $self->{'registry'} = $value if defined $value;
  
  return $self->{'registry'};
}



sub slicecoreadpt{

  my ($self, $value) = @_;
  $self->{'slicecoreadpt'} = $value if defined $value;
  return $self->{'slicecoreadpt'};
}

sub sliceestadpt{

  my ($self, $value) = @_;
  $self->{'sliceestadpt'} = $value if defined $value;
  
  return $self->{'sliceestadpt'};
}


=head2 load

 Title    : load
 Usage    : $obj->load([$newval])
 Function : get/set method for attribute load 
 Returns  : value of load
 Args     : newval of load (optional)

=cut 

sub load{

  my ($self, $value) = @_;
  $self->{'load'} = $value if defined $value;
  
  return $self->{'load'};
}

=head2 region

 Title    : region
 Usage    : $obj->region([$newval])
 Function : get/set method for attribute region 
 Returns  : value of region
 Args     : newval of region (optional)

=cut 

sub region{

  my ($self, $value) = @_;
  $self->{'region'} = $value if defined $value;
  
  return $self->{'region'};
}


=head2 trap_utility

 Title    : trap_utility
 Usage    : $obj->trap_utility([$newval])
 Function : get/set method for attribute trap_utility 
 Returns  : value of trap_utility
 Args     : newval of trap_utility (optional)

=cut 

sub trap_utility{

  my ($self, $value) = @_;
  $self->{'trap_utility'} = $value if defined $value;
  
  return $self->{'trap_utility'};
}

=head2 do_ensgene

 Title    : do_ensgene
 Usage    : $obj->do_ensgene([$newval])
 Function : get/set method for attribute do_ensgene 
 Returns  : value of do_ensgene
 Args     : newval of do_ensgene (optional)

=cut 

sub do_ensgene{

  my ($self, $value) = @_;
  $self->{'do_ensgene'} = $value if defined $value;
  
  return $self->{'do_ensgene'};
}

=head2 Get overlap with ensembl genes 
  
    (both RefSeq and ENSEMBL)
	using the slice to search gene calling annotate_with_ensembl_gene
	resolving ensembl conflicts resolve_ensembl_conflicts
	add strand to the hash
	prepare insert into trapblock_ensg_annotation

=head2 Get overlap with ensembl est-genes
  
  	using the slice to search est
  	resolving ensembl conflicts
  	adding strand to the hash (commented)
  	prepare insert into trapblock_ensestg_annotation

=head2 Get overlap with ensembl genescan_prediction
	
	using the slice to search genescan
	calling annotate_with_ensembl_genescan 
	
=head2 Get overlap with ensembl UniGene
	
=head2 Get overlap with ensembl mouse-cdnas

=head Get overlap with ensembl mouse-ests

=head2 Get overlap with repeats

  	get slices by repeats
  	checking only for more_complexity repeats
 	create new Bio::SeqFeature::Generic
  	create new Bio::SeqFeature::Collection object. This object will efficiently allow one for query subsets of ranges within a large collection of sequence features
  	using each feature to searh it into the collection
	calling relative_positioning
	calling RangeI intersecation method which gives the range that is contained by all ranges
	calling prepare_stmt and insert_set if the new coverage is more large than the previous
	
=cut

sub _get_col{
	my ($self,$feats,$col,$tag,$perc_id) = @_;
	my @feat;
	foreach my $f (@{$feats}) {		
		my $ok = 0;
		my $dis_name;
		if ($tag eq 'dna'){
		if ($f->percent_id >= $perc_id &&  $f->hseqname =~ /Mm/) {$ok = 1}
			$dis_name = $f->hseqname;
		}
		if ($tag eq 'repeat'){
			if ($f->display_id ne 'Low_complexity' && $f->display_id ne 'Simple_repeat') {$ok =1}
			$dis_name = $f->display_id;
		}
		if ($ok){
			my $feature = Bio::SeqFeature::Generic->new (	-display_name => $dis_name,
								-start => $f->hstart,
								-end => $f->hend,
								-strand => $f->hstrand);
			push @feat ,$feature;
		}
	}
	my $totaladded = $col->add_features(\@feat);
	return $col;
}
		
		
sub _get_feat{
	my ($self, $features, $coverage,$col) = @_;
	my %ins;	
	foreach my $f (@${features}) {
		my @subset = $col->features_in_range (-range=>$f, -contain=> 0);			
		foreach my $s (@subset) {				
			my $posrelative = $self->_relative_positioning ($s->display_name, $s->start, $s->end, $f->display_name, $f->start, $f->end);
			my $inter = $f->intersection($s);				
			my $cov = ($inter->length/$s->length)*100;
			unless (defined $coverage){$coverage = $cov};
			if ($cov >= $coverage) { #first time always true
				
				$ins{'annotation'}{'display_name'} = $s->display_name;
				$ins{'annotation'}{'dS_s'} = $posrelative->{'dS_s'};
				$ins{'annotation'}{'dE_e'} = $posrelative->{'dE_e'};
				$ins{'annotation'}{'label'} = $posrelative->{'label'};
				$ins{'annotation'}{'comment'} = $posrelative->{'comment'};
				$ins{'annotation'}{'exon_coverage'} = ($inter->length/$s->length)*100;
				$ins{'annotation'}{'trapblock_id'} = $f->display_name;
				$ins{'annotation'}{'strand'} = $s->strand;	
			}
		}
	}
	return \%ins;
}


sub annotate_with_ensembl_genescan {
	my ($self, $f, $slice) = @_;
	my %ins;
	foreach my $ptrans (@{$slice->get_all_PredictionTranscripts}) {
		my $ptrans_feature = Bio::SeqFeature::Generic->new (	-display_name => $ptrans->stable_id,
									-start => $ptrans->start,
									-end => $ptrans->end,
									-strand => $ptrans->strand);
		my $ig = $f->intersection($ptrans_feature);		
		my $exon_check;		
		foreach my $exon (@{$ptrans->get_all_Exons}) {
			my $exon_feature = Bio::SeqFeature::Generic->new (	-display_name => $exon->display_id,
										-start => $exon->start,
										-end => $exon->end,
										-strand => $exon->strand);
			my $i = $f->intersection($exon_feature);
			if ($i) {
				my $posrelative = $self->_relative_positioning ($exon->display_id, $exon->start, $exon->end, $f->display_name, $f->start, $f->end);				
				my $ilength = $i->length;
				my $ecoverage = ($ilength/$exon->length)*100;
				$ins{'annotation'}{'display_name'} = $ptrans->display_id;
				$ins{'annotation'}{'exon_name'} = $exon->display_id;
				$ins{'annotation'}{'trapblock_id'} = $f->display_name;
				$ins{'trapmap_region'}{'exon_coverage'} = $ecoverage;
				$ins{'annotation'}{'dS_s'} = $posrelative->{'dS_s'};
				$ins{'annotation'}{'dE_e'} = $posrelative->{'dE_e'};
				$ins{'annotation'}{'label'} = $posrelative->{'label'};
				$ins{'annotation'}{'comment'} = $posrelative->{'comment'};
				$ins{'annotation'}{'type'} = $ptrans->analysis->logic_name();
				$ins{'annotation'}{'strand'} = $ptrans->strand;
			}
		}
		unless ($exon_check) {
			$ins{'annotation'}{'display_name'} = $ptrans->display_id;
			$ins{'annotation'}{'trapblock_id'} = $f->display_name;
			$ins{'annotation'}{'label'} = 'intronic';
			$ins{'annotation'}{'strand'} = $ptrans->strand;
			$ins{'annotation'}{'type'} = $ptrans->analysis->logic_name();
		}
	}
	return \%ins;
}


=head2 annotate_with_ensembl_gene

	Title    : annotate_with_ensembl_gene
	Usage    : $obj->annotate_with_ensembl_gene(Args)
 	Function : annotate trap wth ensembl gene
 	Returns  : 
 	Args     : trap href, hit feat, slice on genome	  

=cut

sub annotate_with_ensembl_gene{
	my ($self, $trap,$tmf, $slice) = @_;
	$debug && print STDERR ref $self;
	$debug && print STDERR "->annotate_with_ensembl_gene\n";
	print STDERR Dumper $trap;
	my $trap_id = $trap->{'trap_id'};
	my $type = "ensembl gene";
	my $moltype = $trap->{'mol_type'};
	my $seqdir = $trap->{'sequencing_direction'};
	my $trapname = $trap->{'trap_name'};
	my $hash;
	my $yes_gene = 0;
	my $trapmap_id = $tmf->display_name;
	my $trapblock = $self->load->fetch->get_selected_trapblock_by_trapmap_id($trapmap_id);
			unless (scalar @{$trapblock}){
				print STDERR "-- no blocks were accepted for this map\n";#shouldn't occurr!!!
				die;
			}
	$hash->{'trapmap_region'}{'trapmap_id'} = $trapmap_id;
	foreach my $gene (@{$slice->get_all_Genes}) {
		my $gstart = $gene->seq_region_start;
		my $gend = $gene->seq_region_end;
		my $gstrand = $gene->strand;
		
	    #test if the feature (trapmap) inteserct a gene
		my ($ti) = $self->trap_utility->intersection($gene,$gstart,$gend,$tmf);
		if ($ti->{'i'}) {
			$yes_gene ++;
			my $gene_id = $self->region->gene($slice, $gene); #load gene in region table
			my $amb = $self->trap_utility->check_strand($trap,$tmf,$gene); #used to check if the alignment is on the correct strand of the gene
			if ($debug && $amb == 1){next;}
			#if the strand is not correct it isn't possible to infer the insertion
			$debug && print STDERR "FOUND ".$gene->biotype.": ".$gene->stable_id." and AMB = $amb\n";			
			
			#set generic info
			$hash->{'trapmap_region'}{'region_id'} = $gene_id;
			$hash->{'trapmap_region'}{'overlap'} = $ti->{'coverage'};
			$hash->{'trapmap_region'}{'type'} = $gene->biotype;
			$hash->{'trapmap_region'}{'annotation_ambiguous'} = $amb;				
			$hash->{'annotation'}{'trapmap_id'} = $tmf->display_name;
			$hash->{'trapmap_region'}{'number_trapblocks'} = scalar @{$trapblock};
			my $number_annotate_trapblocks = 0;
			
			my @regions = @{$self->region->rearrange_exons($slice,$gene,$gene_id)};
			my $first = $regions[0]; 
			my $last = $regions[$#regions];
			my $first_exon_id = $first->{'region_id'};
			my $last_exon_id = $last->{'region_id'};
			my $first_start = $first->{'region_start'}; 
			my $last_end = $last->{'region_end'};
			$hash->{'trapmap_region'}{'total_number_exons'} = 1+(int(scalar @regions)/2);
			my ($nti) = $self->trap_utility->intersection($regions[0],$first_start,$last_end,$tmf);			
			unless ($nti->{'i'}) {
				print STDERR "no coding block found for gene: $gene_id\n"; 
				$hash->{'trapmap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
				$self->load->load_trapmap_region($hash->{'trapmap_region'});
				next;
			}
			
			#if the gene is a non coding gene and trap is RACE it is treated like a genomic
			if ($gene->biotype !~ /coding/ && $moltype eq 'mRNA' ){
				print STDERR "Got here because $moltype and ".$gene->biotype."\n";
				$hash->{'additional'}{'processing'} = "annotate with ensembl gene";
				$hash->{'additional'}{'comment'} = "non coding gene";
				$hash->{'additional'}{'label'} = $gene->biotype;
				$hash->{'additional'}{'trap_id'} = $trap_id;
				$hash->{'additional'}{'note'} = $ti->{'label'}.": ".$ti->{'comment'};
				$hash->{'trapmap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
				$hash->{'trapmap_region'}{'annotation_ambiguous'} = 1;
				$self->load->load_trapadditional($hash->{'additional'});
				$self->load->load_trapmap_region($hash->{'trapmap_region'});			
				my $update_trapcheck = qq{UPDATE trapcheck SET annotated = 1, checked = 1, mapped = 1 WHERE trap_id = $trap_id};
				$self->load->fetch->update($update_trapcheck);
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
				
				my ($tif) = $self->trap_utility->intersection($gene,$gstart,$gend,$tbf);				
				if ($tif->{'i'}){$number_annotate_trapblocks ++}
				else{next}#the trapmap find more than one gene but this trapblock could not intersect all the genes
				
				my %check_rank;
				foreach my $region (@regions) {
					my $rstart = $region->{'region_start'};
					my $rend = $region->{'region_end'};
					my ($tir) = $self->trap_utility->intersection($region,$rstart,$rend,$tbf);

					my $region_id = $region->{'region_id'};
					
					
					if ($seen_trapped_region && $post_id == 0){
						# the previous region is trapped						
						$post_id = $region_id; #this is the region after the trapped one
					}

					if ($tir->{'i'}) {
						$debug && print STDERR "line 622 ";
						$debug && print STDERR ref $self;
						$debug && print STDERR "\nfound an intersecant region ".$region->{'region_id'}." -> ".$region->{'seq_id'}."\n";
						next if $check_rank{$region->{'rank'}};
						$check_rank{$region->{'rank'}} ++;
						# found first intersecant region
						$seen_trapped_region ++;
						$trapped{$seen_trapped_region}{'region'} = $region;
						$trapped{$seen_trapped_region}{'posrelative'} = $tir;
						$trapped{$seen_trapped_region}{'trapblock_id'} = $rhref_block->{'trapblock_id'};
						$trapped{$seen_trapped_region}{'ttbf'} = $tbf;
						$post_id = 0; #void in order to get the real post one
						$debug && print STDERR "REGION $region_id RANK ".$region->{'rank'}." SEEN RANK no ".$check_rank{$region->{'rank'}}." no $seen_trapped_region \n";
						
					}#end trapped region
					if ($seen_trapped_region == 0){$prev_id = $region_id;}
				
				}#end region loop
				
				if ($seen_trapped_region == 0){
					print STDERR warn "got a problem with TRAPMAP $trapmap_id: found no region intersecation block for gene ".$gene->stable_id."\n";
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
			$debug && print STDERR "PRE $prev_id POST $post_id  FIRST $first_id LAST $last_id\n";
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
			$debug && print STDERR "PRE $prev_id POST $post_id  FIRST $first_id LAST $last_id\n";
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
				print STDERR "problem with moltype $moltype\n";
				die;
			}
			
			$hash->{'annotation'}{'multiple'} = $seen_trapped_region - 1;
			#defined the flanking region
			$hash = $self->trap_utility->define_internal_trapped_region($hash,$gene,$trapped{$trapped_region_rank}{'tbf'},$first_trapped,$last_trapped,$pre_region,$post_region,$seqdir,$moltype);
			
			$debug && print STDERR "----- no $trapped_region_rank\n";
			$trapblock_id = $trapped{$trapped_region_rank}{'trapblock_id'};
			$debug && print STDERR "----- trapblock_id $trapblock_id\n";
			
			my $posrelative = $trapped{$trapped_region_rank}{'posrelative'};
			
			my $update_trapcheck = qq{UPDATE trapcheck SET annotated = 1, checked = 1, mapped = 1 WHERE trap_id = $trap_id};
			$self->load->fetch->update($update_trapcheck);
			
			#load hash
			
			$hash->{'insertion'}{'trapmap_id'} = $trapmap_id;
			$hash->{'insertion'}{'trapblock_id'} = $trapblock_id;
			$hash->{'insertion'}{'trap_id'} = $trap_id;
			$hash->{'insertion'}{'gene_id'} = $gene_id;
			$hash->{'insertion'}{'tbf'} = $tmf;
			$hash->{'insertion'}{'sequencing_direction'} = $seqdir;						
			$hash->{'insertion'}{'moltype'} = $moltype;
			
			$hash->{'additional'}{'processing'} = "annotate with ensembl gene";
			$hash->{'additional'}{'comment'} = "";
			$hash->{'additional'}{'label'} = $self->region->get_label($gene);
			$hash->{'additional'}{'trap_id'} = $trap_id;
			$hash->{'additional'}{'user'} = "unitrap";
			$hash->{'additional'}{'note'} = $posrelative->{'label'}.": ".$posrelative->{'comment'};
			$self->load->load_trapadditional($hash->{'additional'});
			
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
			
			$hash->{'annotation'}{'display_name'} = $self->trap_utility->find_insertion_site($hash->{'insertion'}) ;
			$hash->{'annotation'}{'flanking_exon_id'} = $hash->{'insertion'}{'flanking_exon_id'};
			$hash->{'annotation'}{'region_id'} = $hash->{'insertion'}{'region_id'};
			$self->load->load_annotation($hash->{'annotation'});			
			
			$hash->{'trapmap_region'}{'trapblock_id'} = $trapblock_id;
			$hash->{'trapmap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
			$self->load->load_trapmap_region($hash->{'trapmap_region'});
		
		}#end intersect gene
	
	}#end gene loop
	
	if ($yes_gene) {
		my $timestmp = `date`;
		print STDERR "-- DONE annotate gene with annotate_with_ensembl_gene\n$timestmp";
		return 1;
	}else{
		print STDERR "no gene found\n";
		print STDERR Dumper $trap;
		my $region_id = $self->region->slice($slice);
		
		$hash->{'trapmap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapmap_region'}{'type'} = "genomic";
		$hash->{'trapmap_region'}{'overlap'} = 100;
		$hash->{'trapmap_region'}{'region_id'} = $region_id;
		
		$self->load->load_trapmap_region($hash->{'trapmap_region'});
		my $update = qq{UPDATE trapcheck SET annotated = 0, checked = 1 WHERE trap_id = $trap_id};
		$self->load->fetch->update($update);
		my $timestmp = `date`;
		print STDERR "-- DONE annotate genomic region with annotate_with_ensembl_gene\n$timestmp";
		return 0;
	}
}

sub annotate {
	my ($self,$id, $idcount) = @_;
	my $lastid = $id + $idcount - 1;
	my $timestmp = `date`;
	print STDERR "-- START annotate for ID $id to $lastid\n";
	foreach (my $trap_id = $id; $trap_id <= $lastid; $trap_id ++) {
		my $trap = $self->load->fetch->get_trap_by_id($trap_id);
		print Dumper $trap;
		my @trapmap = @{$self->load->fetch->get_trapmap_by_trap_id($trap_id)};
		print STDERR "-- anntrap.pl FOUND ".scalar @trapmap." mapping\n";
		foreach my $rhref (@trapmap) {
			my $region;
			if ($rhref->{'hit_id'} =~ /NT/) {$region = "supercontig";} 
			else {$region = "chromosome";}
			$rhref->{'hit_id'} =~ s/chr//;
			my $slice_core_adaptor = $self->slicecoreadpt;
			my $slice_core = $slice_core_adaptor->fetch_by_region($region,$rhref->{'hit_id'}, $rhref->{'start'}, $rhref->{'end'});
			next unless defined $slice_core;
			my $tmfeature = Bio::SeqFeature::Generic->new(	-display_name => $rhref->{'trapmap_id'},
												-start => $rhref->{'start'},
												-end => $rhref->{'end'},
												-strand => $rhref->{'strand'}	);
			
			#$debug && print STDERR "$region,$rhref->{'hit_id'}, $rhref->{'start'}, $rhref->{'end'}\n";
			if ($self->do_ensgene) {
				
				my $found = $self->annotate_with_ensembl_gene($trap,$tmfeature, $slice_core)
				
			}
	# 		if ($self->do_ensestgene && $found == 0){
	# 			$found = $self->annotate_with_ensembl_gene($tmfeature,\@tbfs,$slice_est);
	# 			foreach my $selfhash(@{$selfarray}){
	# 				foreach my $gene_id (keys %{$selfhash}){
	# 					$selfhash->{'annotation'}{'trapblock_id'} = $trapblock_id;
	# 					$selfhash->{'annotation'}{'trapmap_id'} = $trapmap_id;
	# 					$selfhash->{'annotation'}{'trap_id'} = $trap_id;
	# 					print STDERR "line 170 anntrap.pl load->_load_annotation(annhash)\n";
	# 					$self->load->load_annotation($selfhash);
	# 					push(@region_array,$selfhash->{'trapmap_region'});
	# 				}
	# 			}
	# 		}
	# 		
	# 		my $perc_id;
	# 		my $coverage;
	# 		
	# 		if ($self->do_genescan && $found == 0) {
	# 			print STDERR "annotate_with_ensembl_genescan\n";
	# 			$found = $self->annotate_with_ensembl_genescan ($tmfeature,\@tbfs,$slice_core);
	# 			foreach my $gene_id (keys %{$selfhash}){
	# 				$selfhash->{'annotation'}{'logic_name'} = 'GENESCAN';
	# 				$selfhash->{'annotation'}{'trapblock_id'} = $trapblock_id;
	# 				$selfhash->{'annotation'}{'trapmap_id'} = $trapmap_id;
	# 				$selfhash->{'annotation'}{'trap_id'} = $trap_id;
	# 				push(@region_array,$selfhash->{'trapmap_region'});
	# 				$self->load->load_annotation($selfhash);
	# 			}
	# 		}	
	# 		if ($self->do_unigene && $found == 0) {
	# 			print STDERR "annotate_unigene\n";
	# 			my @feats = @{$slice_core->get_all_DnaAlignFeatures ('Unigene')};
	# 			$perc_id = $unigene_perc_id;
	# 			$coverage = $unigene_coverage;
	# 			$col = $self->_get_col(\@feats,$col,'dna',$perc_id);
	# 			foreach my $feature(@tbfs){
	# 				$selfhash = $self->_get_feat($feature, $col);
	# 				foreach my $gene_id (keys %{$selfhash}){
	# 					$selfhash->{'annotation'}{'logic_name'} = 'UNIGENE';
	# 					$selfhash->{'annotation'}{'trapblock_id'} = $trapblock_id;
	# 					$selfhash->{'annotation'}{'trapmap_id'} = $trapmap_id;
	# 					$selfhash->{'annotation'}{'trap_id'} = $trap_id;
	# 					push(@region_array,$selfhash->{'trapmap_region'});
	# 					$self->load->load_annotation($selfhash);
	# 				}
	# 			}
	# 		}
	# 		if ($self->do_mouse_cdna && $found == 0) {
	# 			print STDERR "annotate_cdna\n";
	# 			my @feats = @{$slice_core->get_all_DnaAlignFeatures ('mouse_cdna')};
	# 			$perc_id = $cdna_perc_id;
	# 			$coverage = $cdna_coverage;
	# 			$col = $self->_get_col(\@feats,$col,'dna',$perc_id);
	# 			foreach my $feature(@tbfs){
	# 			$selfhash = $self->_get_feat($feature, $coverage,$col);
	# 				foreach my $gene_id (keys %{$selfhash}){
	# 					$selfhash->{'annotation'}{'logic_name'} = 'CDNA';
	# 					$selfhash->{'annotation'}{'trapblock_id'} = $trapblock_id;
	# 					$selfhash->{'annotation'}{'trapmap_id'} = $trapmap_id;
	# 					$selfhash->{'annotation'}{'trap_id'} = $trap_id;
	# 					push(@region_array,$selfhash->{'trapmap_region'});
	# 					$self->load->load_annotation($selfhash);
	# 				}
	# 			}
	# 		}
	# 		if ($self->do_ensest && $found == 0) {
	# 			print STDERR "annotate_ensest\n";
	# 			my @feats = @{$slice_est->get_all_DnaAlignFeatures ('mouse_est')};
	# 			$perc_id = $est_perc_id;
	# 			$coverage = $est_coverage;
	# 			$col = $self->_get_col(\@feats,$col,'dna',$perc_id);
	# 			foreach my $feature(@tbfs){
	# 				$selfhash = $self->_get_feat($feature, $coverage,$col);
	# 				foreach my $gene_id (keys %{$selfhash}){
	# 					$selfhash->{'annotation'}{'logic_name'} = 'EST';
	# 					$selfhash->{'annotation'}{'trapblock_id'} = $trapblock_id;
	# 					$selfhash->{'annotation'}{'trapmap_id'} = $trapmap_id;
	# 					$selfhash->{'annotation'}{'trap_id'} = $trap_id;
	# 					push(@region_array,$selfhash->{'trapmap_region'});
	# 					$self->load->load_annotation($selfhash);
	# 				}
	# 			}
	# 		}
	# 		if ($self->do_ensrepeats && $found == 0) {
	# 			print STDERR "annotate_repeats\n";
	# 			my @feats = @{$slice_core->get_all_RepeatFeatures()};
	# 			$perc_id = undef;
	# 			$col = $self->_get_col(\@feats,$col,'repeat',$perc_id);
	# 			foreach my $feature(@tbfs){
	# 				$selfhash = $self->_get_feat($feature, $coverage,$col);
	# 				foreach my $gene_id (keys %{$selfhash}){
	# 					$selfhash->{'annotation'}{'logic_name'} = 'REPEAT';
	# 					$selfhash->{'annotation'}{'trapblock_id'} = $trapblock_id;
	# 					$selfhash->{'annotation'}{'trapmap_id'} = $trapmap_id;
	# 					$selfhash->{'annotation'}{'trap_id'} = $trap_id;
	# 					push(@region_array,$selfhash->{'trapmap_region'});
	# 					$self->load->load_annotation($selfhash);
	# 				}
	# 			}
	# 		}
		}
	}
	$timestmp = `date`;
	print STDERR "$timestmp\n-- DONE annotate\n";
	return 1;
}

1;
