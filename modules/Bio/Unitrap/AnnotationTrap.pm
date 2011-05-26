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
	my $do_mouse_cdna = $conf{'annotation'}{'do_mouse_cdna'};
	my $do_ensest = $conf{'annotation'}{'do_ensest'};
	my $do_tclg = $conf{'annotation'}{'do_tclg'};
	my $tclg_coverage = $conf{'annotation'}{'tclg_coverage'};
	my $do_ensrepeats = $conf{'annotation'}{'do_ensrepeats'};
	
	$do_ensgene && $self->do_ensgene($do_ensgene);
	$do_ensestgene && $self->do_ensestgene($do_ensestgene);
	$do_genescan && $self->do_genescan($do_genescan);
	$do_unigene && $self->do_unigene($do_unigene);
	$do_mouse_cdna && $self->do_mouse_cdna($do_mouse_cdna);
	$do_ensest && $self->do_ensest($do_ensest);
	$do_ensrepeats && $self->do_ensrepeats($do_ensrepeats);
	
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





=head2 do_genescan

 Title    : do_genescan
 Usage    : $obj->do_genescan([$newval])
 Function : get/set method for attribute do_genescan 
 Returns  : value of do_genescan
 Args     : newval of do_genescan (optional)

=cut 

sub do_genescan{

  my ($self, $value) = @_;
  $self->{'do_genescan'} = $value if defined $value;
  
  return $self->{'do_genescan'};
}

=head2 do_ensest

 Title    : do_ensest
 Usage    : $obj->do_ensest([$newval])
 Function : get/set method for attribute do_ensest 
 Returns  : value of do_ensest
 Args     : newval of do_ensest (optional)

=cut 

sub do_ensest{

  my ($self, $value) = @_;
  $self->{'do_ensest'} = $value if defined $value;
  
  return $self->{'do_ensest'};
}

=head2 ensrepeats

 Title    : ensrepeats
 Usage    : $obj->do_ensrepeats([$newval])
 Function : get/set method for attribute do_ensrepeats 
 Returns  : value of do_ensrepeats
 Args     : newval of do_ensrepeats (optional)

=cut 

sub do_ensrepeats{

  my ($self, $value) = @_;
  $self->{'do_ensrepeats'} = $value if defined $value;
  
  return $self->{'do_ensrepeats'};
}

=head2 do_mouse_cdna

 Title    : do_mouse_cdna
 Usage    : $obj->do_mouse_cdna([$newval])
 Function : get/set method for attribute do_mouse_cdna 
 Returns  : value of do_mouse_cdna
 Args     : newval of do_mouse_cdna (optional)

=cut 

sub do_mouse_cdna{

  my ($self, $value) = @_;
  $self->{'do_mouse_cdna'} = $value if defined $value;
  
  return $self->{'do_mouse_cdna'};
}


=head2 do_ensestgene

 Title    : do_ensestgene
 Usage    : $obj->do_ensestgene([$newval])
 Function : get/set method for attribute do_ensestgene 
 Returns  : value of do_ensestgene
 Args     : newval of do_ensestgene (optional)

=cut 

sub do_ensestgene{

  my ($self, $value) = @_;
  $self->{'do_ensestgene'} = $value if defined $value;
  
  return $self->{'do_ensestgene'};
}

=head2 do_unigene

 Title    : do_unigene
 Usage    : $obj->do_unigene([$newval])
 Function : get/set method for attribute do_unigene 
 Returns  : value of do_unigene
 Args     : newval of do_unigene (optional)

=cut 

sub do_unigene{

  my ($self, $value) = @_;
  $self->{'do_unigene'} = $value if defined $value;
  
  return $self->{'do_unigene'};
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
	$debug && print STDERR "got in _get_col line 354\n";
	my @feat;
	foreach my $f (@{$feats}) {	
		my $ok = 0;
		my $dis_name;
		if ($tag eq 'dna'){
			$debug && print STDERR "%ID ",$f->percent_id," HSEQNAME ",$f->hseqname,"\n";
			if ($f->percent_id >= $perc_id &&  $f->hseqname =~ /Mm/) {$ok = 1}
			$dis_name = $f->hseqname;
		}
		if ($tag eq 'repeat'){
			if ($f->display_id ne 'Low_complexity' && $f->display_id ne 'Simple_repeat') {$ok =1}
			$dis_name = $f->display_id;
		}
		if ($ok){
			$debug && print STDERR "Found a feat\n";
			
			my $feature = Bio::SeqFeature::Generic->new (	-display_name => $dis_name,
								-start => $f->hstart,
								-end => $f->hend,
								-strand => $f->hstrand);
			push @feat ,$feature;
			$col->add_features(\@feat);
		}
	}
	#my $totaladded = $col->add_features(\@feat);
	
	return $col;
}
		
		



sub annotate_with_ensembl_genescan {
	my ($self, $trap,$tmf, $slice) = @_;
	my $trap_id = $trap->{'trap_id'};
	my $type = "genescan prediction";
	my $moltype = $trap->{'mol_type'};
	my $seqdir = $trap->{'sequencing_direction'};
	my $trapname = $trap->{'trap_name'};
	my $hash;
	my $yes_gene = 0;
	my $trapmap_id = $tmf->display_name;
	my $trapblock = $self->load->fetch->get_selected_trapblock_by_trapmap_id($trapmap_id);
	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		die;
	}
	
	#WHICH IS COLLED GENE HERE IS ACTUALLY A PREDICT TRANSCRIPT
	$hash->{'trapmap_region'}{'trapmap_id'} = $trapmap_id;
	print STDERR "------ got here 706\n";
	foreach my $gene (@{$slice->get_all_PredictionTranscripts}) {
		my $gstart = $gene->seq_region_start;
		my $gend = $gene->seq_region_end;
		my $gstrand = $gene->strand;

	    #test if the feature (trapmap) inteserct a gene
		my ($ti) = $self->trap_utility->intersection($gene,$gstart,$gend,$tmf);
		print STDERR "------ got here 714\n";
		if ($ti->{'i'}) {
			$yes_gene ++;
			
			my $gene_id = $self->region->gene($gene); #load gene in region table
			print STDERR "------ got here 719\n";
			my $amb = $self->trap_utility->check_strand($trap,$tmf,$gene); #used to check if the alignment is on the correct strand of the gene
			#if ($debug && $amb == 1){next;}
			
			#if the strand is not correct it isn't possible to infer the insertion
			$debug && print STDERR "FOUND ".$gene->biotype.": ".$gene->stable_id." and AMB = $amb\n";			

			#set generic info
			$hash->{'trapmap_region'}{'region_id'} = $gene_id;
			$hash->{'trapmap_region'}{'overlap'} = $ti->{'coverage'};
			$hash->{'trapmap_region'}{'type'} = $gene->biotype;
			$hash->{'trapmap_region'}{'annotation_ambiguous'} = $amb;				
			$hash->{'trapmap_region'}{'number_trapblocks'} = scalar @{$trapblock};
			my $number_annotate_trapblocks = 0;

			my @regions = @{$self->region->rearrange_predict_exons($gene,$gene_id)};
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

			my $update_trapcheck = qq{UPDATE trapcheck SET annotated = 1, checked = 1, mapped = 1 WHERE trap_id = $trap_id};
			$self->load->fetch->update($update_trapcheck);

			#load hash

			$hash->{'additional'}{'processing'} = "annotate with genescan prediction";
			$hash->{'additional'}{'comment'} = "";
			$hash->{'additional'}{'label'} = $self->region->get_label($gene);
			$hash->{'additional'}{'trap_id'} = $trap_id;
			$hash->{'additional'}{'user'} = "mysql_dev";
			$hash->{'additional'}{'note'} = 'GENESCAN';
			$self->load->load_trapadditional($hash->{'additional'});

			#$hash->{'trapmap_region'}{'trapblock_id'} = $trapblock_id;
			$hash->{'trapmap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
			$self->load->load_trapmap_region($hash->{'trapmap_region'});

		}#end intersect gene
		print STDERR "------ got here 448\n";
	}#end gene loop

	if ($yes_gene) {
		my $timestmp = `date`;
		print STDERR "------ DONE annotate gene with annotate_with_ensembl_genescan\n";
		return 1;
	}else{
		print STDERR "------ no genescan prediction found\n";
		#print STDERR Dumper $trap;
		my $region_id = $self->region->slice($slice);

		$hash->{'trapmap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapmap_region'}{'type'} = "genomic";
		$hash->{'trapmap_region'}{'overlap'} = 100;
		$hash->{'trapmap_region'}{'region_id'} = $region_id;

		$self->load->load_trapmap_region($hash->{'trapmap_region'});
		my $update = qq{UPDATE trapcheck SET annotated = 0, checked = 1 WHERE trap_id = $trap_id};
		$self->load->fetch->update($update);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_ensembl_genescan\n";
		return 0;
	}
	
	
	
}

sub annotate_with_unigene {
	my ($self, $trap,$tmf, $slice) = @_;
	my $unigene_coverage = $conf{'annotation'}{'unigene_coverage'};
	my $unigene_perc_id = $conf{'annotation'}{'unigene_perc_id'};
	my $trap_id = $trap->{'trap_id'};
	my $type = "unigene";
	my $moltype = $trap->{'mol_type'};
	my $seqdir = $trap->{'sequencing_direction'};
	my $trapname = $trap->{'trap_name'};
	my $hash;
	my $yes_gene = 0;
	my $trapmap_id = $tmf->display_name;
	my $trapblock = $self->load->fetch->get_selected_trapblock_by_trapmap_id($trapmap_id);	
	
	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		die;
	}
	
	$hash->{'trapmap_region'}{'trapmap_id'} = $trapmap_id;
	
	my @dna_dna_align_feats = @{$slice->get_all_DnaAlignFeatures ('Unigene')};
	my @dnadnafeat;
	my $ori_start = $tmf->start;
	my $ori_end = $tmf->end;
	my $s = 1;
	my $e = ($tmf->end - $tmf->start) + 1;
	$tmf->start($s);
	$tmf->end($e);
	
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
			if ($dnadna->hstrand ne $tmf->strand){$amb = 1}
			
			my @subset = $col->features_in_range (-range=>$tmf, -contain=> 0);	
			my $coverage;
			foreach my $s (@subset) {	
				my $inter = $tmf->intersection($s);	
				my $cov = ($inter->length/$s->length)*100;
				unless (defined $coverage){$coverage = $cov};
				if ($cov >= $coverage) { #first time always true
					my $ghash;
					$ghash->{'seq_id'} = $s->display_name;
					my $start = $ori_start - $s->start;
					if ($tmf->start < $s->start){$start = $ori_start + $s->start}
					my $end = $ori_start - $s->end;
					if ($tmf->end < $s->end){$end = $ori_end + $s->end}
					$ghash->{'start'} = $start;
					$ghash->{'end'} = $end;
					$ghash->{'name'} = $tmf->display_name;
					$ghash->{'strand'} = $dnadna->hstrand;
					$ghash->{'parent_id'} = 0;
					$ghash->{'rank'} = 0;
					$ghash->{'description'} = 'Unigene';
					
					my $gene_id = $self->load->load_region($ghash);
					
					$hash->{'trapmap_region'}{'region_id'} = $gene_id;
					$hash->{'trapmap_region'}{'overlap'} = $coverage;
					$hash->{'trapmap_region'}{'type'} = 'Unigene';
					$hash->{'trapmap_region'}{'annotation_ambiguous'} = $amb;				
					$hash->{'trapmap_region'}{'number_trapblocks'} = scalar @{$trapblock};

					my $update_trapcheck = qq{UPDATE trapcheck SET annotated = 1, checked = 1, mapped = 1 WHERE trap_id = $trap_id};
					$self->load->fetch->update($update_trapcheck);
				
					#load hash
				
					$hash->{'additional'}{'processing'} = "annotate with Unigene DnaAlignFeatures";
					$hash->{'additional'}{'comment'} = "";
					$hash->{'additional'}{'label'} = "";
					$hash->{'additional'}{'trap_id'} = $trap_id;
					$hash->{'additional'}{'user'} = "mysql_dev";
					$hash->{'additional'}{'note'} = 'UNIGENE';
					$self->load->load_trapadditional($hash->{'additional'});
				
					#$hash->{'trapmap_region'}{'trapblock_id'} = $trapblock_id;
					$self->load->load_trapmap_region($hash->{'trapmap_region'});
					$yes_gene ++;
				}
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
		#print STDERR Dumper $trap;
		my $region_id = $self->region->slice($slice);

		$hash->{'trapmap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapmap_region'}{'type'} = "genomic";
		$hash->{'trapmap_region'}{'overlap'} = 100;
		$hash->{'trapmap_region'}{'region_id'} = $region_id;

		$self->load->load_trapmap_region($hash->{'trapmap_region'});
		my $update = qq{UPDATE trapcheck SET annotated = 0, checked = 1 WHERE trap_id = $trap_id};
		$self->load->fetch->update($update);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_unigene\n";
		return 0;
	}
}

sub annotate_with_ensest {
	my ($self, $trap,$tmf, $slice) = @_;
	my $est_coverage = $conf{'annotation'}{'est_coverage'};
	my $est_perc_id = $conf{'annotation'}{'est_perc_id'};
	my $trap_id = $trap->{'trap_id'};
	my $type = "est";
	my $moltype = $trap->{'mol_type'};
	my $seqdir = $trap->{'sequencing_direction'};
	my $trapname = $trap->{'trap_name'};
	my $hash;
	my $yes_gene = 0;
	my $trapmap_id = $tmf->display_name;
	my $trapblock = $self->load->fetch->get_selected_trapblock_by_trapmap_id($trapmap_id);	
	
	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		die;
	}
	
	$hash->{'trapmap_region'}{'trapmap_id'} = $trapmap_id;
	
	my @dna_dna_align_feats = @{$slice->get_all_DnaAlignFeatures ('mouse_est')};
	my @dnadnafeat;
	my $ori_start = $tmf->start;
	my $ori_end = $tmf->end;
	my $s = 1;
	my $e = ($tmf->end - $tmf->start) + 1;
	$tmf->start($s);
	$tmf->end($e);
	
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
			if ($dnadna->hstrand ne $tmf->strand){$amb = 1}
			
			my @subset = $col->features_in_range (-range=>$tmf, -contain=> 0);	
			my $coverage;
			foreach my $s (@subset) {	
				my $inter = $tmf->intersection($s);	
				my $cov = ($inter->length/$s->length)*100;
				unless (defined $coverage){$coverage = $cov};
				if ($cov >= $coverage) { #first time always true
					my $ghash;
					$ghash->{'seq_id'} = $s->display_name;
					my $start = $ori_start - $s->start;
					if ($tmf->start < $s->start){$start = $ori_start + $s->start}
					my $end = $ori_start - $s->end;
					if ($tmf->end < $s->end){$end = $ori_end + $s->end}
					$ghash->{'start'} = $start;
					$ghash->{'end'} = $end;
					$ghash->{'name'} = $tmf->display_name;
					$ghash->{'strand'} = $dnadna->hstrand;
					$ghash->{'parent_id'} = 0;
					$ghash->{'rank'} = 0;
					$ghash->{'description'} = 'EST';
					
					my $gene_id = $self->load->load_region($ghash);
					
					$hash->{'trapmap_region'}{'region_id'} = $gene_id;
					$hash->{'trapmap_region'}{'overlap'} = $coverage;
					$hash->{'trapmap_region'}{'type'} = 'EST';
					$hash->{'trapmap_region'}{'annotation_ambiguous'} = $amb;				
					$hash->{'trapmap_region'}{'number_trapblocks'} = scalar @{$trapblock};

					my $update_trapcheck = qq{UPDATE trapcheck SET annotated = 1, checked = 1, mapped = 1 WHERE trap_id = $trap_id};
					$self->load->fetch->update($update_trapcheck);
				
					#load hash
				
					$hash->{'additional'}{'processing'} = "annotate with Mouse Est DnaAlignFeatures";
					$hash->{'additional'}{'comment'} = "";
					$hash->{'additional'}{'label'} = "";
					$hash->{'additional'}{'trap_id'} = $trap_id;
					$hash->{'additional'}{'user'} = "mysql_dev";
					$hash->{'additional'}{'note'} = 'EST';
					$self->load->load_trapadditional($hash->{'additional'});
				
					#$hash->{'trapmap_region'}{'trapblock_id'} = $trapblock_id;
					$self->load->load_trapmap_region($hash->{'trapmap_region'});
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
		#print STDERR Dumper $trap;
		my $region_id = $self->region->slice($slice);

		$hash->{'trapmap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapmap_region'}{'type'} = "genomic";
		$hash->{'trapmap_region'}{'overlap'} = 100;
		$hash->{'trapmap_region'}{'region_id'} = $region_id;

		$self->load->load_trapmap_region($hash->{'trapmap_region'});
		my $update = qq{UPDATE trapcheck SET annotated = 0, checked = 1 WHERE trap_id = $trap_id};
		$self->load->fetch->update($update);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_ensest\n";
		return 0;
	}
	
}


sub annotate_with_ensrepeats {
	my ($self, $trap,$tmf, $slice) = @_;
	my $est_coverage = $conf{'annotation'}{'est_coverage'};
	my $est_perc_id = $conf{'annotation'}{'est_perc_id'};
	my $trap_id = $trap->{'trap_id'};
	my $type = "est";
	my $moltype = $trap->{'mol_type'};
	my $seqdir = $trap->{'sequencing_direction'};
	my $trapname = $trap->{'trap_name'};
	my $hash;
	my $yes_gene = 0;
	my $trapmap_id = $tmf->display_name;
	my $trapblock = $self->load->fetch->get_selected_trapblock_by_trapmap_id($trapmap_id);	
	my $ori_start = $tmf->start;
	my $ori_end = $tmf->end;
	my $s = 1;
	my $e = ($tmf->end - $tmf->start) + 1;
	$tmf->start($s);
	$tmf->end($e);

	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		die;
	}
	
	$hash->{'trapmap_region'}{'trapmap_id'} = $trapmap_id;
	
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
			if ($repeat->hstrand ne $tmf->strand){$amb = 1}
			my $col = new Bio::SeqFeature::Collection();
			my $totaladded = $col->add_features(\@repeatfeat);
			
			foreach my $f (@repeatfeat) {
				$debug && print STDOUT "Feature: ".$f->display_name."\t".$f->start."\t".$f->end."\n";
				
				my @subset = $col->features_in_range(-range=>$f, -contain=> 0);
				my $coverage;
				foreach my $s (@subset) {
					$debug && print STDOUT "In range: " . $s->display_name . "\t". $s->start . "\t". $s->end . "\n";
					my $inter = $tmf->intersection($s);	
					my $cov = ($inter->length/$s->length)*100;
					unless (defined $coverage){$coverage = $cov};
					if ($cov >= $coverage) { #first time always true
						my $ghash;
						$ghash->{'seq_id'} = $s->display_name;
						my $start = $ori_start - $s->start;
						if ($tmf->start < $s->start){$start = $ori_start + $s->start}
						my $end = $ori_start - $s->end;
						if ($tmf->end < $s->end){$end = $ori_end + $s->end}
						$ghash->{'start'} = $start;
						$ghash->{'end'} = $end;
						$ghash->{'name'} = $tmf->display_name;
						$ghash->{'strand'} = $repeat->hstrand;
						$ghash->{'parent_id'} = 0;
						$ghash->{'rank'} = 0;
						$ghash->{'description'} = 'Repeats';
						
						my $gene_id = $self->load->load_region($ghash);
						
						$hash->{'trapmap_region'}{'region_id'} = $gene_id;
						$hash->{'trapmap_region'}{'overlap'} = $coverage;
						$hash->{'trapmap_region'}{'type'} = 'Repeats';
						$hash->{'trapmap_region'}{'annotation_ambiguous'} = $amb;				
						$hash->{'trapmap_region'}{'number_trapblocks'} = scalar @{$trapblock};
	
						my $update_trapcheck = qq{UPDATE trapcheck SET annotated = 1, checked = 1, mapped = 1 WHERE trap_id = $trap_id};
						$self->load->fetch->update($update_trapcheck);
					
						#load hash
					
						$hash->{'additional'}{'processing'} = "annotate with Repeat Features";
						$hash->{'additional'}{'comment'} = "";
						$hash->{'additional'}{'label'} = "";
						$hash->{'additional'}{'trap_id'} = $trap_id;
						$hash->{'additional'}{'user'} = "mysql_dev";
						$hash->{'additional'}{'note'} = 'Repeats';
						$self->load->load_trapadditional($hash->{'additional'});
					
						#$hash->{'trapmap_region'}{'trapblock_id'} = $trapblock_id;
						$self->load->load_trapmap_region($hash->{'trapmap_region'});
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
		#print STDERR Dumper $trap;
		my $region_id = $self->region->slice($slice);

		$hash->{'trapmap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapmap_region'}{'type'} = "genomic";
		$hash->{'trapmap_region'}{'overlap'} = 100;
		$hash->{'trapmap_region'}{'region_id'} = $region_id;

		$self->load->load_trapmap_region($hash->{'trapmap_region'});
		my $update = qq{UPDATE trapcheck SET annotated = 0, checked = 1 WHERE trap_id = $trap_id};
		$self->load->fetch->update($update);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_ensrepeats\n";
		return 0;
	}
}	
	

sub annotate_with_mouse_cdna {
	my ($self, $trap,$tmf, $slice) = @_;
	my $cdna_coverage = $conf{'annotation'}{'cdna_coverage'};
	my $cdna_perc_id = $conf{'annotation'}{'cdna_perc_id'};
	my $trap_id = $trap->{'trap_id'};
	my $type = "cdna";
	my $moltype = $trap->{'mol_type'};
	my $seqdir = $trap->{'sequencing_direction'};
	my $trapname = $trap->{'trap_name'};
	my $hash;
	my $yes_gene = 0;
	my $trapmap_id = $tmf->display_name;
	my $trapblock = $self->load->fetch->get_selected_trapblock_by_trapmap_id($trapmap_id);	
	
	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		die;
	}
	
	$hash->{'trapmap_region'}{'trapmap_id'} = $trapmap_id;
	
	my @dna_dna_align_feats = @{$slice->get_all_DnaAlignFeatures ('mouse_cdna')};
	my @dnadnafeat;
	my $ori_start = $tmf->start;
	my $ori_end = $tmf->end;
	my $s = 1;
	my $e = ($tmf->end - $tmf->start) + 1;
	$tmf->start($s);
	$tmf->end($e);
	
	foreach my $dnadna (@dna_dna_align_feats) {
		
		if ($dnadna->percent_id >= $cdna_perc_id) {
			
			my $feature = Bio::SeqFeature::Generic->new (	-display_name => $dnadna->hseqname,
									-start => $dnadna->hstart,
									-end => $dnadna->hend,
									-strand => $dnadna->hstrand);
			push @dnadnafeat ,$feature;
			
			my $col = new Bio::SeqFeature::Collection();
			my $totaladded = $col->add_features(\@dnadnafeat);
			my $amb = 0;
			if ($dnadna->hstrand ne $tmf->strand){$amb = 1}
			
			my @subset = $col->features_in_range (-range=>$tmf, -contain=> 0);	
			my $coverage;
			foreach my $s (@subset) {	
				my $inter = $tmf->intersection($s);	
				my $cov = ($inter->length/$s->length)*100;
				unless (defined $coverage){$coverage = $cov};
				if ($cov >= $coverage) { #first time always true
					my $ghash;
					$ghash->{'seq_id'} = $s->display_name;
					my $start = $ori_start - $s->start;
					if ($tmf->start < $s->start){$start = $ori_start + $s->start}
					my $end = $ori_start - $s->end;
					if ($tmf->end < $s->end){$end = $ori_end + $s->end}
					$ghash->{'start'} = $start;
					$ghash->{'end'} = $end;
					$ghash->{'name'} = $tmf->display_name;
					$ghash->{'strand'} = $dnadna->hstrand;
					$ghash->{'parent_id'} = 0;
					$ghash->{'rank'} = 0;
					$ghash->{'description'} = 'CDNA';
					
					my $gene_id = $self->load->load_region($ghash);
					
					$hash->{'trapmap_region'}{'region_id'} = $gene_id;
					$hash->{'trapmap_region'}{'overlap'} = $coverage;
					$hash->{'trapmap_region'}{'type'} = 'CDNA';
					$hash->{'trapmap_region'}{'annotation_ambiguous'} = $amb;				
					$hash->{'trapmap_region'}{'number_trapblocks'} = scalar @{$trapblock};

					my $update_trapcheck = qq{UPDATE trapcheck SET annotated = 1, checked = 1, mapped = 1 WHERE trap_id = $trap_id};
					$self->load->fetch->update($update_trapcheck);
				
					#load hash
				
					$hash->{'additional'}{'processing'} = "annotate with Mouse CDNA DnaAlignFeatures";
					$hash->{'additional'}{'comment'} = "";
					$hash->{'additional'}{'label'} = "";
					$hash->{'additional'}{'trap_id'} = $trap_id;
					$hash->{'additional'}{'user'} = "mysql_dev";
					$hash->{'additional'}{'note'} = 'CDNA';
					$self->load->load_trapadditional($hash->{'additional'});
				
					#$hash->{'trapmap_region'}{'trapblock_id'} = $trapblock_id;
					$self->load->load_trapmap_region($hash->{'trapmap_region'});
					$yes_gene ++;
				}
			}
		}
	}
	
	
	if ($yes_gene) {
		my $timestmp = `date`;
		print STDERR "------ DONE annotate gene with annotate_with_mouse_cdna\n";
		return 1;
	}else{
		print STDERR "------ no cdna found\n";
		#print STDERR Dumper $trap;
		my $region_id = $self->region->slice($slice);

		$hash->{'trapmap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapmap_region'}{'type'} = "genomic";
		$hash->{'trapmap_region'}{'overlap'} = 100;
		$hash->{'trapmap_region'}{'region_id'} = $region_id;

		$self->load->load_trapmap_region($hash->{'trapmap_region'});
		my $update = qq{UPDATE trapcheck SET annotated = 0, checked = 1 WHERE trap_id = $trap_id};
		$self->load->fetch->update($update);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_mouse_cdna\n";
		return 0;
	}
	
	
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
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
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
			my $gene_id = $self->region->gene($gene); #load gene in region table
			my $amb = $self->trap_utility->check_strand($trap,$tmf,$gene); #used to check if the alignment is on the correct strand of the gene
			
			if ($debug && $amb == 1){next;}
			#if the strand is not correct it isn't possible to infer the insertion
			$debug && print STDERR "------ FOUND ".$gene->biotype.": ".$gene->stable_id." and AMB = $amb\n";			

			#set generic info
			$hash->{'trapmap_region'}{'region_id'} = $gene_id;
			$hash->{'trapmap_region'}{'overlap'} = $ti->{'coverage'};
			$hash->{'trapmap_region'}{'type'} = $gene->biotype;
			$hash->{'trapmap_region'}{'annotation_ambiguous'} = $amb;				
			$hash->{'annotation'}{'trapmap_id'} = $tmf->display_name;
			$hash->{'trapmap_region'}{'number_trapblocks'} = scalar @{$trapblock};
			my $number_annotate_trapblocks = 0;

			my @regions = @{$self->region->rearrange_exons($gene,$gene_id)};
			my $first = $regions[0]; 
			my $last = $regions[$#regions];
			my $first_exon_id = $first->{'region_id'};
			my $last_exon_id = $last->{'region_id'};
			my $first_start = $first->{'region_start'}; 
			my $last_end = $last->{'region_end'};
			$hash->{'trapmap_region'}{'total_number_exons'} = 1+(int(scalar @regions)/2);
			my ($nti) = $self->trap_utility->intersection($regions[0],$first_start,$last_end,$tmf);			
			unless ($nti->{'i'}) {
				print STDERR "------ no coding block found for gene: $gene_id\n"; 
				$hash->{'trapmap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
				$self->load->load_trapmap_region($hash->{'trapmap_region'});
				next;
			}

			#if the gene is a non coding gene and trap is RACE it is treated like a genomic
			if ($gene->biotype !~ /coding/ && $moltype eq 'mRNA' ){
				print STDERR "------ Got here because $moltype and ".$gene->biotype."\n";
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
					print STDERR warn "------ got a problem with TRAPMAP $trapmap_id: found no region intersecation block for gene ".$gene->stable_id."\n";
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
			$hash = $self->trap_utility->define_internal_trapped_region($hash,$gene,$trapped{$trapped_region_rank}{'tbf'},$first_trapped,$last_trapped,$pre_region,$post_region,$seqdir,$moltype);

			$debug && print STDERR "------ no $trapped_region_rank\n";
			$trapblock_id = $trapped{$trapped_region_rank}{'trapblock_id'};
			$debug && print STDERR "------ trapblock_id $trapblock_id\n";

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
			$hash->{'additional'}{'user'} = "mysql_dev";
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
		print STDERR "------ DONE annotate gene with annotate_with_ensembl_gene\n";
		return 1;
	}else{
		print STDERR "------ no gene found\n";
		#print STDERR Dumper $trap;
		my $region_id = $self->region->slice($slice);

		$hash->{'trapmap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapmap_region'}{'type'} = "genomic";
		$hash->{'trapmap_region'}{'overlap'} = 100;
		$hash->{'trapmap_region'}{'region_id'} = $region_id;

		$self->load->load_trapmap_region($hash->{'trapmap_region'});
		my $update = qq{UPDATE trapcheck SET annotated = 0, checked = 1 WHERE trap_id = $trap_id};
		$self->load->fetch->update($update);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_ensembl_gene\n";
		return 0;
	}
}


=head2 annotate_with_ensembl_estgene

	Title    : annotate_with_ensembl_estgene
	Usage    : $obj->annotate_with_ensembl_gene(Args)
 	Function : annotate trap wth ensembl gene
 	Returns  : 
 	Args     : trap href, hit feat, slice on genome	  

=cut

sub annotate_with_ensembl_estgene{
	my ($self, $trap,$tmf, $slice) = @_;
	my $trap_id = $trap->{'trap_id'};
	my $type = "ensembl est gene";
	my $moltype = $trap->{'mol_type'};
	my $seqdir = $trap->{'sequencing_direction'};
	my $trapname = $trap->{'trap_name'};
	my $hash;
	my $yes_gene = 0;
	my $trapmap_id = $tmf->display_name;
	my $trapblock = $self->load->fetch->get_selected_trapblock_by_trapmap_id($trapmap_id);
	unless (scalar @{$trapblock}){
		print STDERR "------ no blocks were accepted for this map\n";#shouldn't occurr!!!
		die;
	}
	$hash->{'trapmap_region'}{'trapmap_id'} = $trapmap_id;
	print STDERR "------ got here 706\n";
	foreach my $gene (@{$slice->get_all_Genes}) {
		my $gstart = $gene->seq_region_start;
		my $gend = $gene->seq_region_end;
		my $gstrand = $gene->strand;

	    #test if the feature (trapmap) inteserct a gene
		my ($ti) = $self->trap_utility->intersection($gene,$gstart,$gend,$tmf);
		print STDERR "------ got here 714\n";
		if ($ti->{'i'}) {
			$yes_gene ++;
			
			my $gene_id = $self->region->gene($gene); #load gene in region table
			print STDERR "------ got here 719\n";
			my $amb = $self->trap_utility->check_strand($trap,$tmf,$gene); #used to check if the alignment is on the correct strand of the gene
			#if ($debug && $amb == 1){next;}
			
			#if the strand is not correct it isn't possible to infer the insertion
			$debug && print STDERR "FOUND ".$gene->biotype.": ".$gene->stable_id." and AMB = $amb\n";			

			#set generic info
			$hash->{'trapmap_region'}{'region_id'} = $gene_id;
			$hash->{'trapmap_region'}{'overlap'} = $ti->{'coverage'};
			$hash->{'trapmap_region'}{'type'} = $gene->biotype;
			$hash->{'trapmap_region'}{'annotation_ambiguous'} = $amb;				
			$hash->{'trapmap_region'}{'number_trapblocks'} = scalar @{$trapblock};
			my $number_annotate_trapblocks = 0;

			my @regions = @{$self->region->rearrange_est_exons($gene,$gene_id)};
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

			my $update_trapcheck = qq{UPDATE trapcheck SET annotated = 1, checked = 1, mapped = 1 WHERE trap_id = $trap_id};
			$self->load->fetch->update($update_trapcheck);

			#load hash

			$hash->{'additional'}{'processing'} = "annotate with ensembl est gene";
			$hash->{'additional'}{'comment'} = "";
			$hash->{'additional'}{'label'} = $self->region->get_label($gene);
			$hash->{'additional'}{'trap_id'} = $trap_id;
			$hash->{'additional'}{'user'} = "mysql_dev";
			#$hash->{'additional'}{'note'} = $posrelative->{'label'}.": ".$posrelative->{'comment'};
			$self->load->load_trapadditional($hash->{'additional'});

			#$hash->{'trapmap_region'}{'trapblock_id'} = $trapblock_id;
			$hash->{'trapmap_region'}{'number_annotate_trapblocks'} = $number_annotate_trapblocks;
			$self->load->load_trapmap_region($hash->{'trapmap_region'});

		}#end intersect gene
		print STDERR "------ got here 764\n";
	}#end gene loop

	if ($yes_gene) {
		my $timestmp = `date`;
		print STDERR "------ DONE annotate gene with annotate_with_ensembl_estgene\n";
		return 1;
	}else{
		print STDERR "------ no estgene found\n";
		#print STDERR Dumper $trap;
		my $region_id = $self->region->slice($slice);

		$hash->{'trapmap_region'}{'annotation_ambiguous'} = 1;
		$hash->{'trapmap_region'}{'type'} = "genomic";
		$hash->{'trapmap_region'}{'overlap'} = 100;
		$hash->{'trapmap_region'}{'region_id'} = $region_id;

		$self->load->load_trapmap_region($hash->{'trapmap_region'});
		my $update = qq{UPDATE trapcheck SET annotated = 0, checked = 1 WHERE trap_id = $trap_id};
		$self->load->fetch->update($update);
		my $timestmp = `date`;
		print STDERR "------ DONE annotate genomic region with annotate_with_ensembl_estgene\n";
		return 0;
	}
}


sub annotate {
	my ($self,$trap_id, $idcount) = @_;
	
	my $timestmp = `date`;
	print STDERR "---- START annotate for TRAP ID $trap_id \n";
	
	my $trap = $self->load->fetch->get_trap_by_id($trap_id);
	my @trapmap = @{$self->load->fetch->get_trapmap_by_trap_id($trap_id)};
	foreach my $rhref (@trapmap) {
		
		my $region;
		if ($rhref->{'hit_id'} =~ /NT/) {$region = "supercontig";} 
		else {$region = "chromosome";}
		$rhref->{'hit_id'} =~ s/chr//;
		
		my $tmf = Bio::SeqFeature::Generic->new(	-display_name => $rhref->{'trapmap_id'},
											-start => $rhref->{'start'},
											-end => $rhref->{'end'},
											-strand => $rhref->{'strand'}	);
		
		my $hash;
		
		my $trapblock = $self->load->fetch->get_selected_trapblock_by_trapmap_id($rhref->{'trapmap_id'});
		unless (scalar @{$trapblock}){
			print STDERR "ERROR: no blocks were accepted for this map. Exit program\n";#shouldn't occurr!!!
			die;
		}
		$hash->{'trapmap_region'}{'trapmap_id'} = $rhref->{'trapmap_id'};
		$hash->{'tb'} = $trapblock;
		$hash->{'trap'} = $trap;
		$hash->{'tmf'} = $tmf;
		my $slice_core_adaptor = $self->slicecoreadpt;
			my $slice_core = $slice_core_adaptor->fetch_by_region($region,$rhref->{'hit_id'}, $rhref->{'start'}, $rhref->{'end'});
		my $found = 0;
		if ($self->do_ensgene) {
			
			if (defined $slice_core){
				$found = $self->annotate_with_ensembl_gene($trap,$tmf, $slice_core);
			}
		}
		if ($found) {
			my $timestmp = `date`;
			print STDERR "---- DONE annotate gene with annotate_with_ensembl_gene\n";
			#return $found;
		}
		print STDERR "---- FOUND $found try to do EST ".$self->do_ensestgene."\n";
		if ($self->do_ensestgene){
			my $slice_est_adaptor = $self->sliceestadpt;
			my $slice_est = $slice_est_adaptor->fetch_by_region($region,$rhref->{'hit_id'}, $rhref->{'start'}, $rhref->{'end'});
			print STDERR "---- got here 838\n";
			if (defined $slice_est){
				print STDERR "---- SLICE EST defined\n";
				$found = $self->annotate_with_ensembl_estgene($trap,$tmf, $slice_est);
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
		
		print STDERR "---- FOUND $found try to do GENESCAN ".$self->do_genescan."\n";

 		if ($self->do_genescan) {
			print STDERR "annotate_with_ensembl_genescan\n";
			$found = $self->annotate_with_ensembl_genescan($trap,$tmf, $slice_core);
		}
		
		print STDERR "---- FOUND $found try to do UNIGENE ".$self->do_unigene."\n";
		if ($self->do_unigene) {
			print STDERR "annotate_with_unigene\n";
			$found = $self->annotate_with_unigene($trap,$tmf, $slice_core);
			
		}
		print STDERR "---- FOUND $found try to do CDNA ".$self->do_mouse_cdna."\n";
		if ($self->do_mouse_cdna) {
			print STDERR "annotate_with_cdna\n";
			$found = $self->annotate_with_mouse_cdna($trap,$tmf, $slice_core);
		}
		print STDERR "---- FOUND $found try to do ENSEST ".$self->do_ensest."\n";
		if ($self->do_ensest) {
			print STDERR "annotate_ensest\n";
			$found = $self->annotate_with_ensest($trap,$tmf, $slice_core);
		}
		print STDERR "---- FOUND $found try to do REPEATS ".$self->do_ensrepeats."\n";
		if ($self->do_ensrepeats) {
			print STDERR "annotate_repeats\n";
			$found = $self->annotate_with_ensrepeats($trap,$tmf, $slice_core);
		}
		print STDERR "---- FOUND $found. END\n";
	}
	
	$timestmp = `date`;
	print STDERR "---- DONE annotate\n";
	return 1;
}

1;
