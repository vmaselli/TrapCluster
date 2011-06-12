#!/usr/bin/perl
package Bio::MCE::Pipeline::Coding;

use strict;
use warnings;
use Data::Dumper;
use vars qw(@ISA);
use vars qw($CONF);

use Bio::EnsEMBL::Registry;
use Bio::SeqFeature::Collection;
use Bio::SeqFeature::Generic;
use Bio::Range;


=head2 test_coding

 Title   : test_coding

 Usage   : my($coding,$length) = Bio::MCE::Pipeline->test_coding($specie,$coord_sys_name,$chr,$start,$end);

 Function: give the information if the slice coming from supplied coordinates overlap transcribed
           regions

 Returns : number of coding features that overlap with slice and total length of overlap

 Args    : -1 the specie
           -2 the coordinate system name
           -3 the region name
           -4 the start
           -5 the end

 NOTE    : 

=cut

sub test_coding {

	my $class = shift;
	my $specie = shift;
	my $coord_sys_name = shift;
  my $chr = shift;
  my $start = shift;
  my $end = shift;
  my $elem_length = $end-$start+1;
	my $chr_length = Bio::MCE::Pipeline::Coding->get_chr_length($specie,$coord_sys_name,$chr);

	# slice to take in order to avoid the fact to take only the super feature
  # because if only the superfeature overlap our range it will be taken and
  # not dropped because it is one....
	my $slice_start = $start - 50000;
	my $slice_end = $end + 50000;
	$slice_start = 1 if $slice_start < 1;
	$slice_end = $chr_length -1 if $slice_end>=$chr_length;

	my $slice_adaptor_core = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice");
	my $slice_adaptor_other = Bio::EnsEMBL::Registry->get_adaptor($specie,"otherfeatures","slice");

  my $slice_core = $slice_adaptor_core->fetch_by_region($coord_sys_name,$chr,$slice_start,$slice_end);
  my $slice_other = $slice_adaptor_other->fetch_by_region($coord_sys_name,$chr,$slice_start,$slice_end);

	# new coordinate referred to the slice of our region...
	my $new_start = $start-$slice_start+1;
	my $new_end = $new_start+$elem_length-1;

	# costruire la collection con tutti i EST e Prediction
	# ma bisogna prima convertirli in seqfeature generic	
	my $aln_f_c = $slice_core->get_all_DnaAlignFeatures('','75');
	my $aln_pf_c = $slice_core->get_all_ProteinAlignFeatures;
	my $prd_f_c = $slice_core->get_all_PredictionTranscripts;

	my $aln_f_o = $slice_other->get_all_DnaAlignFeatures('','75');
	my $aln_pf_o = $slice_other->get_all_ProteinAlignFeatures;
	my $prd_f_o = $slice_other->get_all_PredictionTranscripts;

	my $aln_g_c = _make_generic_aln($aln_f_c);
	my $aln_pg_c = _make_generic_aln($aln_pf_c);
	my $prd_g_c = _make_generic_prd($prd_f_c);

	my $aln_g_o = _make_generic_aln($aln_f_o);
	my $aln_pg_o = _make_generic_aln($aln_pf_o);
	my $prd_g_o = _make_generic_prd($prd_f_o);

	my $coll = Bio::SeqFeature::Collection->new;

	$coll->add_features($aln_g_c) if $aln_g_c;
	$coll->add_features($aln_pg_c) if $aln_pg_c;
	$coll->add_features($prd_g_c) if $prd_g_c;

	$coll->add_features($aln_g_o) if $aln_g_o;
	$coll->add_features($aln_pg_o) if $aln_pg_o;
	$coll->add_features($prd_g_o) if $prd_g_o;

	if($specie eq 'Mus musculus' or $specie eq 'Homo sapiens') {
		
		my $slice_adaptor_cdna = Bio::EnsEMBL::Registry->get_adaptor($specie,"cdna","slice");
		
		my $slice_cdna = $slice_adaptor_cdna->fetch_by_region($coord_sys_name,$chr,$slice_start,$slice_end);
		
		# costruire la collection con tutti i EST e Prediction
		# ma bisogna prima convertirli in seqfeature generic  
		
		my $aln_f_cdna = $slice_cdna->get_all_DnaAlignFeatures('','75');
		my $aln_pf_cdna = $slice_cdna->get_all_ProteinAlignFeatures;
		my $prd_f_cdna = $slice_cdna->get_all_PredictionTranscripts;
		
		my $aln_g_cdna = _make_generic_aln($aln_f_cdna);
		my $aln_pg_cdna = _make_generic_aln($aln_pf_cdna);
		my $prd_g_cdna = _make_generic_prd($prd_f_cdna);
		
		$coll->add_features($aln_g_cdna) if $aln_g_cdna;
		$coll->add_features($aln_pg_cdna) if $aln_pg_cdna;
		$coll->add_features($prd_g_cdna) if $prd_g_cdna;
	}

	my($coding,$overlap) = _test_coding($new_start,$new_end,$coll);

	return($coding,$overlap);
}


=head2 test_coding_feature_hashref

 Title   : test_coding_feature_hashref

 Usage   : my($coding_href) = Bio::MCE::Pipeline::Coding->test_coding_feature_href($specie,$coord_sys_name,$chr,$start,$end,$slce_strand,$features);

 Function: give the information about the coding overlap of feature or feature_pair coming
           from the valis feature or feature_pair table.
           YOU DON'T HAVE TO SUPPLY CHROMOSOMIC COORDS JUST THAT ONE IN VALIS
           RELATIVE TO THE REFERENCE REGION UNDER CONSIDERATION

 Returns : an hashref of hashref

 Args    : -1 the specie
           -2 the coordinate system name
           -3 the region name
           -4 the start
           -5 the end
           -6 slice_strand           
           -7 feature_href

 NOTE    : YOU DON'T HAVE TO SUPPLY CHROMOSOMIC COORDS JUST THAT ONE IN VALIS
           RELATIVE TO THE REFERENCE REGION UNDER CONSIDERATION


=cut

sub test_coding_feature_hashref {

	my $class = shift;
	my $specie = shift;
	my $coord_sys_name = shift;
  my $chr = shift;
  my $chr_start = shift;
  my $chr_end = shift;
	my $slice_strand = shift;
	my $feature_href = shift;
  #my $elem_length = $chr_end-$chr_start+1;
	my $chr_length = Bio::MCE::Pipeline::Coding->get_chr_length($specie,$coord_sys_name,$chr);

	# slice to take in order to avoid the fact to take only the super feature
  # because if only the superfeature overlap our range it will be taken and
  # not dropped because it is one....
	my $slice_start = $chr_start - 50000;
	my $slice_end = $chr_end + 50000;
	$slice_start = 1 if $slice_start < 1;
	$slice_end = $chr_length - 1 if $slice_end>=$chr_length;

  my $slice_adaptor_core = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice");
  my $slice_adaptor_other = Bio::EnsEMBL::Registry->get_adaptor($specie,"otherfeatures","slice");

  my $slice_core = $slice_adaptor_core->fetch_by_region($coord_sys_name,$chr,$slice_start,$slice_end);
  my $slice_other = $slice_adaptor_other->fetch_by_region($coord_sys_name,$chr,$slice_start,$slice_end);

	# costruire la collection con tutti i EST e Prediction
	# ma bisogna prima convertirli in seqfeature generic	

  my $aln_f_c = $slice_core->get_all_DnaAlignFeatures('','75');
  my $aln_pf_c = $slice_core->get_all_ProteinAlignFeatures;
  my $prd_f_c = $slice_core->get_all_PredictionTranscripts;

  my $aln_f_o = $slice_other->get_all_DnaAlignFeatures('','75');
  my $aln_pf_o = $slice_other->get_all_ProteinAlignFeatures;
  my $prd_f_o = $slice_other->get_all_PredictionTranscripts;

  my $aln_g_c = _make_generic_aln($aln_f_c);
  my $aln_pg_c = _make_generic_aln($aln_pf_c);
  my $prd_g_c = _make_generic_prd($prd_f_c);

  my $aln_g_o = _make_generic_aln($aln_f_o);
  my $aln_pg_o = _make_generic_aln($aln_pf_o);
  my $prd_g_o = _make_generic_prd($prd_f_o);

  my $coll = Bio::SeqFeature::Collection->new;

  $coll->add_features($aln_g_c) if $aln_g_c;
  $coll->add_features($aln_pg_c) if $aln_pg_c;
  $coll->add_features($prd_g_c) if $prd_g_c;

  $coll->add_features($aln_g_o) if $aln_g_o;
  $coll->add_features($aln_pg_o) if $aln_pg_o;
  $coll->add_features($prd_g_o) if $prd_g_o;

	if($specie eq 'Mus musculus') {
		
		my $slice_adaptor_cdna = Bio::EnsEMBL::Registry->get_adaptor($specie,"cdna","slice");
		
		my $slice_cdna = $slice_adaptor_cdna->fetch_by_region($coord_sys_name,$chr,$slice_start,$slice_end);
		
		# costruire la collection con tutti i EST e Prediction
		# ma bisogna prima convertirli in seqfeature generic  
		
		my $aln_f_cdna = $slice_cdna->get_all_DnaAlignFeatures('','75');
		my $aln_pf_cdna = $slice_cdna->get_all_ProteinAlignFeatures;
		my $prd_f_cdna = $slice_cdna->get_all_PredictionTranscripts;
		
		my $aln_g_cdna = _make_generic_aln($aln_f_cdna);
		my $aln_pg_cdna = _make_generic_aln($aln_pf_cdna);
		my $prd_g_cdna = _make_generic_prd($prd_f_cdna);
		
		$coll->add_features($aln_g_cdna) if $aln_g_cdna;
		$coll->add_features($aln_pg_cdna) if $aln_pg_cdna;
		$coll->add_features($prd_g_cdna) if $prd_g_cdna;
	}

	my $coding_href;

	foreach my $f(@$feature_href) {

		my($rel_start,$rel_end,$new_start,$new_end,$feature_id);

		# new coordinate referred to the slice of our region...
		if($f->{'start'}) {
			($rel_start,$rel_end) = Bio::MCE::Pipeline::Postanalysis::Valis->feature_from_ref_to_chr($f,$chr_start,$chr_end,$slice_strand);
			$new_start = $rel_start-$slice_start+1;
			$new_end = $new_start+($f->{'end'}-$f->{'start'}+1)-1;
			$feature_id = $f->{'feature_id'};
		}
		elsif($f->{'q_start'}) {
			($rel_start,$rel_end) = Bio::MCE::Pipeline::Postanalysis::Valis->q_feature_from_ref_to_chr($f,$chr_start,$chr_end,$slice_strand);
			$new_start = $rel_start-$slice_start+1;
			$new_end = $new_start+($f->{'q_end'}-$f->{'q_start'}+1)-1;
			$feature_id = $f->{'feature_pair_id'};
		}

		($coding_href->{$feature_id}->{'coding'},$coding_href->{$feature_id}->{'overlap'}) = _test_coding($new_start,$new_end,$coll);
	}
	return($coding_href);
}


=head2 _make_generic_aln

 Title   : _make_generic_aln

 Usage   : my $alnf = _make_generic_aln($dna_align_feature);

 Function: 

 Returns : 

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub _make_generic_aln{

  my $aln_f = shift;
  my $aln_g;
  my $aln_ok;

  foreach my $f(@$aln_f) {		

		next unless $f->length > 2;

		my $hlength = $f->hend-$f->hstart+1;
		my $slength = $f->end-$f->start+1;

		next if $slength>(2*$hlength);

    my $g = Bio::SeqFeature::Generic->new(-display_name => $f->display_id,
                                          -start => $f->start,
                                          -end => $f->end);
    push(@{$aln_g->{$f->display_id}},$g);
  }
                              
  foreach my $key(keys %$aln_g) {
    my @f = @{$aln_g->{$key}};
    my @sorted = sort {$a->start <=> $b->start} @f;
		foreach my $okf(@sorted) {
      push(@$aln_ok,$okf);
    }
  }
  return($aln_ok);
}


=head2 _make_generic_prd

 Title   : _make_generic_prd

 Usage   : my $prdf = _make_generic_prd($protein_feature);

 Function: 

 Returns : 

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub _make_generic_prd {
   
  my $prd_f = shift;
  my @prd_g;
 
  foreach my $f(@$prd_f) {

	next unless $f->length > 2;

    my @exons = @{$f->get_all_Exons};  
    foreach my $e(@exons) {    
      my $g = Bio::SeqFeature::Generic->new(-start => $e->start,
                                            -end => $e->end);
      push(@prd_g,$g);
    }
  }
  return(\@prd_g); 
}


=head2 _test_coding

 Title   : _test_coding

 Usage   : my($coding,$length) = _test_coding($slice->length,$collection)

 Function: 

 Returns : 

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub _test_coding {
	
	my $start = shift;
	my $end = shift;
  my $coll = shift;

  my @coding = $coll->features_in_range(-start => $start,
                                        -end => $end, 
                                        -contain => 0 );

  my $test = Bio::Range->new(-start => $start,
                             -end => $end );

  if(@coding) {  
    my $length = _calculate_coding(\@coding,$test); 
		return(scalar(@coding),$length);
	}
}


=head2 _calculate_coding

 Title   : _calculate_coding

 Usage   : my $length = _calculate_coding(\@coding,$test);

 Function: 

 Returns :

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub _calculate_coding {

  my $coding = shift;
  my $test = shift;
  my @coding = sort{$a->start<=>$b->start}@$coding;
  my @range;
  my $f1 = shift(@coding);
 
  foreach my $f2(@coding) {  
    if($f1->contains($f2) || $f1->equals($f2)) {
    }  
    elsif($f2->contains($f1)) {
      $f1->start($f2->start);
      $f1->end($f2->end);
    }
    elsif($f1->overlaps($f2)) {
      $f1 = $f1->intersection($f2);
    } 
    else {
      push(@range,$f1);
      $f1->start($f2->start);
      $f1->end($f2->end);
    }
  }

  push(@range,$f1);
  my $length = 0;
  foreach my $r(@range) {
    $length = $length + $test->intersection($r)->length;
  }
  return($length);
}


=head2 _drop_super_feature

 Title   : _drop_super_feature

 Usage   : 

 Function: To solve the ensembl bug that sometimes put also the entire hit with intron as DnaAlign
           feature when mapping EST. This is not rigth so I have created this method that clean
           an array of features from the wrong superfeature

 Returns : NOT MORE USED!

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub _drop_super_feature {

  my $sortedf = shift;
  my @ok;
  my $start = 0;
  my $end = 0;

  foreach(@$sortedf) {
    $start = $_->start if($start>$_->start || $start == 0);
    $end = $_->end if($end<$_->end || $end  == 0);
  }

  my $DROP=1;

  foreach(@$sortedf) {
    if($_->start==$start && $_->end==$end && $DROP==1 && scalar(@$sortedf)>1) {
      $DROP++;
    }
    elsif($_->length > 50000) {
      $DROP++;
    }
    else{
      push(@ok,$_);
    }
  }
  return(\@ok);
}


=head2 get_chr_length

 Title   : get_chr_length

 Usage   : my $length = Bio::MCE::Pipeline->get_chr_length($specie,$coor_sys_name,$seq_region_name);

 Function: return the length of a chromosome

 Returns : a scalar

 Args    : -1 the specie
           -2 the coordinate system name
					 -3 the seq region name

 NOTE    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub get_chr_length {

	my $self = shift;
  my $specie = shift;
	my $coord_sys_name = shift;
	my $seq_region_name = shift;

  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice");

  my $slice = $slice_adaptor->fetch_by_region($coord_sys_name,$seq_region_name);

  return($slice->length);
}


=head2 map_transcript

 Title   : map_transcript

 Usage   : my $collection = Bio::MCE::Pipeline->map_transcript($transcript,$specie);

 Function: buile a map of a transcript returning a collection with all the
           feature: utr3, utr5, utr_int, pre_gene, pre_atg, pre_tss, exon
           intron, post_gene, post_stop, post_end with their coordinate
           relatively to the chromosome on which they are mapped.
           If the transcript contains utr5 then it can have feature with
           primary tag: utr5, utr5_int, pre_tss, pre_atg; if not it contains
           pre_gene. If the transcript contains utr3 then it can contains
           feature with primary tag as utr3, utr_int, post_stop, post_end;
           if not it contains post_gene.

 Returns : A Bio::SeqFeature::Collection object

 Args    : -1 a Bio::EnsEMBL::Transcript object
           -2 The name of the specie on which to work in the notation
              ad example: Mus musculus

 NOTE    : THIS CAN BE BETTER....

=cut

sub map_transcript {

	my $self = shift;
  my $transcript = shift;
	my $specie = shift;
	my $length = Bio::MCE::Pipeline::Coding->get_chr_length($specie,$transcript->coord_system_name,$transcript->seq_region_name);
	my $seq_region_strand = $transcript->seq_region_strand;
  my($exon,$utr) = Bio::MCE::Utils->get_exon_utr($transcript);
	my @map;
	my @sorted_map;

  foreach my $e(@$exon) {
		my $ne;
		$ne->{'primary_tag'} = 'exon';
		$ne->{'start'} = $e->seq_region_start;
		$ne->{'end'} = $e->seq_region_end;
		push(@map,$ne);
	}
	foreach my $e(@$utr) {
		my $ne;
		$ne->{'primary_tag'} = 'utr';
		$ne->{'start'} = $e->seq_region_start;
		$ne->{'end'} = $e->seq_region_end;
		push(@map,$ne);
	}

	my @for_sorted_map = sort {$a->{'start'} <=> $b->{'start'}} @map;
	my @rev_sorted_map = sort {$b->{'start'} <=> $a->{'start'}} @map;
	@sorted_map = @for_sorted_map if $transcript->strand == 1;
	@sorted_map = @rev_sorted_map if $transcript->strand == -1;

	my @pre_exon;
	my @pre_utr5;
	my @pre_utr3;
	my @intron;

	foreach my $feat(@sorted_map) {
		
		if($feat->{'primary_tag'} eq 'exon') {
			print $transcript->stable_id." has UTR in the middle!\n" if @pre_utr3;
			push(@pre_exon,$feat);
		}
		elsif($feat->{'primary_tag'} eq 'utr') {
			if(@pre_exon) {
				push(@pre_utr3,$feat);
			}
			else {
				push(@pre_utr5,$feat);
			}
		}
	}

	my @utr5 = sort {$a->{'start'} <=> $b->{'start'}} @pre_utr5;
	my @utr3 = sort {$a->{'start'} <=> $b->{'start'}} @pre_utr3;
	my @exon = sort {$a->{'start'} <=> $b->{'start'}} @pre_exon;

	my $intron = _make_intron(\@exon);
	my $map_in = _map_intron($intron,$seq_region_strand);
	my $coll;
	
	if(@utr5 and @utr3,) {
		$coll = _make_region_utr_5_3($seq_region_strand,\@exon,\@utr5,\@utr3,$map_in,$length);
	} 
	elsif(@utr5 and not @utr3) {
		$coll = _make_region_utr_5($seq_region_strand,\@exon,\@utr5,$map_in,$length);
	} 
	elsif(@utr3 and not @utr5) {
		$coll = _make_region_utr_3($seq_region_strand,\@exon,\@utr3,$map_in,$length);
	}	
	elsif(not @utr5 and not @utr3) {		
		$coll = _make_region_noutr($seq_region_strand,\@exon,$map_in,$length);
	}
	return($coll);
}


=head2 map_non_protein_coding

 Title   : map_non_protein_coding

 Usage   : my $collection = Bio::MCE::Pipeline::Coding->map_non_protein_coding($gene,$specie);

 Function: build a map of a noncoding gene returning a collection with all the
           exons

 Returns : A Bio::SeqFeature::Collection object

 Args    : -1 a Bio::EnsEMBL::Gene object
           -2 The name of the specie on which to work in the notation
              ad example: Mus musculus

 NOTE    : THIS CAN BE BETTER....

=cut

sub map_non_protein_coding {

	my $self = shift;
  my $transcript = shift;
	my $specie = shift;
	my $length = Bio::MCE::Pipeline::Coding->get_chr_length($specie,$transcript->coord_system_name,$transcript->seq_region_name);
	my $seq_region_strand = $transcript->seq_region_strand;
  my $exon = $transcript->get_all_Exons;
	my @map;
	my @sorted_map;

  foreach my $e(@$exon) {
		my $ne;
		$ne->{'primary_tag'} = 'exon';
		$ne->{'start'} = $e->seq_region_start;
		$ne->{'end'} = $e->seq_region_end;
		push(@map,$ne);
	}

	my @for_sorted_map = sort {$a->{'start'} <=> $b->{'start'}} @map;
	my @rev_sorted_map = sort {$b->{'start'} <=> $a->{'start'}} @map;
	@sorted_map = @for_sorted_map if $transcript->strand == 1;
	@sorted_map = @rev_sorted_map if $transcript->strand == -1;

	my @pre_exon;
	my @intron;

	foreach my $feat(@sorted_map) {		
		if($feat->{'primary_tag'} eq 'exon') {
			push(@pre_exon,$feat);
		}
	}

	my @exon = sort {$a->{'start'} <=> $b->{'start'}} @pre_exon;

	my $intron = _make_intron(\@exon);
	my $map_in = _map_intron($intron,$seq_region_strand);
	
	my $coll = _make_region_noutr($seq_region_strand,\@exon,$map_in,$length);

	return($coll);
}


=head2 _make_intron

 Title   : _make_intron

 Usage   : my $intron = _make_intron(\@exon);

 Function: 

 Returns :

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub _make_intron {

	my $exon = shift;
  my $in;
  my $c = 0;

  foreach my $ex(@$exon) {
    if($in->[$c]->{'start'}) {
      $in->[$c]->{'end'} = $ex->{'start'} - 1;
      $c++;
      $in->[$c]->{'start'} = $ex->{'end'} + 1;
    }

    else {
      $in->[$c]->{'start'} = $ex->{'end'} + 1;
    }
  }
  pop(@$in); 
  return $in;
}


=head2 _map_intron

 Title   : _map_intron

 Usage   : my $numbered_intron= _map_intron(\@exon,$seq_region_strand);

 Function: 

 Returns :

 Args    : -1 the intron arrayref of hashref
           -2 the strand of the seq_region on which the transcript is mapped

 NOTE    : INTERNAL METHOD

=cut

sub _map_intron {
	
	my $intron = shift;
	my $strand = shift;
	my $num = scalar(@$intron);
	my @in;
	my $c;
	
	if($strand == 1) { $c = 0 }
	elsif($strand == -1) { $c = $num }
	else { die "\nProblema di Strand\n" }
	
	foreach my $in(@$intron) {		
		if($strand == 1) { $c++ }
		my $pt = 'intron'.$c;		
		my $fin = Bio::SeqFeature::Generic->new(-primary_tag => $pt,
																						-start => $in->{'start'},
																						-end => $in->{'end'});
		push @in,$fin;
		if($strand == -1) { $c-- }
	}	
	return \@in;
}


=head2 _make_region_utr_5_3

 Title   : _make_region_utr_5_3

 Usage   : my $collection = _make_region_utr_5_3($strand,$exons,$utr5,$utr3,$map_int,$length)

 Function: 

 Returns :

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub _make_region_utr_5_3 {
	
	my $seq_region_strand = shift;
	my $exon = shift;
	my $utr5 = shift;
	my $utr3 = shift;
	my $map_in = shift;
	my $length = shift;	
	my $coll = Bio::SeqFeature::Collection->new;
	
	my $exon_ob = _make_obj($exon,'exon');
	$coll->add_features($exon_ob);
	
	my $utr5_ob = _make_obj($utr5,'utr5');
	my $utr3_ob = _make_obj($utr3,'utr3');
	
	my $utr5_in = _make_intron($utr5);
	my $utr5_in_ob = _make_obj($utr5_in,'utr5_int');
	
	my $utr3_in = _make_intron($utr3);
	my $utr3_in_ob = _make_obj($utr3_in,'utr3_int');

	my $pre_tss = _make_pre_tss($utr5,$length,$seq_region_strand);
	my $post_end = _make_post_end($utr3,$length,$seq_region_strand);
	my $pre_atg = _make_pre_atg($exon,$utr5,$seq_region_strand);
	my $post_stop = _make_post_stop($exon,$utr3,$seq_region_strand);

	my @zone = ($pre_tss,$post_end);
	push(@zone,$pre_atg) if $pre_atg->length > 2;
	push(@zone,$post_stop) if $post_stop->length > 2;
	$coll->add_features(\@zone);
	$coll->add_features($map_in);
	$coll->add_features($utr5_ob);
	$coll->add_features($utr3_ob);
	$coll->add_features($utr5_in_ob) if $utr5_in_ob;
	$coll->add_features($utr3_in_ob) if $utr3_in_ob;
	
	return($coll);
}


=head2 _make_region_utr_5

 Title   : _make_region_utr_5

 Usage   : my $collection = _make_region_utr_5($strand,$exon,$utr5,$map_int,$length);

 Function: 

 Returns :

 Args    : 

 NOTE    : INTERNAL METHOD

=cut
sub _make_region_utr_5 {
	
	my $seq_region_strand = shift;
	my $exon = shift;
	my $utr5 = shift;
	my $map_in = shift;
	my $length = shift;
	my $coll = Bio::SeqFeature::Collection->new;
	
	my $exon_ob = _make_obj($exon,'exon');
	$coll->add_features($exon_ob);
	
	my $utr5_ob = _make_obj($utr5,'utr5');
	
	my $utr5_in = _make_intron($utr5);
	my $utr5_in_ob = _make_obj($utr5_in,'utr5_int');

	my $pre_tss = _make_pre_tss($utr5,$length,$seq_region_strand);
	my $post_gene = _make_post_gene($exon,$length,$seq_region_strand);
	my $pre_atg = _make_pre_atg($exon,$utr5,$seq_region_strand);
	
	my @zone = ($pre_tss,$post_gene);
	push(@zone,$pre_atg) if $pre_atg->length > 2;
	$coll->add_features(\@zone);
	$coll->add_features($map_in);
	$coll->add_features($utr5_ob);
	$coll->add_features($utr5_in_ob) if $utr5_in_ob;
	
	return($coll);
}


=head2 _make_region_utr_3

 Title   : _make_region_utr_3

 Usage   : my $collection = _make_region_utr_3($trand,$exons,$utr3,$map_int,$length)

 Function: 

 Returns :

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub _make_region_utr_3 {
	
	my $seq_region_strand = shift;
	my $exon = shift;
	my $utr3 = shift;
	my $map_in = shift;
	my $length = shift;	
	my $coll = Bio::SeqFeature::Collection->new;
	
	my $exon_ob = _make_obj($exon,'exon');
	$coll->add_features($exon_ob);
	
	my $utr3_ob = _make_obj($utr3,'utr3');
	
	my $utr3_in = _make_intron($utr3);
	my $utr3_in_ob = _make_obj($utr3_in,'utr3_int');

	my $pre_gene = _make_pre_gene($exon,$length,$seq_region_strand);
	my $post_end = _make_post_end($utr3,$length,$seq_region_strand);
	my $post_stop = _make_post_stop($exon,$utr3,$seq_region_strand);
	
	my @zone = ($pre_gene,$post_end);
	push(@zone,$post_stop) if $post_stop->length > 2;
	$coll->add_features(\@zone);
	$coll->add_features($map_in);
	$coll->add_features($utr3_ob);
	$coll->add_features($utr3_in_ob) if $utr3_in_ob;
	
	return($coll);
}


=head2 _make_region_noutr

 Title   : _make_region_noutr

 Usage   : my $collection = _make_region_noutr($strand,$exon,$map_int,$length);

 Function: 

 Returns :

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub _make_region_noutr {
	
	my $seq_region_strand = shift;
	my $exon = shift;
	my $map_in = shift;
	my $length = shift;
	my $coll = Bio::SeqFeature::Collection->new;
	
	my $exon_ob = _make_obj($exon,'exon');
	$coll->add_features($exon_ob);
	
	my $pre_gene = _make_pre_gene($exon,$length,$seq_region_strand);
	my $post_gene = _make_post_gene($exon,$length,$seq_region_strand);
	
	my @zone = ($pre_gene,$post_gene);
	$coll->add_features(\@zone);
	$coll->add_features($map_in);
	
	return($coll);
}


=head2 _make_obj

 Title   : _make_obj

 Usage   : my $seq_features = _make_obj(

 Function: build the generic features that will be inserted into the collection

 Returns : an arrayref of Bio::SeqFeature::Generic

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub _make_obj {
	
	my $list = shift;
	my $pt = shift;	
	my @objs;
	
	foreach my $data(@{$list}) {
		
		my $obj = Bio::SeqFeature::Generic->new(-primary_tag => $pt,
																						-start => $data->{'start'},
																						-end => $data->{'end'});
		push(@objs,$obj);
	}
	return(\@objs);
}


=head2 _make_pre_tss

 Title   : _make_pre_tss

 Usage   : my $collection = _make_pre_tss($utr5,$length,$strand);

 Function: build generic feature that will be added to the collection
           for the pre_tss region

 Returns : Bio::SeqFeature::Generic object

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub	_make_pre_tss {
	
	my $utr5 = shift;
	my $length = shift;
	my $strand = shift;
	my $pre_tss;
	
	if($strand == 1) {
		# from the start of the seq to the start of the first utr5
		$pre_tss = Bio::SeqFeature::Generic->new(-primary_tag => 'pre_tss',
																						 -start => 1,
																						 -end => ($utr5->[0]->{'start'} -1));
	}
	elsif($strand == -1) {
		# from the end of the last utr5 to the end of the seq
		$pre_tss = Bio::SeqFeature::Generic->new(-primary_tag => 'pre_tss',
																						 -start => ($utr5->[$#$utr5]->{'end'} + 1),
																						 -end => $length);
	}
	return($pre_tss);
}
	
sub	_make_pre_atg {

	my $exon = shift;
	my $utr5 = shift;
	my $strand = shift;
	my $pre_atg;
	
	if($strand == 1) {
		# from the end of the last utr5 to the start of the first exon. 
		# if there are more than one utr5 that will be called utr5
		$pre_atg = Bio::SeqFeature::Generic->new(-primary_tag => 'pre_atg',
																						 -start => ($utr5->[$#$utr5]->{'end'} + 1),
																						 -end => ($exon->[0]->{'start'} - 1));
	}
	elsif($strand == -1) {
		# from the end of the last exon to the start of the first utr5.
		# if there are more than one utr5 that will be called utr5
	  $pre_atg = Bio::SeqFeature::Generic->new(-primary_tag => 'pre_atg',
																						 -start => ($exon->[$#$exon]->{'end'} + 1),
																						 -end => ($utr5->[0]->{'start'} - 1));
	}	
	return($pre_atg);
}


=head2 _make_pre_gene

 Title   : _make_pre_gene

 Usage   : my $pre_gene = _make_pre_gene($exons,$length,$strand);

 Function: build generic feature that will be added to the collection
           for the pre_gene region

 Returns : Bio::SeqFeature::Generic object

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub	_make_pre_gene {
	
	my $exon = shift;
	my $length = shift;
	my $strand = shift;
	my $pre_gene;
	
	if($strand == 1) {
		# from the start of the seq to the start of the first exon
		$pre_gene = Bio::SeqFeature::Generic->new(-primary_tag => 'pre_gene',
																							-start => 1,
																							-end => ($exon->[0]->{'start'} -1));
	}
	elsif($strand == -1) {
		# from the end of the last exon to the end of the sequence
		$pre_gene = Bio::SeqFeature::Generic->new(-primary_tag => 'pre_gene',
																							-start => ($exon->[$#$exon]->{'end'} + 1),
																							-end => $length);
	}
	return($pre_gene);
}

sub	_make_post_gene {
	
	my $exon = shift;
	my $length = shift;
	my $strand = shift;
	my $post_gene;
	
	if($strand == 1) {
		# from the end of the last exon to the end of the sequence
		$post_gene = Bio::SeqFeature::Generic->new(-primary_tag => 'post_gene',
																							 -start => ($exon->[$#$exon]->{'end'} + 1),
																							 -end => $length);
	}
	elsif($strand == -1) {
		# from the start of the seq to the start of the first exon
		$post_gene = Bio::SeqFeature::Generic->new(-primary_tag => 'post_gene',
																							 -start => 1,
																							 -end => ($exon->[0]->{'start'} -1));
	}
	return($post_gene);
}	


=head2 _make_post_end

 Title   : _make_post_end

 Usage   : my $collection = _make_post_end($utr3,$length,$strand);

 Function: build generic feature that will be added to the collection
           for the post_end region

 Returns : Bio::SeqFeature::Generic object

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub	_make_post_end {

	my $utr3 = shift;
	my $length = shift;
	my $strand = shift;
	my $post_end;

	if($strand == 1) {
		# from the end of the last utr3 to the end of the seq
		$post_end = Bio::SeqFeature::Generic->new(-primary_tag => 'post_end',
																							-start => ($utr3->[$#$utr3]->{'end'} + 1),
																							-end => $length);
	}
	elsif($strand == -1) {
		# from the start of the seq to the start of the first utr3
		$post_end = Bio::SeqFeature::Generic->new(-primary_tag => 'post_end',
																							-start => 1,
																							-end => ($utr3->[0]->{'start'} -1));
	}
	return($post_end);
}

=head2 _make_post_stop

 Title   : _make_post_stop

 Usage   : my $collection = _make_post_stop($exons,$utr3,$strand);

 Function: build generic feature that will be added to the collection
           for the pre_tss region

 Returns : Bio::SeqFeature::Generic object

 Args    : 

 NOTE    : INTERNAL METHOD

=cut

sub	_make_post_stop {
	
	my $exon = shift;
	my $utr3 = shift;
	my $strand = shift;
	my $post_stop;
	
	if($strand == 1) {
		# from the end of the last exon to the start of the first utr3.
		# if there are more than one utr3 that will be called utr3
		$post_stop = Bio::SeqFeature::Generic->new(-primary_tag => 'post_stop',
																							 -start => ($exon->[$#$exon]->{'end'} + 1),
																							 -end => ($utr3->[0]->{'start'} - 1));
	}
	elsif($strand == -1) {
		# from the end of the last utr3 to the start of the first exon. 
		# if there are more than one utr3 that will be called utr3
		$post_stop = Bio::SeqFeature::Generic->new(-primary_tag => 'post_stop',
																							 -start => ($utr3->[$#$utr3]->{'end'} + 1),
																							 -end => ($exon->[0]->{'start'} - 1));
	}
	return($post_stop);
}

1;
