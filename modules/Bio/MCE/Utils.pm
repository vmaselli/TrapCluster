#!/usr/bin/perl
package Bio::MCE::Utils;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Exon;
#use Bio::Factory::EMBOSS;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Range;
use Bio::AlignIO;

use Bio::MCE::Pipeline::Coding;

=head2 get_all_homologous

 Title   : get_all_homologous

 Usage   : my $hom = Bio::MCE::Utils->get_all_homologous($id,$specie);

 Function: retrieve all homologous genes from compara

 Returns : listref of [Bio::EnsEMBL::Gene,Bio::EnsEMBL::Compara::Homology,string $specie]

 Args    : -1 gene_stable_id
           -2 specie

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub get_all_homologous {

  my $class = shift;
  my $id = shift;
  my $specie = shift;

	print "class $class id $id specie $specie\n";

  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","gene") or die "ERROR 2 $!\n";

  my $gene = $gene_adaptor->fetch_by_stable_id($id) or die "ERROR 3 $!\n";

  my $hom = $gene->get_all_homologous_Genes;

	return($hom);
}


=head2 _regex_id_genome

 Title   : _regex_id_genome

 Usage   : my $id_genome = _regex_id_genome();

 Function: it is a function in order to have a regular expression to
           link a gene stable id to a species !!!INTERNAL METHOD!!!!

           ACTUAL HASH:

  my %hash_id = ('ENSG0' => 'Homo sapiens',
                 'ENSMUSG' => 'Mus musculus',
                 'ENSRNOG' => 'Rattus norvegicus',
#                 'NEWSINFRUG' => 'Takifugu rubripes',
#                 'SINFRUG' => 'Takifugu rubripes',
                 'ENSTRUG' => 'Takifugu rubripes',
                 'ENSDARG' => 'Danio rerio',
                 'ENSPTRG' => 'Pan troglodytes',
                 'ENSGALG' => 'Gallus gallus',
                 'ENSCAFG' => 'Canis familiaris',
                 'ENSCING' => 'Ciona intestinalis',
                 'ENSXETG' => 'Xenopus tropicalis',
                 'AY' => 'Tetraodon nigroviridis',
                 'EVX' => 'Tetraodon nigroviridis',
                 'GSTENG' => 'Tetraodon nigroviridis',
                 'HOX' => 'Tetraodon nigroviridis',
                 'ENSBTAG' => 'Bos taurus',
                 'ENSMODG' => 'Monodelphis domestica',
                 'ENSMMUG' => 'Macaca mulatta',
                 'ENSGACG' => 'Gasterosteus aculeatus',
                 'ENSCSAVG0' => 'Ciona savignyi',
                 'ENSORLG0' => 'Oryzias latipes');

 Returns : hasref of id_regex => genome

 Args    : no

 NOTE    : INTERNAL METHOD
					 TO BE UPDATED WHEN ADD OR DELETE SPECIES FROM ENSEMBL

=cut

sub _regex_id_genome {

	my %hash_id = ('ENSG0' => 'Homo sapiens',
								 'ENSMUSG' => 'Mus musculus',
								 'ENSRNOG' => 'Rattus norvegicus',
#								 'NEWSINFRUG' => 'Takifugu rubripes',
#                 'SINFRUG' => 'Takifugu rubripes',
                 'ENSTRUG' => 'Takifugu rubripes',
								 'ENSDARG' => 'Danio rerio',
								 'ENSPTRG' => 'Pan troglodytes',
								 'ENSGALG' => 'Gallus gallus', 
								 'ENSCAFG' => 'Canis familiaris',
								 'ENSCING' => 'Ciona intestinalis',
								 'ENSXETG' => 'Xenopus tropicalis',
								 'AY' => 'Tetraodon nigroviridis',
								 'EVX' => 'Tetraodon nigroviridis',
								 'GSTENG' => 'Tetraodon nigroviridis',
								 'HOX' => 'Tetraodon nigroviridis',
								 'ENSBTAG' => 'Bos taurus',
								 'ENSMODG' => 'Monodelphis domestica',
								 'ENSMMUG' => 'Macaca mulatta',
								 'ENSGACG' => 'Gasterosteus aculeatus',
								 'ENSCSAVG0' => 'Ciona savignyi',
                 'ENSORLG0' => 'Oryzias latipes');

	return(\%hash_id);
}


=head2 get_genomes

 Title   : get_genomes

 Usage   : my $genomes = Bio::MCE::Utils->get_genomes();

 Function: If you need to know all the genomes in reg_ex hash
           that are the only that you can use in valis

 Returns : array of genomes strings with '_' separating
           genus from specie

 Args    : no

=cut

sub get_genomes {

	my $class = shift;
	my $res;
  my $hash_id = _regex_id_genome();
	my @genomes = values(%$hash_id);
	foreach my $ele(@genomes) {		
		my $o = lc($ele);
		$o =~ s/ /\_/;
		push(@$res,$o);
	}
	return $res;
}


=head2 regex_id_genome

 Title   : regex_id_genome

 Usage   : my $id_genome = Bio::MCE::Utils->regex_id_genome();

 Function: it is a function in order to have a regular expression to
           link a gene stable id to a species

 Returns : hasref of id_regex => genome

 Args    : no

=cut

sub regex_id_genome {

	my $class = shift;
	
	my $hash_id = _regex_id_genome();

	return($hash_id);
}

	
=head2 get_genome_from_id

 Title   : get_genome_from_id

 Usage   : my $id_genome = Bio::MCE::Utils->get_genome_from_id($id);

 Function: Returns the genome from the suffix of the stable_id

 Returns : a string

 Args    : ensembl stable gene id 

=cut

sub get_genome_from_id {

	my $class = shift;
	my $id = shift;
	my $c = 0;

	my $hash_id = _regex_id_genome();

	foreach my $key(keys %$hash_id) {
		if($id =~ /$key/) {

			return($hash_id->{$key});
			my $c ++;
		}
	}
	#print STDERR "\n$id not in my book!\n" and return undef unless $c;
	return undef unless $c;
}


=head2 get_prefix_from_genome

 Title   : get_prefix_from_genome

 Usage   : my $id_genome = Bio::MCE::Utils->get_prefix_from_genome($genome);

 Function: Returns the prefix of the stable_id carachteristic of the genome

 Returns : a string

 Args    : a string of the genome

 Note    : AT RISK!

           There is at least one genome in reg_ex hash in which
           the prefix are more than one. This function give back
           only one prefix.

           TO CHECK WHERE THIS IS USED AND TO CHANGE!!!!

=cut

sub get_prefix_from_genome {

	my $class = shift;
  my $genome = shift;
	my $hash = _regex_id_genome;

	foreach my $prefix(keys %$hash) {
		if($genome eq $hash->{$prefix}) {
			return $prefix;
		}
	}
#  die "\n$genome not in my book!\n";
}


=head2 seq_perc

 Title   : seq_perc

 Usage   : my $perc = Bio::MCE::Utils->seq_perc($string,$letter);

 Function: it is a function to have the percentage of 
           composition of $letter in the $string

 Returns : a number

 Args    : -1 a string representing a sequence
           -2 a single char representing the letter of wich
              to know the composition

=cut

sub seq_perc {

  my $class = shift;
	my $pre_string = shift;
	my $pre_letter = shift;

	my $string = uc($pre_string);
	my $letter = uc($pre_letter);
	
	my @n = ($string =~ m/$letter/g);
	my $n = scalar(@n);
	my $n_perc;
    
	if($n) {
		$n_perc = (($n/(length($string)))*100);
  }
  else {
		$n_perc = 0;
  }
	return $n_perc;
}


=head2 seq_perc_2_chr

 Title   : seq_perc_2_chr

 Usage   : my $perc = Bio::MCE::Utils->seq_perc_2_chr($string,$letter1,$letter2);

 Function: it is a function to have the percentage of
           composition of $letter1 + $letter2 in the $string

 Returns : a number

 Args    : -1 a string representing a sequence
           -2 a single char representing the letter of wich
              to know the composition
           -3 a single char representing the letter of wich
              to know the composition

=cut

sub seq_perc_2_chr {

  my $class = shift;
  my $pre_string = shift;
  my $pre_letter1 = shift;
  my $pre_letter2 = shift;

  my $string = uc($pre_string);
  my $letter1 = uc($pre_letter1);
  my $letter2 = uc($pre_letter2);
 
  my @n1 = ($string =~ m/$letter1/g);
  my @n2 = ($string =~ m/$letter2/g);
  my $n = scalar(@n1) + scalar(@n2);
  my $n_perc;
 
  if($n) {
    $n_perc = (($n/(length($string)))*100);
  }
  else {
    $n_perc = 0;
  }
  return $n_perc;
}


=head2 get_exon_utr

 Title   : get_exon_utr

 Usage   : my($exon,$utr) = Bio::MCE::Utils->get_exon_utr($transcript);

 Function: map and create two array containing the exons object
           representing the exons in the first and the exon 
           objects representing the utr in the second

 Returns : two array_ref (\@EXONS,\@UTR)

 Args    : an ensembl transcript

=cut

sub get_exon_utr {

	my $class = shift;	# Transcript
  my $self = shift;
	my $end_exon_check = 0;
	my(@EXON);
	my(@UTR);

	my $translation = $self->translation
	or $self->warning("No translation attached to transcript object");

	my $start_exon = $translation->start_Exon;
	my $end_exon = $translation->end_Exon;
	my $t_start = $translation->start;
	my $t_end = $translation->end;

	foreach my $ex(@{$self->get_all_Exons}) {

		if($ex ne $start_exon and ! @EXON) {
			push(@UTR,$ex);	# UTR Before the first translated exon. Not yet in translated region
			next;
		}
		
		if ($end_exon_check) {
			push(@UTR, $ex);	# UTR After the last translated exon. We are after translated region
			next;
		}
		
		my $length = $ex->length;

		my $adjust_start = 0;
		my $adjust_end = 0;
		
		my $utr5_start = 0;
		my $utr5_end = 0;
		my $utr3_start = 0;
		my $utr3_end = 0;

		### HERE IN THE NEXT LINES THERE IS STRANGE HACK
		### NOT SURE BUT SHOULD BE DUE TO THE FACT THAT IN THE
		### ADJUST_START_END FUNCTION OF THE EXON MODULE THERE
		### ARE OPERATION WITHOUT '+1' OR '-1'....
		### SO HERE YOU WILL FIND '0' OR '-1' BUT NOT AS YOU EXPECT
		### AND NO '+1'

		# Adjust to translation start if this is the start exon
		if ($ex == $start_exon) {
			if ($t_start < 1 or $t_start > $length) {
				$self->throw("Translation start '$t_start' is outside exon $ex length=$length");
			}
			$adjust_start = $t_start - 1; #+ 1; #QUI C'ERA ANCHE - 1 ERA SBAGLIATO!?!?;

			if ($t_start > 1) {
				$utr5_start = 0;
				$utr5_end = $adjust_start - $length;
			}
		}

		# Adjust to translation end if this is the end exon
		if ($ex == $end_exon) {
			if ($t_end < 1 or $t_end > $length) {
				$self->throw("Translation end '$t_end' is outside exon $ex length=$length");
			}	
			$adjust_end = $t_end - $length; #QUI +1 NON C'ERA..

			if ($t_end < $length) {
				$utr3_start = $t_end + 1;# - 1;
				$utr3_end = 0;
			}
		}

		# Make a truncated exon if the translation start or
		# end causes the coordinates to be altered.
		if ($adjust_end || $adjust_start) {
			my $newex = exon_adjust_start_end($ex,$adjust_start,$adjust_end);
			push(@EXON, $newex) unless $end_exon_check;
			$end_exon_check ++ if $ex == $end_exon;
		}
		else {
			push(@EXON, $ex) unless $end_exon_check;
			$end_exon_check ++ if $ex == $end_exon;
		}

		if ($utr5_start || $utr5_end) {
			my $newutr = exon_adjust_start_end($ex,$utr5_start,$utr5_end);
			push(@UTR, $newutr);
		} 
		if ($utr3_start || $utr3_end) {
			my $newutr = exon_adjust_start_end($ex,$utr3_start,$utr3_end);
			push(@UTR, $newutr);
		} 
	}

	return(\@EXON, \@UTR);
}


=head2 exon_adjust_start_end

 Title   : exon_adjust_start_end

 Usage   : my $e = Bio::MCE::Utils->exon_adjust_start_end($exon,$start_adj,$end_adj);

 Function: 

 Returns : a new exon with adjusted coordinates

 Args    : -1 an ensembl exon object
           -2 start adjust
           -3 end adjust

 Note    : This functions has been just cutted and pasted from ensembl
           in order to mantain it stable because is used by get_exon_utr

=cut

sub exon_adjust_start_end  {

  my ( $self, $start_adjust, $end_adjust ) = @_;

  my $new_exon = Bio::EnsEMBL::Exon->new();
  %{$new_exon} = %{$self};

  #invalidate the sequence cache
  delete $new_exon->{'_seq_cache'};

  if( $self->strand() == 1 ) {
    $new_exon->start( $self->start() + $start_adjust );
    $new_exon->end( $self->end() + $end_adjust )
  } else {
    $new_exon->start( $self->start() - $end_adjust );
    $new_exon->end( $self->end() - $start_adjust )
  }
  return $new_exon;
}


=head2 limit_to_slice

 Title    : limit_to_slice

 Usage    : limit_to_slice($slice,$exon)

 Function : map exon or every type of obj with start and end function on the slice 
            $slice and return start and end adjusted to the start and end of the slice.
            This is why we need all coords between start-end of slice.
 
 Example  : my($start,$end) = limit_to_slice($slice,$exon);
 
 Returns  : $start,$end are scalar adjusted to be between slice start-end coords
 
 Args     : -1 slice obj
            -2 exon obj
 
 Note     : THIS SUB IS USED BY MAKE_SEQOBJ in ortoinspector
            Here should consider some other cases.....

=cut

sub limit_to_slice {
	
	my $class = shift;
	my $slice = shift;
	my $f = shift;
	my $start;
	my $end;
	
	if(($f->start <= 1 && $f->end <= 1) || ($f->start >= $slice->length && $f->end >= $slice->length)) {
		$start = 0;
		$end = 0;
	}
	else{
		$start = $f->start;
		$end = $f->end;
		if($start<1){
			$start=1;
		}
		if($end>$slice->length){
			$end=$slice->length;
		}
	}
	return ($start,$end);
}


=head2 run_needle

 Title    : run_needle

 Usage    : Bio::MCE::Utils->run_needle($string1,$string2)

 Function : Run the program needle from the EMBOSS package

 Example  : my $percentage = run_needle($string1,$string2);

 Returns  : a scalar indicating the percentage identity

 Args     : -1 Bio::Seq object
            -2 Bio::Seq object

=cut

sub run_needle {
	
	my $class = shift;
	my $needle1 = shift;
	my $needle2 = shift;
	my $factory = Bio::Factory::EMBOSS->new;
	my $needle = $factory->program('needle') or die "WHERE IS NEEDLE?\n";
	my @arrayall;
	push @arrayall,$needle2;
	my $needleout = $needle1->display_id.$needle2->display_id.".needle";

	$needle->run({'-asequence' => $needle1,
								'-bsequence' => $needle2,
								'-outfile'   => $needleout});

	# now you might want to get the alignment
	open(OUT,"<$needleout") or die "\nPRTOBLEM: $!\n";

	my $perc;

	while (my $line = <OUT>) {
		next unless $line =~ /Identity\:.+\/.+\((.+\..+)\%\)/;
		$perc = $1;
	}
	system 'rm *.needle';
	return $perc;
}


=head2 get_overlapping_features

 Title    : get_overlapping_features

 Usage    : Bio::MCE::Utils->get_overlapping_features($dbh,table,$ref_id,$start,$end)

 Function : Perform a query in db in table $table and returns
            arrayref of fetchrow_hashref of the results

 Example  : my @$overlapping_features = @{get_overlapping_features($dbh,$table,$ref_id,$start,$end)};

 Returns  : arrayref of hashref

 Args     : -1 db handle
            -2 table
            -3 ref_id
            -4 start
            -5 end

=cut

sub get_overlapping_features {

	my $class = shift;
	my $dbh = shift;
	my $table = shift;
	my $ref_id = shift;
	my $start = shift;
	my $end = shift;
	my @f;

	my $sth = $dbh->prepare('SELECT * FROM '.$table.
													' WHERE ref_id = '.$dbh->quote($ref_id).
													' AND start <= '.$dbh->quote($end).' AND end >= '.$dbh->quote($start));
	$sth->execute;

	while(my $f = $sth->fetchrow_hashref) {
		push(@f,$f);
	}
	return(\@f);
}


=head2 get_overlapping_features_by_loc

 Title    : get_overlapping_features_by_loc

 Usage    : Bio::MCE::Utils->get_overlapping_features_by_loc($dbh,table,$chr,$chr_start,$chr_end)

 Function : Perform a query in db in table $table and returns
            arrayref of fetchrow_hashref of the results

 Example  : my @$overlapping_features = @{get_overlapping_features_by_loc($dbh,$table,$chr,$chr_start,$chr_end)};

 Returns  : arrayref of hashref

 Args     : -1 db handle
            -2 table
            -3 chr
            -4 chr_start 
            -5 chr_end

=cut

sub get_overlapping_features_by_loc {

  my $class = shift;
  my $dbh = shift;
  my $table = shift;
  my $chr = shift;
  my $chr_start = shift;
  my $chr_end = shift;
  my @f;

  my $sth = $dbh->prepare('SELECT * FROM '.$table.
                          ' WHERE chr = '.$dbh->quote($chr).
                          ' AND start <= '.$dbh->quote($chr_end).' AND end >= '.$dbh->quote($chr_start));
  $sth->execute;

  while(my $f = $sth->fetchrow_hashref) {
    push(@f,$f);
  }
  return(\@f);
}


=head2 get_overlapping_feature_pairs

 Title    : get_overlapping_feature_pairs

 Usage    : Bio::MCE::Utils->get_overlapping_feature_pairs($dbh,table,$q_ref_id,$q_start,$q_end)

 Function : Perform a query in db in table $table and returns
            arrayref of fetchrow_hashref of the results

 Example  : my @$overlapping_feature_pairs = @{get_overlapping_feature_pairs($dbh,$table,$q_ref_id,$q_start,$q_end)};

 Returns  : arrayref of hashref

 Args     : -1 db handle
            -2 table
            -3 q_ref_id
            -4 q_start
            -5 q_end

=cut

sub get_overlapping_feature_pairs {

	my $class = shift;
	my $dbh = shift;
	my $table = shift;
	my $q_ref_id = shift;
	my $start = shift;
	my $end = shift;
	my @f;

	my $sth = $dbh->prepare('SELECT * FROM '.$table.
													' WHERE q_ref_id = '.$dbh->quote($q_ref_id).
													' AND q_start <= '.$dbh->quote($end).' AND q_end >= '.$dbh->quote($start));
	$sth->execute;

	while(my $f = $sth->fetchrow_hashref) {
		push(@f,$f);
	}
	return(\@f);
}


=head2 ref_to_chr

 Title    : ref_to_chr

 Usage    : Bio::MCE::Utils->ref_to_chr($dbh,$ref_id,$feature)

 Function : change coordinates from the ref_id to the chromosomic

 Example  : my($chr,$chr_start,$chr_end) = Bio::MCE::Utils->ref_to_chr($dbh,$ref_id,$feature);

 Returns  : -1 chr (string)
            -2 chr_start (scalar)
            -3 chr_end (scalar)

 Args     : -1 dbh
            -2 ref_if
            -3 feature obj to move

=cut

sub ref_to_chr {

  my $class = shift;
	my $dbh = shift;
	my $ref_id = shift;
	my $feature = shift;

	my $sth = $dbh->prepare('SELECT ref_id, chr, chr_start, chr_end, slice_strand '.
													'FROM seq WHERE ref_id = '.$dbh->quote($ref_id));
	$sth->execute;

	my $row = $sth->fetchrow_hashref;
	my $chr = $row->{'chr'};
	my $start = $row->{'chr_start'};
	my $end = $row->{'chr_end'};
	my $slice_strand = $row->{'slice_strand'};
	my $region_length = ($end-$start+1);
	my $rel_start;

	if($slice_strand == -1) {
		$rel_start = $region_length - $feature->end +1;
	}
	else {
		$rel_start = $feature->start;
	}

	my $chr_start = $start + $rel_start - 1;
	my $chr_end = $chr_start + $feature->length - 1;

	return($chr,$chr_start,$chr_end);
}


=head2 get_overlapping_ref

 Title    : get_overlapping_ref

 Usage    : Bio::MCE::Utils->get_overlapping_ref($dbh,$chr,$chr_start,$chr_end)

 Function : returns ref_id overlapping with the submitted cohordinates

 Example  : my($dbh,$chr,$chr_start,$chr_end) = get_overlapping_ref($dbh,$chr,$chr_start,$chr_end);

 Returns  : arrayref: ref_id

 Args     : -1 dbh
            -2 chr
            -3 chr_start
            -4 chr_end

=cut

sub get_overlapping_ref {

  my $class = shift;
	my $dbh = shift;
	my $chr = shift;
	my $chr_start = shift;
	my $chr_end = shift;
	my @ref;
	
	my $sth = $dbh->prepare('SELECT * FROM seq WHERE chr = '.$dbh->quote($chr).
													' AND chr_start <= '.$dbh->quote($chr_end).' AND chr_end >= '.$dbh->quote($chr_start));
	$sth->execute;

	while(my $row = $sth->fetchrow_hashref){
		push(@ref,($row->{'ref_id'}));
	}

	return(\@ref);
}


=head2 chr_to_ref

 Title    : chr_to_ref

 Usage    : Bio::MCE::Utils->chr_to_ref($dbh,$ref_id,$chr_start,$chr_end,$chr_strand)

 Function : change coordinates from the chormosomic to the ref_id

 Example  : my($start,$end,$strand) = ref_to_chr($dbh,$ref_id,$chr_start,$chr_end,$chr_strand);

 Returns  : feature obj

 Args     : -1 dbh
            -2 ref_if
            -3 start
            -4 end
            -5 strand
 
 NOTE     : TO CHECK!!!!!
            NOT STABLE!!!

=cut

sub chr_to_ref {

  my $class = shift;
  my $dbh = shift;
  my $ref_id = shift;
  my $start = shift;
	my $end = shift;
	my $chr_strand = shift;

  my $sth = $dbh->prepare('SELECT ref_id, chr, chr_start, chr_end, slice_strand '.
                          'FROM seq WHERE ref_id = '.$dbh->quote($ref_id));
  $sth->execute;

  my $row = $sth->fetchrow_hashref;
  my $chr = $row->{'chr'};
  my $ref_start = $row->{'chr_start'};
  my $ref_end = $row->{'chr_end'};
	my $slice_strand = $row->{'slice_strand'};
	my $length = $end - $start + 1;
  my $f_start;
  my $f_end;

	# CASE 1
	#	ref		|--------------------|
	#	feat      |---|
	if($ref_start<=$start and $ref_end>=$end){
		$f_start = $start - $ref_start + 1;
		$f_end = $end - $ref_start + 1;
	}

  # CASE 2
  # ref      |--------------------|
  # feat   |---|
	if($start<=$ref_start and $end>=$ref_start and $end<=$ref_end){
		$f_start = $ref_start;
		$f_end = $end - $ref_start + 1;
    print "CASE 2 REDUCTION!\n";
	}

  # CASE 3
  # ref   |--------------------|
  # feat      							 |---|
	if($start>=$ref_start and $start<=$ref_end and $end>=$ref_end){
		$f_start = $start - $ref_start + 1;
		$f_end = $ref_end;
    print "CASE 3 REDUCTION!\n";
	}

  # CASE 4
  # ref   		|------------|
  # feat  |---------------------|
	if($ref_start>=$start and $ref_end<=$end){
		$f_start = $ref_start;
		$f_end = $ref_end;
		print "CASE 4 REDUCTION!\n";
	}

	my $new_strand = ($slice_strand * $chr_strand);

	if($slice_strand == -1) {
		my $flipped_start = $length - $f_end + 1;
		my $flipped_end = $flipped_start + ($f_end - $f_start + 1) - 1;
		return($flipped_start,$flipped_end,$new_strand);
	}
	else {
	  return($f_start,$f_end,$new_strand);
	}
}


=head2 write_cst_organism_fasta

 Title    : write_cst_organism_fasta

 Usage    : Bio::MCE::Utils->write_cst_organism_fasta($dbh,$feature_id,$number,$type)

 Function : write a fasta file in which there are all the fragment
            in a final cst for all considered organism according to the
             forward strand from the table 'cst_organism'.

 Example  : my $fasta = write_cst_organism($dbh,$feature_id,'1','final');

 Returns  : name of the generated file

 Args     : -1 dbh
            -2 feature_id of a final_cst
            -3 'number' of seq foreach specie in this case you 
               will get the ones with the higer identity*length
            -4 type of cst eg 'final' or 'vista'

 NOTE     : Tha names of the sequences are in the format: 'cst_organism_id-ref_id'

=cut

sub write_cst_organism_fasta {

  my $class = shift;
	my $dbh = shift;
	my $cst = shift;
	my $singleseq = shift;
  my $type = shift;
	my $seqref;
	my $valref;
	my $fasta = Bio::SeqIO->new(-file => ">$type\_cst_$cst\.fa",
															-format => 'fasta');

	my $sth = $dbh->prepare('SELECT c.*, SUBSTRING(sequence,start,end-start+1) string '.
													'FROM cst_organism c, seq s '.
													'WHERE c.ref_id=s.ref_id '.
													'AND feature_id = '.$dbh->quote($cst).
													' AND conserved = 1');
	$sth->execute;

	while(my $row = $sth->fetchrow_hashref) {

		if($row->{'feature_pair'}){
			my $fpid = $row->{'fature_pair'};
			my $query = $dbh->prepare('SELECT ((identity/100)*(q_end-q_start+1)) val FROM feature_pair '.
																'WHERE feature_pair_id = '.$dbh->quote($fpid));
			$query->execute;
			my $fh = $query->fetchrow_hashref;
			$row->{'val'} = $fh->{'val'};
			push(@{$seqref->{$row->{'ref_id'}}},$row);
		}
		else{
			$row->{'val'} = 100;
			push(@{$seqref->{$row->{'ref_id'}}},$row);
		}
	}			

	foreach my $k(keys %$seqref){
		my @ordered = sort{$b->{'val'}<=>$a->{'val'}}@{$seqref->{$k}};

		my $x;
		for($x=1;$x<=$singleseq;$x++){
			next unless $ordered[$x-1];
			my $row = $ordered[$x-1];
			my $string = $row->{'string'};
			my $strand = $row->{'strand'};
			my $cst_id = $row->{'cst_organism_id'};
			my $ref_id = $row->{'ref_id'};

			my $seq = Bio::Seq->new(-seq => $string,
															-id => "$cst_id\-$ref_id");
			my $seqobj;

			if($strand == '-1'){
				$seqobj = $seq->revcom;
			}
			else{
				$seqobj = $seq;
			}
			$fasta->write_seq($seqobj);
		}
	}	
	return("$type\_cst_$cst\.fa")
}


=head2 write_cst_organism_fasta_from_hit

 Title    : write_cst_organism_fasta_from_hit

 Usage    : Bio::MCE::Utils->write_cst_organism_fasta_from_hit($dbh,$feature_id,$singleseq,$type)

 Function : write a fasta file in which there are all the hits (feature_pair)
            overlapping a cst from the table 'cst_organism'

 Example  : my $fasta = write_cst_organism_from_hit($dbh,$feature_id,'1','final');

 Returns  : name of the generated file

 Args     : -1 dbh
            -2 feature_id of a final_cst
            -3 'number' of seq foreach specie in this case you
               will get the ones with the higer identity*length
            -4 type of cst eg 'final' or 'vista'
 
 NOTE     : -1 The names of the sequences are in the format: 'cst_organism_id-ref_id'
            -2 ADD the possibility to exclude species

=cut

sub write_cst_organism_fasta_from_hit {

  my $class = shift;
	my $dbh = shift;
	my $cst = shift;
	my $singleseq = shift;
  my $type = shift;
	#	my $left = shift;
	#	my $rigth = shift;
	my $seqref;
	my $valref;
	my $fasta = Bio::SeqIO->new(-file => ">$type\_cst_$cst\.fa",
															-format => 'fasta');

	my $sth = $dbh->prepare('SELECT c.*, SUBSTRING(sequence,start,end-start+1) string '.
													'FROM cst_organism c, seq s '.
													'WHERE c.ref_id=s.ref_id '.
													'AND feature_id = '.$dbh->quote($cst).
													' AND conserved = 1');
	$sth->execute;

	while(my $row = $sth->fetchrow_hashref) {

		if($row->{'feature_pair_id'}){
			my $fpid = $row->{'feature_pair_id'};
			my $query = $dbh->prepare('SELECT ((identity/100)*(q_end-q_start+1)) val, '.
																'substring(sequence,t_start,t_end-t_start+1) se '.
																'FROM feature_pair fp, seq s '.
																'WHERE feature_pair_id = '.$dbh->quote($fpid).
																' AND t_ref_id=ref_id');
			$query->execute;

			my $fh = $query->fetchrow_hashref;

			$row->{'val'} = $fh->{'val'};
			$row->{'string'} = $fh->{'se'};
			push(@{$seqref->{$row->{'ref_id'}}},$row);
		}
		else{
			$row->{'val'} = 100;
			push(@{$seqref->{$row->{'ref_id'}}},$row);
		}
	}			

	foreach my $k(keys %$seqref){
		my @ordered = sort{$b->{'val'}<=>$a->{'val'}}@{$seqref->{$k}};

		my $x;
		for($x=1;$x<=$singleseq;$x++){
			next unless $ordered[$x-1];
			my $row = $ordered[$x-1];
			my $string = $row->{'string'};
			my $strand = $row->{'strand'};
			my $cst_id = $row->{'cst_organism_id'};
			my $ref_id = $row->{'ref_id'};

			my $seq = Bio::Seq->new(-seq => $string,
															-id => "$cst_id\-$ref_id");
			my $seqobj;

			if($strand == '-1'){
				$seqobj = $seq->revcom;
			}
			else{
				$seqobj = $seq;
			}
			$fasta->write_seq($seqobj);
		}
	}	
	return("$type\_cst_$cst\.fa")
}


=head2 get_chr

 Title   : get_chr

 Usage   : my $length = Bio::MCE::Utils->get_chr($chr,$specie);

 Function: a slice of the entire chromosome

 Returns : a slice object

 Args    : -1 chr
           -2 specie

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!
					 TO CHANGE: YOU HAVE TO PASS ALSO THE COORD_SYS_NAME
=cut

sub get_chr {

  my $class = shift;
  my $chr = shift;
  my $specie = shift;

  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice") or die "ERROR 4 $!\n";

  my $slice = $slice_adaptor->fetch_by_region('',$chr) or die "ERROR 5 $!\n";

  return($slice);
}


=head2 get_locals

 Title   : get_locals

 Usage   : my $pw = Bio::MCE::Utils->get_locals($feature_id);

 Function: retrieve all the promoterwise or chaos hits from a vista_cst

 Returns : arrayref of hashref from the table feature_pair

 Args    : -1 database handler
           -2 feature_id

=cut

sub get_locals {

  my $class = shift;
	my $dbh = shift;
	my $fid = shift;
	my @fp;

	my $sth = $dbh->prepare('SELECT * FROM feature_pair '.
													' WHERE feature_id = '.$dbh->quote($fid).
													' AND primary_tag = '.$dbh->quote('pw'));
	$sth->execute;

	while(my $row = $sth->fetchrow_hashref) {

		push(@fp,$row);
	}

	return(\@fp);
}


=head2 get_pos_from_string

 Title   : get_pos_from_string

 Usage   : my $pos = Bio::MCE::Utils->get_pos_from_string($seq,$string);

 Function: return the position at which START a string into a sequence

 Returns : an arrayref of digit

 Args    : -1 the seq object on which to search
           -2 string to search

=cut

sub get_pos_from_string {

  my $class = shift;
  my $seqobj = shift;
	my $string = shift;
	my @pos;

	my $seq = $seqobj->seq;
	
	while($seq =~ m/($string)/gi){
		#my $pos = (pos($seq)-length($string)+1);
		my $pos = ($-[1]+1);
		push(@pos,$pos);
	} 
  return(\@pos);
}


=head2 get_seq_from_feature

 Title   : get_seq_from_feature

 Usage   : my $pos = Bio::MCE::Utils->get_seq_from_feature($dbh,$id);

 Function: return the seqobj from a feature

 Returns : seqobj

 Args    : -1 the dbh
           -2 the id of the feature

=cut

sub get_seq_from_feature {

  my $class = shift;
  my $dbh = shift;
  my $id = shift;

	my $sth = $dbh->prepare('SELECT s.ref_id, s.sequence, start, end '.
													'FROM seq s, feature f '.
													'WHERE feature_id = '.$dbh->quote($id).
													' AND s.ref_id=f.ref_id');
	$sth->execute;

	my $row = $sth->fetchrow_hashref;

	my $offset = $row->{'start'}-1;
	my $length = $row->{'end'}-$row->{'start'}+1;
	
	my $seq = Bio::Seq->new(-id => $row->{'ref_id'},
													-seq => substr($row->{'sequence'},$offset,$length));

	return $seq;
}


=head2 get_seq_from_feature_pair

 Title   : get_seq_from_feature_pair

 Usage   : my $pos = Bio::MCE::Utils->get_seq_from_feature_pair($dbh,$id);

 Function: return the seqobj from a feature_pair

 Returns : seqobj

 Args    : -1 the seqobj
           -2 the id of the feature_pair

=cut

sub get_seq_from_feature_pair {

  my $class = shift;
  my $dbh = shift;
  my $id = shift;

  my $sth = $dbh->prepare('SELECT s.ref_id, s.sequence, q_start, q_end '.
                          'FROM seq s, feature_pair fp '.
                          'WHERE feature_pair_id = '.$dbh->quote($id).
                          ' AND s.ref_id=fp.q_ref_id');
  $sth->execute;

  my $row = $sth->fetchrow_hashref;

  my $offset = $row->{'q_start'}-1;
  my $length = $row->{'q_end'}-$row->{'q_start'}+1;

  my $seq = Bio::Seq->new(-id => $id,
                          -seq => substr($row->{'sequence'},$offset,$length));

  return $seq;
}


=head2 get_target_seq_from_feature_pair

 Title   : get_target_seq_from_feature_pair

 Usage   : my $pos = Bio::MCE::Utils->get_target_seq_from_feature_pair($dbh,$id);

 Function: return the seqobj from a feature_pair

 Returns : seqobj

 Args    : -1 the seqobj
           -2 the id of the feature_pair

=cut

sub get_target_seq_from_feature_pair {

  my $class = shift;
  my $dbh = shift;
  my $id = shift;

  my $sth = $dbh->prepare('SELECT s.ref_id, s.sequence, t_start, t_end '.
                          'FROM seq s, feature_pair fp '.
                          'WHERE feature_pair_id = '.$dbh->quote($id).
                          ' AND s.ref_id=fp.t_ref_id');
  $sth->execute;

  my $row = $sth->fetchrow_hashref;

  my $offset = $row->{'t_start'}-1;
  my $length = $row->{'t_end'}-$row->{'t_start'}+1;

  my $seq = Bio::Seq->new(-id => $id,
                          -seq => substr($row->{'sequence'},$offset,$length));

  return $seq;
}


=head2 set_lagan_tree

 Title   : set_lagan_tree

 Usage   : my $tree = Bio::MCE::Utils->set_lagan_tree($hashref);

 Function: write the string with the tree for running mlagan

 Returns : a string

 Args    : -1 an hashref with the genome as key

 Note    : TO UPDATE EACH TIME ACCORDING TO THE USED GENOMES!

=cut

sub set_lagan_tree {

  my $class = shift;
  my $work = shift;
  my @strings;

  my $tree;
  @{$tree->{'rodents'}} = ('Mus musculus', 'Rattus norvegicus', 'Monodelphis domestica');
  @{$tree->{'bipeds'}} = ('Homo sapiens', 'Pan troglodytes', 'Macaca mulatta');
  @{$tree->{'others'}} = ('Canis familiaris', 'Bos taurus');
  @{$tree->{'birds'}} = ('Gallus gallus');
  @{$tree->{'fishes'}} = ('Takifugu rubripes', 'Tetraodon nigroviridis', 'Oryzias latipes', 'Gasterosteus aculeatus',
                           'Danio rerio');
  @{$tree->{'anphs'}} = ('Xenopus tropicalis');
  @{$tree->{'cionas'}} = ('Ciona intestinalis', 'Ciona savignyi');

  my @order = qw(rodents bipeds others birds fishes anphs cionas);

  foreach my $group(@order) {

    my @got = keys(%$work);
    my $string;
    my $count = 0;

    foreach my $s(@{$tree->{$group}}) {
      if(grep(/$s/,@got)){
        $count++;
      }
    }

    my $c = 0;
    next unless $count;

    foreach my $s(@{$tree->{$group}}) {
      if(grep(/$s/,@got)){
        $c++;
        if($count == 1) {
          $string .= $work->{$s};
        }
        elsif($count > 1 and $c == 1) {
          $string .= '('.$work->{$s};
        }
        elsif($count > 1 and $c > 1) {
          $string .= $work->{$s}.')';
        }
        $string .= ' ';
      }
    }
    my $bra_to_add = '(' x ($count-2);
    $string = $bra_to_add.$string;
    $string =~ s/ $//g if $count > 1;
    $string =~ s/ $//g if $string =~ /\) +$/;
    push(@strings,$string);
  }
  _finalize(\@strings);
}


=head2 _finalize

 Title   : _finalize

 Usage   : 

 Function: Internal used by get_lagan_tree

 Returns : A string: the lagan tree

 Args    : 

=cut

sub _finalize {

  my $strings = shift;
  my $c = 0;
  my $tree;
  foreach my $string(@$strings) {
    $c++;
    if($c == 1) {
      $tree = $string;
    }
    else {
      $tree = '('.$tree.$string.')';
    }
  }
	$tree =~ s/ +\)/\)/g;
	$tree =~ s/ +\(/\(/g;
	$tree =~ s/\) +/\)/g;
	$tree =~ s/\( +/\(/g;
  return $tree;
}


=head2 get_coord_intergenic

 Title   : get_coord_intergenic

 Usage   : my($coord_sys_name,seq_region_name,$start,$end,$strand,$name) = 
           Bio::MCE::Utils->get_coord_intergenic($stable_id,$slice_adaptor);

 Function: retrieve the coordinates of the slice containing the gene plus
           flanking up- and down-stream till the next not-overlapping gene.
           The strand is returned in the way that the gene is forward on this
           slice

 Returns : coordinate system name,
           seq region name,
           start,
           end,
           strand
           name of the ref gene

 Args    : -1 gene_stable_id
           -2 slice adaptor
           -3 the length of the flanking region or 'intergenic'

=cut

sub get_coord_intergenic {

	my $self = shift;
  my $gene_id = shift;
  my $slice_adaptor = shift;
	my $flanking = shift;
	
	my $pre_slice = $slice_adaptor->fetch_by_gene_stable_id($gene_id);
	my $coord_sys_name = $pre_slice->coord_system->name;
	my $seq_region_name = $pre_slice->seq_region_name;
	my $slice = $slice_adaptor->fetch_by_region($coord_sys_name,$seq_region_name);

	my $gene;
  foreach my $slice_gene(@{$slice->get_all_Genes}) {
    $gene = $slice_gene if $slice_gene->stable_id eq $gene_id;
  }

	my $name = $gene->external_name;

  my $start = $slice->start;
  my $end = $slice->end;

  ($start,$end) = _adjust_start_end($start,$end,$slice,$gene_id,$gene);
  
	my $strand = _get_gene_strand_on_slice($slice_adaptor,$gene_id,$coord_sys_name,$seq_region_name,$start,$end);
	
  return($coord_sys_name,$seq_region_name,$start,$end,$strand,$name) if $flanking eq 'intergenic';
  return($coord_sys_name,$seq_region_name,($gene->start-$flanking+1),($gene->end+$flanking-1),$strand,$name);	
}


=head2 _adjust_start_end

 Title   : _adjust_start_end

 Usage   : 

 Function: Internal

 Returns : 

 Args    : 

=cut

sub _adjust_start_end {

  # LE COORDINATE VENGONO AGGIUSTATE SOLO SE IL GENE
  # CONFINANTE E' TOTALMENTE ESTERNO A QUELLO CONSIDERATO...
  # SFRUTTANDO IL FATTO CHE LO START E' SEMPRE MINORE DELL'END

  my $start = shift;
  my $end = shift;
  my $slice = shift;
  my $gene_id = shift;
  my $ref = shift;
  my @coord_start;
  my @coord_end;
  
  my @genes = @{$slice->get_all_Genes};

  foreach my $gene(@genes) {

    if($gene->stable_id ne $gene_id) {
      push(@coord_start,($slice->start+$gene->start-1));
      push(@coord_end,($slice->start+$gene->end-1));
    }
  }
 
  my $new_start = $start;
  foreach my $coord(sort{$a<=>$b}@coord_end) {
    $new_start = $coord if ($coord > $start && $coord < $start+$ref->start-1);
  }
  
  my $new_end = $end;
  foreach my $coord(sort{$b<=>$a}@coord_start) {
    $new_end = $coord if ($coord < $end && $coord > $start+$ref->end-1);
  }

  # TO AVOID GENES OF 1 BP ON THE EDGE OF THE SLICE:
  $new_start += 1;
  $new_end -= 1;
   
  return($new_start,$new_end);
}


=head2 _get_gene_strand_on_slice

 Title   : _get_gene_strand_on_slice

 Usage   :

 Function: Internal

 Returns :

 Args    :

=cut

sub _get_gene_strand_on_slice {
	
	my $slice_adaptor = shift;
	my $gene_id = shift;
	my $coord_sys_name = shift;
	my $seq_region_name = shift;
	my  $start = shift;
	my $end = shift;
	my $strand;
	
	my $slice = $slice_adaptor->fetch_by_region($coord_sys_name,$seq_region_name,$start,$end);
	
	my @genes = @{$slice->get_all_Genes};

  foreach my $gene(@genes) {

    if($gene->stable_id eq $gene_id) {
			
			$strand = $gene->strand;
		}
	}
	
	return $strand;
}


=head2 mask_seq

 Title   : mask_seq

 Usage   : my $masked_seq = Bio::MCE::Utils->mask_seq($seq,$start,$end);

 Function: mask a sequence given the start, end coordinates of the 
           fragment to mask using 'X'

 Returns : a string

 Args    : -1 a sequence (string)
           -2 start
           -3 end

=cut

sub mask_seq {
	
	my $class = shift;
	my $seq = shift;
	my $start = shift;
	my $end = shift;
	
	my $length = ($end - $start + 1);
	
	substr($seq,($start-1),$length,'X' x $length);

	return $seq;
}


=head2 mask_features

 Title   : mask_seq

 Usage   : my $masked_seq = Bio::MCE::Utils->mask_features($seq,$features);

 Function: replaces string positions described in the Features with 'N'

 Returns : a string

 Args    : -1 a sequence (string)
           -2 array_ref of features

=cut

sub mask_features {

  my $self = shift;
	my $seq = shift;
	my $features = shift;
	
  # explicit CORE::length call, to avoid any confusion with the Slice
  # length method
  my $dnalen = CORE::length($seq);

  REP:foreach my $f (@{$features}) {
    my $start  = $f->start;
    my $end    = $f->end;
    my $length = ($end - $start) + 1;

    # check if we get feature completely outside of expected slice range
    if ($end < 1 || $start > $dnalen) {
       print STDERR "Unexpected: Feature completely outside slice coordinates";
      next REP;
    }

    # feature partly outside slice range, so correct
    # the repeat start and length to the slice size if needed
    if ($start < 1) { 
      $start = 1;
      $length = ($end - $start) + 1;
    }

    # feature partly outside slice range, so correct
    # the repeat end and length to the slice size if needed
    if ($end > $dnalen) {
      $end = $dnalen;
      $length = ($end - $start) + 1;
    }

    $start--;

    my $padstr;
    $padstr = 'N' x $length;
    substr ($seq,$start,$length) = $padstr;
  }
	return $seq;
}


=head2 get_bioseq_with_coding_features

 Title   : get_bioseq_with_coding_features

 Usage   : my $ann_slice = Bio::MCE::Utils->get_bioseq_with_coding_features($slice,
                                                                            $gene_id,1,1);

 Function: add to the slice the feature of 'gene', 'mRNA', 'exon', 'utr'

 Returns : a seq object with features

 Args    : -1 a slice
           -2 stable_id of the ref gene
           -3 '1' if you want to mask repeats
           -4 '1' if you want to mask coding exons

=cut

sub get_bioseq_with_coding_features {
	
	my $class = shift;
	my $slice = shift;
	my $id = shift;
	my $mask_repeat = shift;
	my $mask_coding_exon = shift;

	my $seq;	
	if($mask_repeat){
		$seq = $slice->get_repeatmasked_seq->seq;
	}
	else{
		$seq = $slice->seq();
	}
	my $seq_obj = Bio::Seq->new(-id => $id);

	foreach my $slice_gene (@{$slice->get_all_Genes()}) {
		
		my $transcript_count = 0;
				
		my $feat_name = ($slice_gene->external_name || $slice_gene->stable_id);
    
		my $feat = Bio::SeqFeature::Generic->new(-seq_id => $slice_gene->stable_id,
																						 -start => $slice_gene->start,
																						 -end => $slice_gene->end,
																						 -strand => $slice_gene->strand,
																						 -primary_tag => 'gene',
																						 -source => $slice_gene->stable_id,
																						 -tag => {'Gene' => $feat_name});
    
		$seq_obj->add_SeqFeature($feat);
    
		foreach my $slice_transcript(@{$slice_gene->get_all_Transcripts}) {
			
			$transcript_count++;   
						
			my $featT = Bio::SeqFeature::Generic->new(-seq_id => $slice_gene->stable_id,
																								-start => $slice_transcript->start,
																								-end => $slice_transcript->end,
																								-strand => $slice_transcript->strand,
																								-primary_tag => 'mRNA',
																								-source => $slice_gene->stable_id,
																								-tag => {'mRNA' => ($feat_name.'.'.$transcript_count),
																												 'transcript' => $slice_transcript->stable_id});
			
			$seq_obj->add_SeqFeature($featT);
			
			# now work on translated regions but first test if transcript is translated!!!
			if ($slice_transcript->translation) {
				
				my($translated_regions,$utr_regions) = Bio::MCE::Utils->get_exon_utr($slice_transcript);

				foreach my $exon(@{$translated_regions}) {
					
					my($start,$end) = Bio::MCE::Utils->limit_to_slice($slice,$exon);
					my $length = ($end - $start + 1);
					next if($length <= 1);
					
					my $featE = Bio::SeqFeature::Generic->new(-seq_id => $slice_gene->stable_id,
																										-start => $start, 
																										-end => $end,
																										-strand => $slice_gene->strand,
																										-primary_tag => 'exon',
																										-source => $slice_gene->stable_id,
																										-tag => {'mRNA' => ($feat_name.'.'.$transcript_count),
																														 'transcript' => $slice_transcript->stable_id});
					
					$seq_obj->add_SeqFeature($featE);

					substr($seq,($start-1),$length,'N' x $length) if $mask_coding_exon;
				}

				foreach my $utr(@{$utr_regions}) {

					my($start,$end) = Bio::MCE::Utils->limit_to_slice($slice,$utr);
					my $length = ($end - $start + 1);
					next if ($length <= 1);

					my $featU = Bio::SeqFeature::Generic->new(-seq_id => $slice_gene->stable_id,
																										-start => $start, 
																										-end => $end,
																										-strand => $slice_gene->strand,
																										-primary_tag => 'utr',
																										-source => $slice_gene->stable_id,
																										-tag => {'mRNA' => ($feat_name.'.'.$transcript_count),
																														 'transcript' => $slice_transcript->stable_id});
					
					$seq_obj->add_SeqFeature($featU);
				}
			}

      else {

				my $non_coding_exons = $slice_transcript->get_all_Exons;
        foreach my $exon(@$non_coding_exons) {

          my($start,$end) = Bio::MCE::Utils->limit_to_slice($slice,$exon);
          my $length = ($end - $start + 1);
          next if($length <= 1);

          my $featE = Bio::SeqFeature::Generic->new(-seq_id => $slice_gene->stable_id,
                                                    -start => $start, 
                                                    -end => $end,
                                                    -strand => $slice_gene->strand,
                                                    -primary_tag => 'exon',
                                                    -source => $slice_gene->stable_id,
                                                    -tag => {'mRNA' => ($feat_name.'.'.$transcript_count),
                                                             'transcript' => $slice_transcript->stable_id});

          $seq_obj->add_SeqFeature($featE);

          substr($seq,($start-1),$length,'N' x $length) if $mask_coding_exon;
				}
			}
		}
	}
	$seq_obj->seq($seq);
	return $seq_obj;
}


=head2 get_seq_by_id_ses

 Title   : get_seq_by_id_ses

 Usage   : my $seq = Bio::MCE::Utils->get_seq_by_id_ses($dbh,$ref_id,$start,$end,$strand);
           print $seq."\n";

 Function: Fetch a sequence from seq table

 Returns : string

 Args    : 1 - database handler
           2 - ref_id
           3 - start
           4 - end

=cut

sub get_seq_by_id_ses {

  my $self = shift;
  my $dbh = shift;
  my $id = shift;
  my $start =shift;
  my $end = shift;
	my $strand = shift;
	my $seqobj;
		
  my $off = ($start-1);
  my $length = ($end-$start+1);

  my $sth = $dbh->prepare('SELECT * FROM seq '.
  									      'WHERE ref_id = '.$dbh->quote($id));
  $sth->execute;

  my $s_row = $sth->fetchrow_hashref;

  my $seq = substr($s_row->{'sequence'},$off,$length);

	$seq =~ s/X/N/g;
	$seq =~ s/x/n/g;

	unless($seq) {
		print "\nSeq ".$id.':'.$start.'-'.$end." is EMPTY!\n";
		my $not = Bio::Seq->new(-id =>  $id.':'.$start.'-'.$end, 
                            -seq => 'undef');
		return $not;
	}

	my $preseqobj = Bio::Seq->new(-id => $id.':'.$start.'-'.$end,
															  -seq => $seq);
	if($strand == -1) {
		$seqobj = $preseqobj->revcom;
	}
	else {
		$seqobj = $preseqobj;
	}
  if($length>length($s_row->{'sequence'})) {
		print "\nLength of subseq greater than sequence $id in db!\n";
		$seqobj = Bio::Seq->new(-id => $id.':'.$start.'-'.$end,  
                            -seq => 'undef');
	}
  return $seqobj;
}


=head2 generate_random_string

 Title   : generate_random_string

 Usage   : my $id = Bio::MCE::Utils->generate_random_string('10');

 Function: generrate a random string of the required length

 Returns : a string

 Args    : -1 length of the string needed

 Note    :

=cut

sub generate_random_string {

	my $class = shift;
  my $length_of_randomstring = shift;

  my @chars=('a'..'z','A'..'Z','_');
  my $random_string;
  foreach (1..$length_of_randomstring) {
    # rand @chars will generate a random 
    # number between 0 and scalar @chars
    $random_string .= $chars[rand @chars];
  }
  return $random_string;
}


=head2 longest_translation

 Title   : longest_translation

 Usage   : my $transcript = Bio::MCE::Utils->longest_translation($id,$specie);

 Function: retrieve the transcript object belonging to the longest translation
           from the required gene in the required genome

 Returns : a Bio::EnsEMBL::Transcript object

 Args    : -1 gene_stable_id
           -2 specie

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub longest_translation {

  my $class = shift;
  my $id = shift;
  my $specie = shift;

  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","gene");

  my $gene = $gene_adaptor->fetch_by_stable_id($id);

	my $transcripts = $gene->get_all_Transcripts;
	my %transcript_hash;
	foreach my $transcript(@{$transcripts}) {
		# WE NEED THE LONGEST TRANSLATION IN ORDER TO BE COMPARA COMPATIBLE!
		$transcript_hash{$transcript->translate->length} = $transcript if $transcript->translation;
	}
	my @length = keys %transcript_hash;
	my @descending_length = sort {$b<=>$a} @length;
	my $longest = $transcript_hash{$descending_length[0]};
 
	return $longest;
}


=head2 get_genomic_coord

 Title   : get_genomic_coord

 Usage   : my $start = Bio::MCE::Utils->get_genomic_coord($dbh,$ref_id,$pre_start,$pre_end,$job_id);

 Function: give the genomic coord from a feature on a ref_id

 Returns : integer

 Args    : -1 dbh of db with 'seq' table containing ref_id and genomic location
              of extracted sequence
           -2 ref_id
           -3 start coord
           -4 end coord
           -5 job_id or '' (empty)

 Note    : 

=cut

sub get_genomic_coord {

	my $class = shift;
	my $dbh = shift;
	my $ref_id = shift;
	my $pre_start = shift;
	my $pre_end = shift;
	my $job_id = shift;

	my $string = 'SELECT chr_start, chr_end, slice_strand FROM seq WHERE ref_id = ';
	$string .= "\'$ref_id\'";
	$string .= " AND job_id \= \'$job_id\'" if $job_id;

	my $sth = $dbh->prepare($string);

	$sth->execute;

	my $row = $sth->fetchrow_hashref;

	my $chr_start = $row->{'chr_start'};
	my $chr_end = $row->{'chr_end'};
	
	my $start;
	my $end;

	if($row->{'slice_strand'} == 1) {
		$start = $chr_start + $pre_start - 1;
		$end = $chr_start + $pre_end - 1;
	}
	else {
		$start = $chr_end - $pre_end + 1;
		$end = $chr_end - $pre_start + 1;
	}

	return($start,$end);
}


=head2 get_nearest_gene_from_region

 Title   : get_nearest_gene_from_region

 Usage   : my $gene_stable_id = Bio::MCE::Utils->get_nearest_gene_from_region();

 Function: 

 Returns : the gene stable id of the gene nearest to the supplied coord in term of start/end
           of the entire gene. DO NOT TAKE EXON BOUNDUARIES IN CONSIDERATION!!!

 Args    : -1 the name of the specie
           -2 coord system name
           -3 seq region name
           -4 the start of a feature on the same seq region
           -5 the end of a feature on the same seq region

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!
           DO NOT TAKE EXON BOUNDUARIES IN CONSIDERATION!!!
           THIS FUNCTION JUST ASSOCIATE THE FEATURE TO THE
           START-END OF THE GENES

=cut

sub get_nearest_gene_from_region {

  my $class = shift;
	my $specie = shift;
  my $coord_system_name = shift;
  my $seq_region_name = shift;
	my $seq_region_start = shift;
	my $seq_region_end = shift;

  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice");

  my $slice = $slice_adaptor->fetch_by_region($coord_system_name,$seq_region_name);

	my $slice_genes = $slice->get_all_Genes();	

	my $smaller_dist;
	my $candidate;

	foreach my $gene(@$slice_genes) {
		my $start = $gene->seq_region_start;
		my $end = $gene->seq_region_end;
		my $distance = _get_smaller_distance($seq_region_start,$seq_region_end,$start,$end);

		unless($candidate) {
			$candidate = $gene->stable_id;
			$smaller_dist = $distance;
		}

		if($distance < $smaller_dist) {
			$candidate = $gene->stable_id;
			$smaller_dist = $distance;
		}
	}
	return($candidate);
}

=head2 _get_smaller_distance

 Title   : _get_smaller_distance

 Usage   : It take two pair in this order: start1, end1 and start2, end2
           of coord and return the minimum distance between the two segments
           even if they overlap in fact it sorts the coord before of all

 Function: Internal used by get_nearest_from_region

 Returns : THE MINIMUM DISTANCE BETWEEN THE TWO SETS OF COORDINATES

 Args    : -1 start1
           -2 end1
           -3 start2
           -4 end2

 NOTE    : THE START HAVE TO BE ALWAYS SMALLER THAN END

=cut

sub _get_smaller_distance {

	my $s1 = shift;
	my $e1 = shift;
	my $s2 = shift;
	my $e2 = shift;

	next unless($s1 && $e1 && $s2 && $e2);

	my @pair;
	my @dis;

  @{$pair[0]} =  sort {$b <=> $a} ($s1,$s2);
  @{$pair[1]} =  sort {$b <=> $a} ($s1,$e2);
  @{$pair[2]} =  sort {$b <=> $a} ($e1,$s2);
  @{$pair[3]} =  sort {$b <=> $a} ($e1,$e2);

  foreach my $edge(@pair) {
    push(@dis,($$edge[0]-$$edge[1]));
  }

  my @last = sort {$a <=> $b} @dis;

  return(shift(@last));
}


=head2 get_nearest_external_left_gene

 Title   : get_nearest_external_left_gene

 Usage   : my $gene_stable_id = Bio::MCE::Utils->get_nearest_external_left_gene(.......);

 Function: 

 Returns : the gene stable id of the gene nearest to the supplied coord in term of start/end
           of the entire gene. DO NOT TAKE EXON BOUNDUARIES IN CONSIDERATION!!!

           IN RELATION TO THE LEFT PART OF YOUR COORDINATES

           IF THE FEATURE IS INTRONIC THE GENE IN WHICH IT IS INTRONIC
           IS NOT TAKEN IN CONSIDERATION

 Args    : -1 the name of the specie
           -2 coord system name
           -3 seq region name
           -4 the start of a feature on the same seq region
           -5 the end of a feature on the same seq region

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!
           DO NOT TAKE EXON BOUNDUARIES IN CONSIDERATION!!!
           THIS FUNCTION JUST ASSOCIATE THE FEATURE TO THE
           START-END OF THE GENES

=cut

sub get_nearest_external_left_gene {

  my $class = shift;
	my $specie = shift;
  my $coord_system_name = shift;
  my $seq_region_name = shift;
	my $seq_region_start = shift;
	my $seq_region_end = shift;

  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice");

  my $slice = $slice_adaptor->fetch_by_region($coord_system_name,$seq_region_name);

	my $slice_genes = $slice->get_all_Genes();	

	my $smaller_dist;
	my $candidate;

	foreach my $gene(@$slice_genes) {

		next if $gene->seq_region_start > $seq_region_start;

		my $start = $gene->seq_region_start;
		my $end = $gene->seq_region_end;
		my $distance = _get_smaller_distance($seq_region_start,$seq_region_end,$start,$end);

		unless($candidate) {
			$candidate = $gene->stable_id;
			$smaller_dist = $distance;
		}

		if($distance < $smaller_dist) {
			my $region = Bio::MCE::Utils->relate_feature_to_gene($coord_system_name,$seq_region_name,$seq_region_start,$seq_region_end,$gene->stable_id);

			next if $region =~ /int/;

			$candidate = $gene->stable_id;
			$smaller_dist = $distance;
		}
	}
	return($candidate);
}


=head2 get_nearest_external_right_gene

 Title   : get_nearest_external_right_gene

 Usage   : my $gene_stable_id = Bio::MCE::Utils->get_nearest_external_right_gene(.......);

 Function: 

 Returns : the gene stable id of the gene nearest to the supplied coord in term of start/end
           of the entire gene. DO NOT TAKE EXON BOUNDUARIES IN CONSIDERATION!!!

           IN RELATION TO THE RIGHT PART OF YOUR COORDINATES

           IF THE FEATURE IS INTRONIC THE GENE IN WHICH IT IS INTRONIC
           IS NOT TAKEN IN CONSIDERATION

 Args    : -1 the name of the specie
           -2 coord system name
           -3 seq region name
           -4 the start of a feature on the same seq region
           -5 the end of a feature on the same seq region

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!
           DO NOT TAKE EXON BOUNDUARIES IN CONSIDERATION!!!
           THIS FUNCTION JUST ASSOCIATE THE FEATURE TO THE
           START-END OF THE GENES

=cut

sub get_nearest_external_right_gene {

  my $class = shift;
	my $specie = shift;
  my $coord_system_name = shift;
  my $seq_region_name = shift;
	my $seq_region_start = shift;
	my $seq_region_end = shift;

  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice");

  my $slice = $slice_adaptor->fetch_by_region($coord_system_name,$seq_region_name);

	my $slice_genes = $slice->get_all_Genes();	

	my $smaller_dist;
	my $candidate;

	foreach my $gene(@$slice_genes) {

		next if $gene->seq_region_start < $seq_region_start;

		my $start = $gene->seq_region_start;
		my $end = $gene->seq_region_end;
		my $distance = _get_smaller_distance($seq_region_start,$seq_region_end,$start,$end);

		unless($candidate) {
			$candidate = $gene->stable_id;
			$smaller_dist = $distance;
		}

		if($distance < $smaller_dist) {
			my $region = Bio::MCE::Utils->relate_feature_to_gene($coord_system_name,$seq_region_name,$seq_region_start,$seq_region_end,$gene->stable_id);

			next if $region =~ /int/;

			$candidate = $gene->stable_id;
			$smaller_dist = $distance;
		}
	}
	return($candidate);
}


=head2 is_esonic

 Title   : is_esonic

 Usage   : my $gene_stable_id = Bio::MCE::Utils->is_esonic(.......);

 Function: 

 Returns : the gene stable id of the gene in which your coordinated are located
           or 'undef'

 Args    : -1 the name of the specie
           -2 coord system name
           -3 seq region name
           -4 the start of a feature on the same seq region
           -5 the end of a feature on the same seq region

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub is_esonic {

  my $class = shift;
	my $specie = shift;
  my $coord_system_name = shift;
  my $seq_region_name = shift;
	my $seq_region_start = shift;
	my $seq_region_end = shift;

  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice");

  my $slice = $slice_adaptor->fetch_by_region($coord_system_name,$seq_region_name,$seq_region_start,$seq_region_end);

	my $slice_genes = $slice->get_all_Genes();	

	my $candidate;

	foreach my $gene(@$slice_genes) {
		my $start = $gene->seq_region_start;
		my $end = $gene->seq_region_end;
		my $region = Bio::MCE::Utils->relate_feature_to_gene($coord_system_name,$seq_region_name,$seq_region_start,$seq_region_end,$gene->stable_id);
		next unless $region =~ /exon/;
		$candidate = $gene->stable_id;
	}
	return($candidate || undef);
}


=head2 is_intronic

 Title   : is_intronic

 Usage   : my $gene_stable_id = Bio::MCE::Utils->is_intronic(.......);

 Function: 

 Returns : the gene stable id of the gene in which your coordinated are located
           or 'undef'

 Args    : -1 the name of the specie
           -2 coord system name
           -3 seq region name
           -4 the start of a feature on the same seq region
           -5 the end of a feature on the same seq region

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub is_intronic {

  my $class = shift;
	my $specie = shift;
  my $coord_system_name = shift;
  my $seq_region_name = shift;
	my $seq_region_start = shift;
	my $seq_region_end = shift;

  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice");

  my $slice = $slice_adaptor->fetch_by_region($coord_system_name,$seq_region_name,$seq_region_start,$seq_region_end);

	my $slice_genes = $slice->get_all_Genes();	

	my $smaller_dist;
	my $candidate;

	foreach my $gene(@$slice_genes) {

		my $start = $gene->seq_region_start;
		my $end = $gene->seq_region_end;
		my $distance = _get_smaller_distance($seq_region_start,$seq_region_end,$start,$end);

		unless($candidate) {
			$candidate = $gene->stable_id;
			$smaller_dist = $distance;
		}

		if($distance < $smaller_dist) {
			my $region = Bio::MCE::Utils->relate_feature_to_gene($coord_system_name,$seq_region_name,$seq_region_start,$seq_region_end,$gene->stable_id);

			next unless $region =~ /int/;

			$candidate = $gene->stable_id;
			$smaller_dist = $distance;
		}
	}
	return($candidate || undef);
}


=head2 relate_feature_to_gene

 Title   : relate_feature_to_gene

 Usage   : my $ergion = Bio::MCE::Utils->relate_feature_to_gene($coord_sys_name,$seq_region_name,$chr_start,$chr_end,$gene_id);

 Function: take 2 genomic coordinates and the name of the gene to which
           to relate this coordinate and return the relative position of the
					 feature in respect to the  gene

 Returns : a string

 Args    : -1 coordinate system name
           -2 seq region name
           -3 chr start
           -4 chr end
					 -5 gene stable id

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub relate_feature_to_gene {

  my $class = shift;
	my $coord_sys_name = shift;
	my $seq_region_name = shift;
	my $start = shift;
	my $end = shift;
  my $gene_id = shift;

  my $specie = Bio::MCE::Utils->get_genome_from_id($gene_id);

  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","gene") or die "ERROR $!\n";
    
  my $gene = $gene_adaptor->fetch_by_stable_id($gene_id) or die "ERROR $!\n";

#	if($gene->coord_system_name ne $coord_sys_name || $gene->seq_region_name ne $seq_region_name) {
	if($gene->seq_region_name ne $seq_region_name) {

#####
print $gene->coord_system_name." doveva essere $coord_sys_name ed ".$gene->seq_region_name." doveva essere $seq_region_name \n";
####
		$gene->throw('mapping not possible');
		return undef;
	}

	my $transcript;
	my $map_collection;

	if($gene->biotype eq 'protein_coding') {
	  $transcript = Bio::MCE::Utils->longest_translation($gene_id,$specie);
		$map_collection = Bio::MCE::Pipeline::Coding->map_transcript($transcript,$specie);
	}
		
	else {
		$map_collection = Bio::MCE::Pipeline::Coding->map_non_protein_coding($gene,$specie);
	}

	my $range = Bio::Range->new(-start => $start,
															-end => $end);

	my @over_zone;
	foreach my $ele($map_collection->get_all_features) {
		push (@over_zone,($ele->primary_tag)) if $range->overlaps($ele);
	}

	my $zone = join(', ',@over_zone);	    
	$zone =~ s/, $//;
	return($zone);
}

=head2 _regex_id_group

 Title   : _regex_id_group

 Usage   : my $group = _regex_id_group();

 Function: it is a function in order to have a regular expression to
           link a gene stable id to a group !!!INTERNAL METHOD!!!!

           ACTUAL HASH:

             my %hash_id = ('ENSG0' => 'Mammal',
                            'ENSMUSG' => 'Mammal',
                            'ENSRNOG' => 'Mammal',
               #             'NEWSINFRUG' => 'Fish',
#						                'SINFRUG' => 'Fish',
														 'ENSTRUG' => 'Fish',
                            'ENSDARG' => 'Fish',
                            'ENSPTRG' => 'Mammal',
                            'ENSGALG' => 'Bird',
                            'ENSCAFG' => 'Mammal',
                            'ENSCING' => 'Tunicate',
                            'ENSCSAVG0' => 'Tunicate',
                            'ENSXETG' => 'Anfibia',
                            'AY' => 'Fish',
                            'EVX' => 'Fish',
                            'GSTENG' => 'Fish',
                            'HOX' => 'Fish',
                            'ENSBTAG' => 'Mammal',
                            'ENSMODG' => 'Mammal',
                            'ENSMMUG' => 'Mammal',
                            'ENSGACG' => 'Fish',
                            'ENSORLG0' => 'Fish');

 Returns : hasref of id_regex => group

 Args    : no

 NOTE    : INTERNAL METHOD
					 TO BE UPDATED WHEN ADD OR DELETE SPECIES FROM ENSEMBL

=cut

sub _regex_id_group {

	my %hash_id = ('ENSG0' => 'Mammal',
								 'ENSMUSG' => 'Mammal',
								 'ENSRNOG' => 'Mammal',
				#				 'NEWSINFRUG' => 'Fish',
#                 'SINFRUG' => 'Fish',
									'ENSTRUG' => 'Fish',
								 'ENSDARG' => 'Fish',
								 'ENSPTRG' => 'Mammal',
								 'ENSGALG' => 'Bird', 
								 'ENSCAFG' => 'Mammal',
								 'ENSCING' => 'Tunicate',
								 'ENSCSAVG0' => 'Tunicate',
								 'ENSXETG' => 'Anfibia',
								 'AY' => 'Fish',
								 'EVX' => 'Fish',
								 'GSTENG' => 'Fish',
								 'HOX' => 'Fish',
								 'ENSBTAG' => 'Mammal',
								 'ENSMODG' => 'Mammal',
								 'ENSMMUG' => 'Mammal',
								 'ENSGACG' => 'Fish',
                 'ENSORLG0' => 'Fish');

	return(\%hash_id);
}


=head2 regex_id_group

 Title   : regex_id_group

 Usage   : my $id_genome = Bio::MCE::Utils->regex_id_group();

 Function: it is a function in order to have a regular expression to
           link a gene stable id to a group

 Returns : hasref of id_regex => group

 Args    : no

=cut

sub regex_id_group {

	my $class = shift;
	
	my $hash_id = _regex_id_group();

	return($hash_id);
}

	
=head2 get_group_from_id

 Title   : get_group_from_id

 Usage   : my $id_genome = Bio::MCE::Utils->get_group_from_id($id);

 Function: Returns the group from the suffix of the stable_id

 Returns : a string

 Args    : ensembl stable gene id 

=cut

sub get_group_from_id {

	my $class = shift;
	my $id = shift;
	my $c = 0;

	my $hash_id = _regex_id_group();

	foreach my $key(keys %$hash_id) {
		if($id =~ /$key/) {

			return($hash_id->{$key});
			my $c ++;
		}
	}
#	print STDERR "\n$id not in my book!\n" and return undef unless $c;
}


=head2 relative_strand

 Title   : relative_strand

 Usage   : Bio::PMCE::Utils->relative_strand($stable_id1,$stable_id2)

 Function: return 1 if the two gene are returned with the same strand
           by using Bio::EnsEMBL::Registry->fetch_by_stable_id()
           or -1 if they don't...

 Returns : 1 or -1

 Args    : -1 gene stable id
           -2 gene stable id

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!!!

=cut

sub relative_strand {

  my $class = shift;
  my $id1 = shift;
  my $id2 = shift;
  my $specie1 = Bio::MCE::Utils->get_genome_from_id($id1);
  my $specie2 = Bio::MCE::Utils->get_genome_from_id($id2);

  my $gene_adaptor1 = Bio::EnsEMBL::Registry->get_adaptor($specie1,"core","gene");
  my $gene_adaptor2 = Bio::EnsEMBL::Registry->get_adaptor($specie2,"core","gene");

  my $gene1 = $gene_adaptor1->fetch_by_stable_id($id1);
  my $gene2 = $gene_adaptor2->fetch_by_stable_id($id2);

  my $strand1 = $gene1->strand;
  my $strand2 = $gene2->strand;

  my $rel_strand = $strand1 * $strand2;

  return($rel_strand);
}


=head2 get_xref_from_id_constraint

 Title   : get_xref_from_id_constraint

 Usage   : Bio::PMCE::Utils->get_xref_from_id_constraint($ref_id,'aniseed')

 Function: return the xref as an array of Bio::EnsEMBL::DBEntry object for all the
           entry coming from the database matching the constraint into the dbname of
           Bio::EnsEMBL::DBEntry object. Usually you would need to call from the
           Bio::EnsEMBL::DBEntry object:
           $entry->dbname
           $entry->display_id

 Returns : an array of objects

 Args    : -1 gene stable id
           -2 a constraint to match the dbname as Marker (MGI) or GO or other...

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!!!

=cut

sub get_xref_from_id_constraint {

	my $self = shift;

  my $ID = shift;

  my $CONSTRAINT = shift;

  my @res;

  my $ORGANISM = Bio::MCE::Utils->get_genome_from_id($ID);

  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($ORGANISM,'core','gene');

  my $gene = $gene_adaptor->fetch_by_stable_id($ID);

  my $dbentries = $gene->get_all_DBLinks;

  foreach my $xref(@$dbentries) {
    next unless $xref->dbname =~ "$CONSTRAINT";
    push(@res,$xref);
  }
  return(\@res) if scalar(@res);
	return undef;
}

=head2 get_aln_string_from_simple_align

 Title   : get_aln_string_from_simple_align

 Usage   : Bio::PMCE::Utils->get_aln_string_from_simple_align($simple_align)

 Function: converter to obtain a string with the alignment from a SimpleAlign
           object in ordeer to put this into a database

 Returns : a string

 Args    : a SimpleAalign object

 Note    : 

=cut

sub get_aln_string_from_simple_align {

	my $self = shift;
  my $simple_align = shift;
  my $alnio = '';

  open(MEMORY,'>', \$alnio);
  
  my $fh = Bio::AlignIO->new(-format => 'clustalw',
                             -fh => \*MEMORY);
  
  $fh->write_aln($simple_align);
  
  close(MEMORY);

  return($alnio);
}

=head2 get_coord_sys_from_ref  

 Title   : get_coord_sys_from_ref
 
 Usage   : my $csn = Bio::MCE::Utils->get_coord_sys_from_ref($ref,$DBH);
 
 Function: 
 
 Returns : a string
           
 Args    : -1 a ref_id
           -2 DBH to valis

=cut

sub get_coord_sys_from_ref {

  my $class = shift;
  my $ref = shift;
  my $DBH = shift;

  my $sth = $DBH->prepare('SELECT coord_sys_name FROM seq WHERE ref_id = '.$DBH->quote($ref));
  $sth->execute;
  my $row = $sth->fetchrow_hashref;
  my $csn = $row->{'coord_sys_name'};

  return $csn;
}

1;
