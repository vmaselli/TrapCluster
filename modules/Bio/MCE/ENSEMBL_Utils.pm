#!/usr/bin/perl
package Bio::MCE::Ensembl::Utils;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::MCE::Utils;

=head2 get_all_domains

 Title   : get_all_domains

 Usage   : my $domains = Bio::MCE::Ensembl::Utils->get_all_domains($id,$specie);

 Function: take a gene stable id, the specie and returns all the domains of all the translations

 Returns : a listref whith strings

 Args    : -1 a gene stable id
           -2 the name of the specie

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub get_all_domains {
	my $class = shift;
	my $id = shift;
	my $specie = shift;
	my $DOM;
	my @DOM;

  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","gene") or die "ERROR ADAPTOR $!\n";
  my $gene = $gene_adaptor->fetch_by_stable_id($id) or die "ERROR FETCH $!\n";

  my @transcript = @{ $gene->get_all_Transcripts };
  foreach my $ele(@transcript) {
    my $translation = $ele->translation;
		next unless $translation;
    my @dom = @{$translation->get_all_DomainFeatures()};
    foreach my $dom(@dom) {
			next unless $dom->idesc =~ /.+/;
			$DOM->{$dom->idesc} = 1;
    }
  }
	foreach my $key(keys %$DOM) {
		push(@DOM,$key);
	}
	return(\@DOM);
}

=head2 get_common_domains

 Title   : get_common_domains

 Usage   : my $domains = Bio::MCE::Ensembl::Utils->get_common_domains(\@stable_ids);

 Function: returns the names of the common domains among the genes you are asking for

 Returns : an hashref with the names of common domain as keys

 Args    : -1 an array of gene stable id

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub get_common_domains {

	my $class = shift;
	my $ref = shift;

	my @id;
	foreach my $id(@$ref) {
		push(@id,$id);
	}

  my $href;
  my $res;

  foreach my $id(@id) {
    my $specie = Bio::MCE::Utils->get_genome_from_id($id);
    my $domains = Bio::MCE::Ensembl::Utils->get_all_domains($id,$specie);
    foreach my $domain(@$domains) {
      $href->{$id}->{$domain} ++;
    }
  }

  my $query = shift(@id);
  foreach my $key(keys %{$href->{$query}}) {
    my $c;
    foreach my $id(@id) {
      my $used = 'NO';
      if(exists $href->{$id}->{$key} && $used eq 'NO') {
        $c->{$key} ++;
        $used = 'YES';
      }
    }
    if($c->{$key}) {
      $res->{$key} ++ if $c->{$key} == scalar(@id);
    }
  }
  return $res;
}

=head2 get_common_domains_from_slices

 Title   : get_common_domains_from_slices

 Usage   : my $domains = Bio::MCE::Ensembl::Utils->get_common_domains_from_slices($slice1,$slice2);

 Function: returns the names of the common domains among the genes on the slices you are asking for

 Returns : an hashref with the names of common domain as keys

 Args    : -1 the first slice
           -2 the second slice

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub get_common_domains_from_slices {

  my $class = shift;
  my $slice1 = shift;
  my $slice2 = shift;

  my $gene_dom;

  # this can be simply updated to work with more than two slices simply using an arrayref
  # instead of '1' and '2' and modifing a little the loop... BUT NOW I'M IN HURRY!!!!!!!!
  $gene_dom->{'1'}->{'genes'} = $slice1->get_all_Genes_by_type('protein_coding');
  $gene_dom->{'2'}->{'genes'} = $slice2->get_all_Genes_by_type('protein_coding');

  my $href;
  my $res;

  foreach my $number(keys %$gene_dom) {
    my $genes = $gene_dom->{$number}->{'genes'};
    foreach my $gene(@$genes) {
      my @transcript = @{ $gene->get_all_Transcripts };
      foreach my $ele(@transcript) {
        my $translation = $ele->translation;
        my @dom = @{$translation->get_all_DomainFeatures()};
        foreach my $dom(@dom) {
          next unless $dom->idesc =~ /.+/;
          $gene_dom->{$number}->{'domains'}->{$dom->idesc} = 1;
        }
      }
    }
  }

  my $query = $gene_dom->{'1'}->{'domains'};
  my $target = $gene_dom->{'2'}->{'domains'};

  foreach my $key(keys %$query) {
    my $c;
    my $used = 'NO';
    if(exists $target->{$key} && $used eq 'NO') {
      $c->{$key} ++;
      $used = 'YES';
    }
    if($c->{$key}) {
      $res->{$key} ++;
    }
  }
  return $res;
}

=head2 get_all_rank1_coord_length_hashref
 
 Title   : get_all_rank1_coord_length_hashref
    
 Usage   : my $hashref = Bio::MCE::Ensembl::Utils->get_all_rank1_coord_length_hashref($specie);

 Function: give you an hashref of hashref. The keys are the name of the rank 1 coordinate system
           and in each element you have other two keys: length and coord_system_name.
           Here is the print Dumper for the human genome as example:

           $VAR1 = {
                      '11' => {
                      'length' => 134452384,
                      'coord_sys_name' => 'chromosome'
                              },
                      '21' => {
                      'length' => 46944323,
                      'coord_sys_name' => 'chromosome'
                             },
                      '7' => {
                      'length' => 158821424,
                      'coord_sys_name' => 'chromosome'
                             },
                      'Y' => {
                      'length' => 54733917,
                      'coord_sys_name' => 'chromosome'
                    }................

  
 Returns : hashref of hashref

 Args    : -1 a string: the name of the specie
    
 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub get_all_rank1_coord_length_hashref {

	my $class = shift;
	my $specie = shift;

	my $csa = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","coordsystem") or die "ERROR COORDSYS ADAPTOR $!\n";
	my $cs = $csa->fetch_by_rank(1);

	my $sa = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice") or die "ERROR SLICE ADAPTOR $!\n";
	my $s = $sa->fetch_all($cs->name);

	my $res;

	foreach my $chr(@$s) {
	  next if $chr->seq_region_name eq 'MT';
	  $res->{$chr->seq_region_name}->{'length'} = $chr->length;
	  $res->{$chr->seq_region_name}->{'coord_sys_name'} = $chr->coord_system_name;
	}

	return $res;
}

=head2 get_random_chr_length
           
 Title   : get_random_chr_length
                      
 Usage   : my $domains = Bio::MCE::Ensembl::Utils->get_random_chr_length($specie);
                                
 Function: return an hashref with 3 keys: chr, coord_sys_name and length
                      
 Returns : an hashref
  
 Args    : -1 the name of the specie

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!
   
=cut                  

sub get_random_chr_length {

	my $class = shift;
	my $specie = shift;

	my $hashref = Bio::MCE::Ensembl::Utils->get_all_rank1_coord_length_hashref($specie);

	my @regions = keys %{$hashref};

	# it's ok because from perldoc:
	# int(rand(10)) returns a random integer between 0 and 9, inclusive.
	my $inx = int(rand(scalar(@regions)));

	my $chr = $regions[$inx];
	
	my $rand_sel = $hashref->{$chr};

	$rand_sel->{'chr'} = $chr;

	return $rand_sel;
}


=head2 get_random_noncoding_norepeat_noN

 Title   : get_random_noncoding_norepeat_noN

 Usage   : my($slice,$chr,$start,$end) = Bio::MCE::Ensembl::Utils->get_random_noncoding_norepeat_noN($specie,$length);
 
 Function: look for random region of required length from a random chromosome of the required specie
           and check if there are repeats, proteins, ESTs, prediction, 'N's. If not it returns the 
           selected random region.

 Returns : -1 the slice
           -2 the chr
           -3 the genomic start
           -4 the genomic end

 Args    : -1 the name of the specie
           -2 the length of the desidered random region

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!
           THE FUNCTION LOOK BOTH IN ENSEMBL CORE AND IN ENSEMBL OTHERFEATURES SO REMEMBER TO LOAD
           BOTH THE DATABASES WHEN YOU LOAD THE REGISTRY

=cut

sub get_random_noncoding_norepeat_noN {

	my $class = shift;
	my $specie = shift;
	my $length = shift;
	my $res = 'NO';

	while($res eq 'NO') {

		my $hashref = Bio::MCE::Ensembl::Utils->get_random_chr_length($specie);

		my $srn = $hashref->{'chr'};
		my $csn = $hashref->{'coord_sys_name'};
		my $srl = $hashref->{'length'};
	
		my $limit = $srl-$length;

		my $rand_start = int(rand($limit)+1);
		my $rand_end = $rand_start+$length;

		my $sa = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice") or die "ERROR SLICE ADAPTOR $!\n";
		my $slice = $sa->fetch_by_region($csn,$srn,$rand_start,$rand_end);

		my $seq = $slice->seq;

		my $pn = int(Bio::MCE::Utils->seq_perc($seq,'N'));
		if($pn > 0) {
			next;
		}

		my $ra = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","repeatfeature") or die "ERROR REPEAT ADAPTOR $!\n";
		my $rf = $ra->fetch_all_by_Slice($slice);
		if(scalar(@$rf) > 0) {
			next;
		}

		my($coding,$l) = Bio::MCE::Pipeline::Coding->test_coding($specie,$csn,$srn,$rand_start,$rand_end);
		if($coding > 0) {
			next;
		} 

		else {
			$res = 'YES'; # not really useful...
			return($slice,$srn,$rand_start,$rand_end);
		}
	}
}

=head2 get_window

 Title   : get_window

 Usage   : my($start,$end) = Bio::MCE::Ensembl::Utils->get_window($specie,$chr,$start,$end,$wl);

 Function: calculate and returns the start and end of a window of length wl around the central point of 
           coordinate that you ask for

 Returns : two scalars: start and end

 Args    : -1 the name of the specie
           -2 the name of the chr (seq_region_name)
           -3 the start (seq_region_start)
           -4 the end (seq_region_end)

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub get_window {

  my $class = shift;
  my $specie = shift;
  my $chr = shift;
  my $pre_start = shift;
  my $pre_end = shift;
  my $window_length = shift;

  my $chr_slice = Bio::MCE::Utils->get_chr($chr,$specie);
  my $chr_length = $chr_slice->length;

  my $sce_length = $pre_end-$pre_start+1;

  my $flanking = int($window_length / 2);
  $flanking = 1 if $flanking < 1;

	my $middle_point = int($pre_start + ($sce_length/2));

	my $start = $middle_point - $flanking + 1;
  $start = 1 if $start < 1;

	my $end = $middle_point + $flanking - 1;
  $end = $chr_length if $end > $chr_length;

  return($start,$end);
}  

=head2 get_smaller_distance

 Title   : get_smaller_distance

 Usage   : It take two pair in this order: start1, end1 and start2, end2
           of coord and return the minimum distance between the two segments
           even if they overlap in fact it sorts the coord before of all
 
 Function: 

 Returns : THE MINIMUM DISTANCE BETWEEN THE TWO SETS OF COORDINATES

 Args    : -1 start1
           -2 end1
           -3 start2
           -4 end2

 NOTE    : THE START HAVE TO BE ALWAYS SMALLER THAN END

=cut

sub get_smaller_distance {

	my $class = shift;
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


=head2 get_random_slice

 Title   : get_random_slice

 Usage   : my($slice,$chr,$start,$end) = Bio::MCE::Ensembl::Utils->get_random_slice($specie,$length);
 
 Function: look for random region of required length from a random chromosome of the required specie
           without checking if there are repeats, proteins, ESTs, prediction, 'N's.

 Returns : -1 the slice
           -2 the chr
           -3 the genomic start
           -4 the genomic end

 Args    : -1 the name of the specie
           -2 the length of the desidered random region

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub get_random_slice {

	my $class = shift;
  my $specie = shift;
  my $length = shift;

  my $hashref = Bio::MCE::Ensembl::Utils->get_random_chr_length($specie);

  my $srn = $hashref->{'chr'};
  my $csn = $hashref->{'coord_sys_name'};
  my $srl = $hashref->{'length'};

  my $limit = $srl-$length;

  my $rand_start = int(rand($limit)+1);
  my $rand_end = $rand_start+$length;

  my $sa = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice") or die "ERROR SLICE ADAPTOR $!\n";
  my $slice = $sa->fetch_by_region($csn,$srn,$rand_start,$rand_end);

  return($slice,$srn,$rand_start,$rand_end);
}
 

=head2 get_seq_region_length

 Title   : get_seq_region_length

 Usage   : my $length = Bio::MCE::Ensembl::Utils->get_seq_region_length($specie,$coord_sys_name,$seq_region_name);

 Function: returns the length of the requested region

 Returns : a scalar

 Args    : -1 the name of the specie
           -2 the coordinate system name
           -3 the seq region name

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub get_seq_region_length {

  my $class = shift;
  my $specie = shift;
  my $coord_system_name = shift;
  my $srn = shift;

  my $sa = Bio::EnsEMBL::Registry->get_adaptor($specie,"core","slice");
  my $s = $sa->fetch_by_region($coord_system_name,$srn);
     
  return $s->length;
}


=head2 get_other_db_ids

 Title   : get_other_db_ids

 Usage   : my $ids = Bio::MCE::Ensembl::Utils->get_other_db_ids('Danio rerio','ENSDARG00000037777','ZFIN_ID');

 Function: take a gene stable id, the specie, the db and returns the ids of the gene in the other db

 Returns : a hashref with keys of ids

 Args    : -1 the name of the specie
           -2 a gene stable id
           -3 the name of the db (as in dbname() of the DBEntry object)

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!

=cut

sub get_other_db_id {
	my $class = shift;
  my $organism = shift;
  my $gene_id = shift;
  my $db = shift;
  my $res;
  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($organism,"core","gene");
  my $gene = $gene_adaptor->fetch_by_stable_id($gene_id);
  my @all = @{$gene->get_all_DBLinks};
  foreach my $entry(@all) {
    #print $entry->dbname."\n";
    if($entry->dbname eq $db) {
      $res->{$entry->primary_id} ++;
    }
  }
  return($res);
}

=head2 get_family

 Title   : get_family

 Usage   : my $ids = Bio::MCE::Ensembl::Utils->get_family('ENSDARG00000037777');

 Function: take a gene stable id and returns the id and the description of the family

 Returns : a hashref of id => description

 Args    : -1 the gene stable_id

 Note    : THE REGISTRY MUST HAVE BEEN LOADED IN YOUR SCRIPT!
           THE REGISTRY MUST CONTAINS THE COMPARA DB!!!!

=cut

sub get_family {
  my $class = shift;
  my $id = shift;
  my $res = {};
  my $member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Member');
  my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',$id);
	return undef unless $member;
  my $family_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Family');
  my $families = $family_adaptor->fetch_all_by_Member($member);
	return undef unless $families;
  foreach my $family (@{$families}) {
    $res->{$family->stable_id} = $family->description;
  }
  return $res;
}

1;
