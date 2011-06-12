#!/usr/bin/perl -w

=head1 Bio::TrapCluster::Region

=head2 Authors

=head3 Created by

             Vincenza Maselli
             v.maselli@ucl.ac.uk


=head2 Description
        
             
=head2 Usage

	    my $obj = Bio::TrapCluster::Region->new;
            
=cut

package Bio::TrapCluster::Region;

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
use Bio::EnsEMBL::Intron;

@ISA = qw(Bio::Root::Root );

#setting global variables

require "$ENV{'TrapCluster'}/trapcluster_conf.pl";
my %conf =  %::conf;
my $debug = $conf{'global'}{'debug'};


sub new{
  	my $caller = shift;

 	my $class = ref($caller) || $caller;
  	my $self = $class->SUPER::new(@_);
	
 my (
   		$load,
   		$registry
    )
    = $self->_rearrange( [
      	'LOAD',
      	'REGISTRY'
    ],
    @_
    );
	
	
	$load && $self->load($load);
	$registry && $self->registry($registry);
  	return $self;
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

sub registry{

  my ($self, $value) = @_;
  $self->{'registry'} = $value if defined $value;
  
  return $self->{'registry'};
}



sub slice{
	my ($self, $slice) = @_;
	my $region_hash;
	$region_hash->{'name'} = $slice->seq_region_name;
	$region_hash->{'seq_id'} = $slice->seq_region_name;
	$region_hash->{'strand'} = $slice->strand;
	$region_hash->{'start'} = $slice->start;
	$region_hash->{'end'} = $slice->end;
	$region_hash->{'parent_id'} = 0;
	$region_hash->{'rank'} = 0;
	$region_hash->{'description'} = "genomic";
	$region_hash->{'biotype'} = "genomic";
	my $region_id = $self->load->load_region($region_hash);
	my $failure = "insert failure for ".$region_hash->{'name'};
	unless($region_id){my $error_id =$self->load->load_error($failure,"region");warning("$failure")}
			
	return $region_id
}

sub gene{
	my ($self, $gene) = @_;
	
	my $hash = $self->_prepare_region_hash($gene);
	$hash->{'refseq'} = $self->_get_refseq_by_gene($gene);
	$hash->{'parent_id'} = 0;
	$hash->{'rank'} = 0;
	$hash->{'description'} = "ensembl gene";
	$hash->{'biotype'} = $gene->biotype;

	my $external_name = $gene->external_name;
	unless (defined $external_name){
		$external_name = $hash->{'refseq'} if $hash->{'refseq'};
		$external_name = $gene->stable_id unless $hash->{'refseq'};
	}
	$hash->{'external_name'} = $external_name;
	my $gene_id = $self->load->load_region($hash);
	my $failure = "insert failure for ".$hash->{'name'};
	unless($gene_id){my $error_id =$self->load->load_error($failure,"region");warning("$failure")}
	
	return $gene_id;
}

sub up_downstream{
	my ($self, $slice,$gene,$f) = @_;
	
	my $hash = $self->_prepare_region_hash($slice,$gene);
	$hash->{'parent_id'} = $self->load->fetch->get_region_by_seq_id($gene->stable_id)->{'region_id'};
	$hash->{'rank'} = 0;
	$hash->{'description'} = $f->display_name;
	$hash->{'start'} = $f->start;
	$hash->{'end'} = $f->end;
	
	my $utr_id = $self->load->load_region($hash);
	my $failure = "insert failure for ".$hash->{'name'};
	unless($utr_id){my $error_id =$self->load->load_error($failure,"region");warning("$failure")}
	
	return $utr_id;
}

sub exon_array {
	my ($self, $exarray, $gene_id) = @_;
	my $count = 1;
	foreach my $exon (@{$exarray}){
		my $slice = $exon->slice;
		$self->exon($slice, $exon,$count,$gene_id);
		$count ++;
	}
	return 1;
}

sub exon {
	my ($self,$exon,$exoncount,$gene_id) = @_;
								
	my $hash = $self->_prepare_region_hash($exon);
	$hash->{'parent_id'} = $gene_id;
	$hash->{'rank'} = $exoncount;
	#defining coding start and coding end
	my $transcript_adaptor = $self->registry->get_adaptor( 'Mouse', 'Core', 'Transcript' ) || die $!;
	my ($transcript) = @{$transcript_adaptor->fetch_all_by_exon_stable_id($exon->stable_id)};
	my $coding_start = $transcript->coding_region_start;
	my $coding_end = $transcript->coding_region_end;
	my $exon_type;
	my $estart = $hash->{'start'};
	my $eend = $hash->{'end'};
	if ($coding_start && $coding_end){
        if ($eend < $coding_start){$exon_type = "5' UTR"}
        elsif ($estart < $coding_start && $coding_start < $eend){$exon_type = "5' UTR-EXON"}
        elsif($estart > $coding_end){$exon_type = "3' UTR"}
        elsif($estart < $coding_end && $coding_end < $eend){$exon_type = "3' EXON-UTR"}
        else{$exon_type = "CODING-EXON"}
    }
	else{$exon_type = "PSEUDOGENE-EXON"}
	
	$hash->{'description'} = $exon_type;
	
	my $exon_id = $self->load->load_region($hash);
	my $failure = "insert failure for exon ".$hash->{'name'};
	unless($exon_id){my $error_id =$self->load->load_error($failure,"region");die("$failure\n")}
	$hash->{'region_id'} = $exon_id;
	
	return $exon_id;
}

sub est_exon {
	my ($self,$exon,$exoncount,$gene_id) = @_;
								
	my $hash = $self->_prepare_region_hash($exon);
	$hash->{'parent_id'} = $gene_id;
	$hash->{'rank'} = $exoncount;
	#defining coding start and coding end
	print "got here est_exon 206\n";
	my $transcript_adaptor = $self->registry->get_adaptor( 'Mouse', 'Otherfeatures', 'Transcript' ) || die $!;
	
	my $transcript;
	if (defined $exon->stable_id){
		($transcript) = @{$transcript_adaptor->fetch_all_by_exon_stable_id($exon->stable_id)};
	}
	else{
		($transcript) = @{$transcript_adaptor->fetch_all};
	}	
	
	my $coding_start = $transcript->coding_region_start;
	my $coding_end = $transcript->coding_region_end;
	my $exon_type;
	my $estart = $hash->{'start'};
	my $eend = $hash->{'end'};
	if ($coding_start && $coding_end){
        if ($eend < $coding_start){$exon_type = "5' UTR"}
        elsif ($estart < $coding_start && $coding_start < $eend){$exon_type = "5' UTR-EXON"}
        elsif($estart > $coding_end){$exon_type = "3' UTR"}
        elsif($estart < $coding_end && $coding_end < $eend){$exon_type = "3' EXON-UTR"}
        else{$exon_type = "CODING-EXON"}
    }
	else{$exon_type = "PSEUDOGENE-EXON"}
	
	$hash->{'description'} = $exon_type;

	my $exon_id = $self->load->load_region($hash);
	my $failure = "insert failure for exon ".$hash->{'name'};
	unless($exon_id){my $error_id =$self->load->load_error($failure,"region");die("$failure\n")}
	$hash->{'region_id'} = $exon_id;
	die unless $exon_id;
	return $exon_id;
}


sub intron{
	my ($self,$intron,$start,$end,$exoncount,$gene_id) = @_;

	my $hash = $self->_prepare_region_hash($intron);
	
	my $prev_st_id = $intron->prev_Exon->stable_id;
	unless (defined $prev_st_id){$prev_st_id = "Exon.".$gene_id.".".$exoncount;}
	my $next_st_id = $intron->next_Exon->stable_id;
	unless (defined $next_st_id){$next_st_id = "Exon.".$gene_id.".".($exoncount+1);}
	$hash->{'seq_id'} = $prev_st_id."_".$next_st_id;
	$hash->{'start'} = $start;
	$hash->{'end'} = $end;
	$hash->{'parent_id'} = $gene_id;
	$hash->{'rank'} = $exoncount;
	$hash->{'description'} = "minimal intron";
	my $intron_id = $self->load->load_region($hash);	
	
	my $failure = "insert failure for intron ".$hash->{'name'};
	unless($intron_id){my $error_id =$self->load->load_error($failure,"region");die("$failure\n")}
	$hash->{'region_id'} = $intron_id;
	return $intron_id;
}


sub _prepare_region_hash{
	my ($self,  $obj) = @_;
	
	my $hash;
	
	$hash->{'name'} = $obj->seq_region_name;
	$hash->{'seq_id'} = $obj->stable_id unless $obj->isa("Bio::EnsEMBL::Intron");
	$hash->{'strand'} = $obj->strand;
	$hash->{'start'} = $obj->seq_region_start;
	$hash->{'end'} = $obj->seq_region_end;
		
	return $hash;

}




sub rearrange_exons{
	my ($self, $gene, $gene_id) = @_;
	$debug && print STDERR ref $self;
	$debug && print STDERR "->rearrange_exons\n";
	my $gene_stable_id = $gene->stable_id;
	my $sql = qq{select * from region where parent_id = $gene_id };
	#check for already stored exons
	my $res = $self->load->fetch->select_many_from_table($sql);
	my @regions = sort {$b->{'region_end'} <=> $a->{'region_end'} && $a->{'region_start'} <=> $b->{'region_start'} && $a->{'rank'} <=> $b->{'rank'}} @{$res};
	
	if (scalar @regions){return \@regions}
 	
	$debug && print STDERR "NEW GENE ";
	my @exons = @{$gene->get_all_Exons};
	$debug && print STDERR "with ", scalar @exons," exons\n";
	my ($exon1,$exon2);
	my ($intron_start, $intron_end);
	my $min_start;
	my $max_end;
	my $prev_end;
	my $exoncount = 1;
	my $chosen;
	my $test = 1;
	my @noncoding;
	my @coding;
	
	foreach my $exon (@exons){		
		my $transcript_adaptor = $self->registry->get_adaptor( 'Mouse', 'Core', 'Transcript' ) || die $!;
		my @transcripts = @{$transcript_adaptor->fetch_all_by_exon_stable_id($exon->stable_id)};
		my $transcript = $transcripts[0];
		if ($gene->biotype =~ /coding/ && $transcript->biotype !~ /coding/ && scalar @transcripts > 1){
			push (@noncoding, $exon);
		}
		else{push (@coding,$exon)}
	
	}
	#$self->exon_array(\@noncoding, $gene_id);
	
	my @sorted_exons = sort {$a->start <=> $b->start} @coding;
	#my @sorted_exons = @coding;
	$debug && print STDERR "saved ".scalar @sorted_exons."\n";
	foreach my $exon (@sorted_exons){
		my $estart = $exon->seq_region_start;
		my $eend = $exon->seq_region_end;
		
		#first exon

		if ($test == 1){
			$min_start = $estart;
			$max_end = $eend;
			$chosen->{$exoncount}{'exon'} = $exon;
			my $id = $self->exon($exon,$exoncount,$gene_id);
			$test ++;
		}
		else{
			if ($estart < $max_end){
				#if is still the same slot
				if ($eend > $max_end){
					#the next exon is more oustanding
					$max_end = $eend; 
					$chosen->{$exoncount}{'exon'} = $exon;
				}
				my $id = $self->exon($exon,$exoncount,$gene_id);
			}
			else{
				#next slot
				$min_start = $estart;
				$exoncount ++;
				$max_end = $eend;
				$chosen->{$exoncount}{'exon'} = $exon;
				my $id = $self->exon($exon,$exoncount,$gene_id);
			}	
		}
		$chosen->{$exoncount}{'start'} = $min_start;
		$chosen->{$exoncount}{'end'} = $max_end;
	}
	
	foreach my $rank (keys %{$chosen}){	
		my $exon1 = $chosen->{$rank}{'exon'};
		my $exon2 = $chosen->{$rank + 1}{'exon'};
		if ($exon2){
			$intron_start = $chosen->{$rank}{'end'} + 1;
			$intron_end = $chosen->{$rank + 1}{'start'} - 1;
			if ($intron_start < $intron_end){
				my $intron = Bio::EnsEMBL::Intron->new($exon1,$exon2);
				my $id = $self->intron($intron,$intron_start,$intron_end,$rank,$gene_id);
			}
			else{
				$debug && print STDERR "SKIPPED because $intron_start > $intron_end\n";
			}
 		}
	}
	$res = $self->load->fetch->select_many_from_table($sql);
	@regions = sort {$a->{'rank'} <=> $b->{'rank'} && $a->{'region_start'} <=> $b->{'region_start'}} @{$res};
	
	if ($gene->biotype =~ /coding/){
		my $first_start = $regions[0]->{'region_start'}; 
		my $last_end = $regions[$#regions]->{'region_end'};
		my $update = qq{UPDATE region SET region_start = $first_start, region_end = $last_end WHERE seq_id = \"$gene_stable_id\"};
		$self->load->fetch->update($update);
	}
	return \@regions;
}


sub rearrange_est_exons{
	my ($self, $gene, $gene_id) = @_;
	$debug && print STDERR ref $self;
	$debug && print STDERR "->rearrange_est_exons\n";
	
	my $gene_stable_id = $gene->stable_id;
	$debug && print STDERR "got here 348\n";
	my $sql = qq{select * from region where parent_id = $gene_id };
	$debug && print STDERR "$sql\n";
	#check for already stored exons
	my $res = $self->load->fetch->select_many_from_table($sql);
	my @regions = sort {$b->{'region_end'} <=> $a->{'region_end'} && $a->{'region_start'} <=> $b->{'region_start'} && $a->{'rank'} <=> $b->{'rank'}} @{$res};
	
	if (scalar @regions){return \@regions}
 	
	$debug && print STDERR "NEW GENE ";
	my @exons = @{$gene->get_all_Exons};
	$debug && print STDERR "with ", scalar @exons," exons\n";
	my ($exon1,$exon2);
	my ($intron_start, $intron_end);
	my $min_start;
	my $max_end;
	my $prev_end;
	my $exoncount = 1;
	my $chosen;
	my $test = 1;
	
	my @sorted_exons = sort {$a->start <=> $b->start} @exons;
	$debug && print STDERR "saved ".scalar @sorted_exons."\n";
	foreach my $exon (@sorted_exons){
		my $estart = $exon->seq_region_start;
		my $eend = $exon->seq_region_end;
		$debug && print "S $estart E $eend\n"; 
		#first exon

		if ($test == 1){
			$min_start = $estart;
			$max_end = $eend;
			$chosen->{$exoncount}{'exon'} = $exon;
			$debug && print "got here 381\n";
			my $id = $self->est_exon($exon,$exoncount,$gene_id);
			$test ++;
		
		}
		else{
			if ($estart < $max_end){
				#if is still the same slot
				if ($eend > $max_end){
					#the next exon is more oustanding
					$max_end = $eend; 
					$chosen->{$exoncount}{'exon'} = $exon;
				}
				my $id = $self->est_exon($exon,$exoncount,$gene_id);
			}
			else{
				#next slot
				$min_start = $estart;
				$exoncount ++;
				$max_end = $eend;
				$chosen->{$exoncount}{'exon'} = $exon;
				my $id = $self->est_exon($exon,$exoncount,$gene_id);
			}	
		}
		$chosen->{$exoncount}{'start'} = $min_start;
		$chosen->{$exoncount}{'end'} = $max_end;
	}
	
	foreach my $rank (keys %{$chosen}){	
		$debug && print STDERR "EXON RANK $rank\n";
		my $exon1 = $chosen->{$rank}{'exon'};
		my $exon2 = $chosen->{$rank + 1}{'exon'};
		if ($exon2){
			$intron_start = $chosen->{$rank}{'end'} + 1;
			$intron_end = $chosen->{$rank + 1}{'start'} - 1;
			if ($intron_start < $intron_end){
				my $intron = Bio::EnsEMBL::Intron->new($exon1,$exon2);
				my $id = $self->intron($intron,$intron_start,$intron_end,$rank,$gene_id);
			}
			else{
				$debug && print STDERR "SKIPPED because $intron_start > $intron_end\n";
			}
 		}
	}
	$res = $self->load->fetch->select_many_from_table($sql);
	@regions = sort {$a->{'rank'} <=> $b->{'rank'} && $a->{'region_start'} <=> $b->{'region_start'}} @{$res};
	
	return \@regions;
}

sub rearrange_predict_exons{
	my ($self, $gene, $gene_id) = @_;
	$debug && print STDERR ref $self;
	$debug && print STDERR "->rearrange_est_exons\n";
	
	my $gene_stable_id = $gene->stable_id;
	$debug && print STDERR "got here 348\n";
	my $sql = qq{select * from region where parent_id = $gene_id };
	$debug && print STDERR "$sql\n";
	#check for already stored exons
	my $res = $self->load->fetch->select_many_from_table($sql);
	my @regions = sort {$b->{'region_end'} <=> $a->{'region_end'} && $a->{'region_start'} <=> $b->{'region_start'} && $a->{'rank'} <=> $b->{'rank'}} @{$res};
	
	if (scalar @regions){return \@regions}
 	
	$debug && print STDERR "NEW GENE ";
	my @exons = @{$gene->get_all_Exons};
	$debug && print STDERR "with ", scalar @exons," exons\n";
	my ($exon1,$exon2);
	my ($intron_start, $intron_end);
	my $min_start;
	my $max_end;
	my $prev_end;
	my $exoncount = 1;
	my $chosen;
	my $test = 1;
	
	my @sorted_exons = sort {$a->start <=> $b->start} @exons;
	$debug && print STDERR "saved ".scalar @sorted_exons."\n";
	foreach my $exon (@sorted_exons){
		my $estart = $exon->seq_region_start;
		my $eend = $exon->seq_region_end;
		$debug && print "S $estart E $eend\n"; 
		#first exon

		if ($test == 1){
			$min_start = $estart;
			$max_end = $eend;
			$chosen->{$exoncount}{'exon'} = $exon;
			$debug && print "got here 381\n";
			my $id = $self->est_exon($exon,$exoncount,$gene_id);
			$test ++;
		
		}
		else{
			if ($estart < $max_end){
				#if is still the same slot
				if ($eend > $max_end){
					#the next exon is more oustanding
					$max_end = $eend; 
					$chosen->{$exoncount}{'exon'} = $exon;
				}
				my $id = $self->est_exon($exon,$exoncount,$gene_id);
			}
			else{
				#next slot
				$min_start = $estart;
				$exoncount ++;
				$max_end = $eend;
				$chosen->{$exoncount}{'exon'} = $exon;
				my $id = $self->est_exon($exon,$exoncount,$gene_id);
			}	
		}
		$chosen->{$exoncount}{'start'} = $min_start;
		$chosen->{$exoncount}{'end'} = $max_end;
	}
	
	foreach my $rank (keys %{$chosen}){	
		$debug && print STDERR "EXON RANK $rank\n";
		my $exon1 = $chosen->{$rank}{'exon'};
		my $exon2 = $chosen->{$rank + 1}{'exon'};
		if ($exon2){
			$intron_start = $chosen->{$rank}{'end'} + 1;
			$intron_end = $chosen->{$rank + 1}{'start'} - 1;
			if ($intron_start < $intron_end){
				my $intron = Bio::EnsEMBL::Intron->new($exon1,$exon2);
				my $id = $self->intron($intron,$intron_start,$intron_end,$rank,$gene_id);
			}
			else{
				$debug && print STDERR "SKIPPED because $intron_start > $intron_end\n";
			}
 		}
	}
	$res = $self->load->fetch->select_many_from_table($sql);
	@regions = sort {$a->{'rank'} <=> $b->{'rank'} && $a->{'region_start'} <=> $b->{'region_start'}} @{$res};
	
	return \@regions;
}

sub _get_refseq_by_gene{
    my ($self, $gene) = @_;	
	my @dblinks = @{$gene->get_all_DBLinks()};
	my $refseq_id;
	foreach my $link (@dblinks) {
		if ($link->database eq 'RefSeq_dna') {
			$refseq_id = $link->primary_id;
			if ($refseq_id =~ /NM/) {return $refseq_id;}
		} 
		elsif ($link->database eq 'RefSeq_peptide') {$refseq_id = $link->primary_id;}
	}
	unless (defined $refseq_id){$refseq_id = "NA"}
	return $refseq_id;
        
}

sub get_label{
	my ($self,$gene) = @_;
	die unless $gene;
	my $refseq_id = $self->_get_refseq_by_gene($gene);
	my $logic_name;
	if ($gene->stable_id =~ "ENSMUSG") {
		$logic_name = 'ENSG';
		if ($refseq_id ne "NA"){$logic_name = 'REFSEQ'}	
	}
	elsif ($gene->stable_id =~ "ENSMUSESTG") {
		$logic_name = 'ENSESTG';
	}
	return $logic_name;
}


1;