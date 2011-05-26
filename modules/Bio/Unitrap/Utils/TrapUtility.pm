#!/usr/bin/perl -w

=head1 Bio::Unitrap::Utils::TrapUtility

=head2 Authors

=head3 Created by

             Vincenza Maselli
             v.maselli@ucl.ac.uk


=head2 Description
        
             
=head2 Usage

	    my $obj = Bio::Unitrap::Utils::TrapUtility->new;
            
=cut

package Bio::Unitrap::Utils::TrapUtility;

require "$ENV{'Unitrap'}/unitrap_conf.pl";
my %conf =  %::conf;
my $debug = $conf{'global'}{'debug'};

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
use Bio::SeqIO;

use Bio::Unitrap::Region;

@ISA = qw(Bio::Root::Root Bio::SeqIO);

#setting global variables

require "$ENV{'Unitrap'}/unitrap_conf.pl";

sub new{
  	my $caller = shift;

 	my $class = ref($caller) || $caller;
  	my $self = $class->SUPER::new(@_);
	
	my (
   		$load
    )
    = $self->_rearrange( [
      	'LOAD'
    ],
    @_
    );
	
	
	$load && $self->load($load);
	my $region = Bio::Unitrap::Region->new(-LOAD=>$load);
	$region && $self->region($region);
  	return $self;
}

sub load{
  my ($self, $value) = @_;
  $self->{'load'} = $value if defined $value;
  
  return $self->{'load'};
}

sub region{
  my ($self, $value) = @_;
  $self->{'region'} = $value if defined $value;
  
  return $self->{'region'};
}

sub bio_intersection{
	my ($self,$slice,$obj, $f) = @_;
	
	my $start = $slice->start + $obj->start - 1;
	my $end = $start + $obj->length - 1;
	my $display_name;
	if($obj->isa("Bio::EnsEMBL::Intron")){
		$display_name = $obj->prev_Exon->stable_id."-".$obj->next_Exon->stable_id;
	}
	else{
		$display_name = $obj->stable_id;
	}
	
	
	my $feature = Bio::SeqFeature::Generic->new (	-display_name => $display_name,
														-start => $start,
														-end => $end,
														-strand => $obj->strand);
	
	my $i = $f->intersection($feature);	

	return $i;
}

sub intersection{
	my ($self,$obj,$obj_start,$obj_end, $f) = @_;
	my %hash;
	
	my $block_start = $f->start;
	my $block_end = $f->end;
	my $block_length = ($block_end - $block_start + 1);
	my $obj_length = ($obj_end - $obj_start + 1);
	if($obj_start > $block_start && $obj_end > $block_end && $obj_length==$block_length ){$block_length = ($obj_end - $block_start + 1)}
    elsif($obj_start < $block_start && $obj_end < $block_end && $obj_length==$block_length){$block_length = ($block_end - $obj_start + 1)}
    
	my $dE_e = $obj_end - $block_end; 
    my $dS_s = $obj_start - $block_start;
    my $i = 0;
    my $coverage;
	#$debug && print STDERR "INSERSECTION BTW OBJ ";
	#$debug && print STDERR ref $obj;
	#$debug && print STDERR " $obj_start - $obj_end BLOCK $block_start - $block_end\n";
	 #the impossible cases
    if ($obj_start > $obj_end || $block_start > $block_end) {
		print STDERR "line 127 impossible case OBJ $obj_start > $obj_end || BLOCK $block_start > $block_end ";
		print STDERR ref $self;
		print STDERR "->intersection\n";
		exit;
    }	
   	
   	if ($obj_end < $block_start || $block_end < $obj_start){
   	
   		$hash{'i'} = 0;
		return \%hash
   	
   	}
   	
	if ($obj_start == $block_start){
    	#  trapblock s|-----------|e
	    #  obj       S|----------------|E
	    #  obj       S|-------|E
	     #  obj      S|-----------|E
    	if ($obj_end > $block_end){$coverage = $block_length/$block_length * 100}
    	elsif ($obj_end < $block_end){$coverage = $obj_length/$block_length * 100}
    	elsif ($obj_end == $block_end){$coverage = $obj_length/$block_length * 100}
    }
    elsif($obj_end == $block_end){
    	#  trapblock      s|-----------|e
	    #  obj       S|----------------|E
	    #  obj                S|-------|E
    	if ($obj_start < $block_start){$coverage = $block_length/$block_length * 100}
    	elsif ($obj_start > $block_start){$coverage = $obj_length/$block_length * 100}
    }
    else{
    	#  trapblock      s|-----------|e
    	#  obj         S|---------|E
	    #  obj         S|----------------|E
	    #  obj           	S|-------|E
	    #  obj            S|----------------|E
    	if ($obj_length > $block_length){$coverage = $block_length/$block_length * 100}
    	elsif ($obj_length <= $block_length){$coverage = $obj_length/$block_length * 100}
    	else{
    		print STDERR "line 163 odd case OBJ $obj_length BLOCK $block_length\n ";
			print STDERR ref $self;
			print STDERR "->intersection\n";
			exit;
		}
    }
	
	unless (defined $coverage){
		print STDERR "line 171 odd case no coverage detected: OBJ $obj_start > $obj_end || BLOCK $block_start > $block_end ";
		print STDERR ref $self;
		print STDERR "->intersection\n";
		exit;
	}
	
    my $rank = 0;
	#rank the exons, give points for coverage and sharing of splice sites
	$rank  += $coverage;
	if ($dE_e == 0) {$rank += 100;}							
	if ($dS_s == 0) {$rank += 100;}
	$hash{'rank'} = $rank;
    
    
    $hash{'dE_e'} = $dE_e;
    $hash{'dS_s'} = $dS_s;
    $hash{'tb_start'} = $block_start;
    $hash{'tb_end'} = $block_end;
    $hash{'feat_start'}=  $obj_start;
    $hash{'feat_end'} = $obj_end;
	$hash{'coverage'} = $coverage;
	
    
	if ($obj_start == $block_start && $obj_end == $block_end && $dE_e == 0 && $dS_s ==0){
	#  trapblock s|-----------|e
	#  feature   S|-----------|E
	$hash{'label'} = 'perfect';
	$hash{'comment'} = 'the trapblock matchs the feature perfectly';
	$i = 1;
    }elsif ($obj_start >  $block_start && $obj_end < $block_end && $dE_e < 0 && $dS_s > 0){
	#  trapblock s|-----------|e
	#  feature      S|----|E
	$hash{'label'} = 'trapblock_contains_feature';
	$hash{'comment'} = 'the trapblock contains the feature';
	$i = 1;
    } elsif ($obj_start <  $block_start && $obj_end > $block_end && $dE_e > 0 && $dS_s < 0){
	#  block trapblock     s|------|e
	#  obj feature S|-----------------|E
	$hash{'label'} = 'trapblock_contained_in_feature';
	$hash{'comment'} = 'the trapblock is contained in the feature';
	$i = 1;
    }elsif ($obj_start ==  $block_start && $obj_end > $block_end && $dE_e > 0 && $dS_s == 0){
	#  block trapblock s|------------|e
	#  obj feature     S|-----------------|E
	$hash{'label'} = 'share_start_tb_end_whitin';
	$hash{'comment'} = 'the trapblock and the feature share start, but the trapblock end is whitin feature';
	$i = 1;
    }elsif ($obj_start ==  $block_start && $obj_end < $block_end && $dE_e < 0 && $dS_s == 0){
	#  block trapblock s|------------------|e
	#  obj feature     S|-----------|E
	$hash{'label'} = 'share_start_tb_end_extends';
	$hash{'comment'} = 'the trapblock and the feature share start, but the trapblock end extends outside of the feature';
	$i = 1;
    } elsif ($obj_start >  $block_start && $obj_end == $block_end && $dE_e == 0 && $dS_s >  0){
	#  block trapblock s|------------------|e
	#  obj feature            S|-----------|E
	$hash{'label'} = 'share_end_tb_start_whitin';
	$hash{'comment'} = 'the trapblock and the feature share end, but the feature start is whitin the trapblock';
    $i = 1;
	} elsif ($obj_start <  $block_start && $obj_end == $block_end && $dE_e == 0 && $dS_s <  0){
	#  block trapblock       s|------------|e
	#  obj feature    S|-------------------|E
	$hash{'label'} = 'share_end_tb_start_extends';
	$hash{'comment'} = 'the trapblock and the feature share end, but the feature start extends outside the trapblock';
    $i = 1;
	} elsif ($obj_start <  $block_start && $obj_end < $block_end && $dE_e < 0 && $dS_s <  0){
	#  block trapblock       s|------------|e
	#  obj feature  S|----------|E
	$hash{'label'} = 'overlap in tb start';
	$hash{'comment'} = 'the trapblock overlaps the feature at start of trapblock';
    $i = 1;
	} elsif ($obj_start >  $block_start && $obj_end > $block_end && $dE_e > 0 && $dS_s >  0){
	#  block trapblock       s|------------|e
	#  obj feature               S|----------|E
	$hash{'label'} = 'overlap at tb end';
	$hash{'comment'} = 'the trapblock overlaps the feature at end of trapblock';
    $i = 1;
	} elsif ($obj_start >  $block_start && $obj_end < $block_end && $dE_e < 0 && $dS_s >  0){
	#  block trapblock       s|------------------|e
	#  obj feature               S|----------|E
	$hash{'label'} = 'feature_contained_in_trapblock';
	$hash{'comment'} = 'the feature is contained in trapblock';
    $i = 1;
	}elsif ($obj_end == $block_start){
	#  block trapblock            s|------------|e
	#  obj feature   S|----------|E
	$hash{'label'} = 'tb_start_matches_feat_end';
	$hash{'comment'} = 'the trapblock start is the end of the feature';
    $i = 1;
	} elsif ($obj_start == $block_end){
	#  block trapblock s|------------|e
	#  obj feature                S|----------|E
	$hash{'label'} = 'tb_end_matches_feat_start';
	$hash{'comment'} = 'the trapblock end is the start of the feature';
    $i = 1;
	} 
    else {
	print STDERR "\t\tPOSITION: something is fishy  exit\n";
		exit;
    }   
	$hash{'i'} = 1;
	return \%hash;
}

sub find_insertion_site{
	my ($self, $hash) = @_;
	my $fetch = $self->load->fetch;
	my $case;

	my $case1 ="trapped at least one exon and insertion is in upstream intron";
	my $case2 ="trapped at least one exon and insertion is in downstream intron";
	my $case3 ="trapped at least one exon and insertion is in downstream utr";
	my $case4 ="trapped at least one exon and insertion is in upstream utr"; 
	my $case5 ="trapped an intron";
	my $case6 ="trapped an exon";
	my $case7 ="insertion is in 5'UTR";
	my $case8 ="insertion is in 3'UTR";
	
	
	my $region_id = $hash->{'region_id'};
	my $flanking_exon_id = $hash->{'flanking_exon_id'};
	
	my $seqdir = $hash->{'sequencing_direction'};
	my $moltype = $hash->{'moltype'};
		
	my $region = $fetch->get_region_by_id($region_id);
	my $flanking_region = $fetch->get_region_by_id($flanking_exon_id);
	my $strand = $region->{'region_strand'};
	my $gene = $fetch->get_region_by_id($region->{'parent_id'});
	unless ($gene){$gene = $fetch->get_region_by_id($region->{'region_id'});}
	
	my $putative_insertion_start = 0;
	my $putative_insertion_end = 0;
	
	my $rank = $region->{'rank'};
	my $name = $region->{'seq_id'};
	
	
	
	if($region->{'description'} eq "5' UPSTREAM"){
		$debug && print STDERR "\t\t-- Trap outside transcript \n";
		$putative_insertion_start = $region->{'region_start'};
		$putative_insertion_end = $region->{'region_end'};
		$case = 7;
	}
	elsif($region->{'description'} eq "3' UPSTREAM"){
		$debug && print STDERR "\t\t-- Trap outside transcript \n";
		$putative_insertion_start = $region->{'region_start'};
		$putative_insertion_end = $region->{'region_end'};
		$case = 8;
	}
	
	elsif($region->{'description'} =~ /intron/){#GOT AN INTRON;	
		$debug && print STDERR "\t\t-- Trap Intronic \n";
		$putative_insertion_start = $region->{'region_start'};
		$putative_insertion_end = $region->{'region_end'};
		
		$case = 5;
	}
	else{
		if($moltype eq "genomic DNA"){
			$putative_insertion_start = $region->{'region_start'};
			$putative_insertion_end = $region->{'region_end'};
			if ($region->{'description'} =~ /intron/){$case = 5}
			elsif ($region->{'description'} =~ /EXON/){$case = 6}
			elsif ($region->{'description'} =~ /5' UTR/){$case = 7}
			elsif ($region->{'description'} =~ /3' UTR/){$case = 8}
			else{print STDERR " GOT HERE BECAUSE PROBLEM WITH REGION: $region->{'description'}\n"}
		}
		else{
			if ($region_id == $flanking_exon_id){
				if (($strand eq "1" && $seqdir eq "3") || ($strand eq "-1" && $seqdir eq "5")){
					# GOT A 5'UTR;
					$case = 4;
					$putative_insertion_end = $region->{'region_start'} - 1;
					$putative_insertion_start = $region->{'region_start'} - 100;
				}
				if (($strand eq "1" && $seqdir eq "5") || ($strand eq "-1" && $seqdir eq "3")){
					# GOT A 3'UTR;
					$case = 3;
					$putative_insertion_end = $region->{'region_end'} + 100;
					$putative_insertion_start = $region->{'region_end'} + 1;
				}
			}
			else{
				my ($intron,$intron_rank);
				
				if (($strand eq "1" && $seqdir eq "3") || ($strand eq "-1" && $seqdir eq "5")){
					$debug && print STDERR "unitrap in upstream intron\n";
					$case = 1;
					$debug && print STDERR "STRAND $strand  Rank  $rank\t";
					$intron_rank = $region->{'rank'} - 1;				
				}
				if (($strand eq "1" && $seqdir eq "5") || ($strand eq "-1" && $seqdir eq "3")){
					$debug && print STDERR "unitrap in downstream intron\n";
					$case = 2;
					$intron_rank = $region->{'rank'};
					
				}
				if ($flanking_region->{'description'} =~ /intron/){
					$intron = $flanking_region;	
				}
				
				unless (defined $intron){
					$intron = $fetch->get_intron_by_rank_parent_id($intron_rank,$region->{'parent_id'});
					unless (defined $intron){
						if ($case == 1){
							# GOT A 5'UTR;
							$case = 4;
							$putative_insertion_end = $region->{'region_start'} - 1;
							$putative_insertion_start = $gene->{'region_start'};
						}
						if ($case == 2){
							# GOT A 3'UTR;
							$case = 3;
							$putative_insertion_end = $gene->{'region_end'};
							$putative_insertion_start = $region->{'region_end'} + 1;
						}					
					}
					else{
						$name = $intron->{'seq_id'};
						$putative_insertion_end = $intron->{'region_end'};
						$putative_insertion_start = $intron->{'region_start'};
					}
				}
				
				
			}
		}
	}
	
	$debug && print STDERR "NAME $name CASE $case";
	$debug && print STDERR "\n";
	#if ($putative_insertion_end == $putative_insertion_start){$putative_insertion_start = $putative_insertion_end - 100;}

	my $mycase;
	if ($case == 1){$mycase = $case1}
	if ($case == 2){$mycase = $case2}
	if ($case == 3){$mycase = $case3}
	if ($case == 4){$mycase = $case4}
	if ($case == 5){$mycase = $case5}
	if ($case == 6){$mycase = $case6}
	if ($case == 7){$mycase = $case7}
	if ($case == 8){$mycase = $case8}
	if (!$putative_insertion_start && !$putative_insertion_end){

		print STDERR "Error: insertion range not defined  case $mycase\n";
		return undef;
	}
	
	$hash->{'region_start'} = $region->{'region_start'};
	$hash->{'region_end'} = $region->{'region_end'};
	$hash->{'putative_insertion_start'} = $putative_insertion_start;
	$hash->{'putative_insertion_end'} = $putative_insertion_end;
	$hash->{'insertion_case'} = $mycase;
	$hash->{'trapped_exon_rank'} = $rank;
	$hash->{'insertion_ambiguous'} = 1 if $case != 5;
	$hash->{'display_name'} = $name;
	
	$self->load->load_insertion($hash);
	
	return $name;
}

sub _old_find_insertion_site {
	my ($self, $hash) = @_;
	
	my @blocks = @{$hash->{'blocks'}};
	my $exs = $hash->{'exons'};
	my $race = $hash->{'sequencing_direction'};
	#working on transcript instead of gene
	my $gene = $hash->{'transcript'};
	my $moltype = $hash->{'moltype'};
	
	my %insert;
	$insert{'insertion'}{'transcript_id'} = $gene->{'region_id'};
	my $case1 ="trapped_exon and upstream_exon"; # (i.e., polyA trap - 3' race)
	my $case2 ="trapped_exon and downstream_exon"; # (i.e., 5' race with reversed trap)
	my $case3 ="trapped_exon and downstream_exon";# (i.e., 5' race with not-reversed trap)
	
	if ($race eq 'na') {$race = "3";}### Default race is 3'
	$insert {'insertion'}{'putative_insertion_start'} = 0;
	$insert {'insertion'}{'putative_insertion_end'} = 0;
	my $nblocks = scalar (@blocks);		
	my $last_block_num = $nblocks - 1;
	my $trap_strand = $blocks[0]->{'strand'};
	
	my @exons = @{$exs};
	my $exons = scalar (@exons);
	my $gene_strand = $gene->{'region_strand'};
	if ($gene_strand == 1) {@exons = sort {$a->{'rank'} <=> $b->{'rank'}} @exons;} 
	else {@exons = sort {$b->{'rank'} <=> $a->{'rank'}} @exons;}
	
	my $exon_id = 0;
	my $flanking_exon_id = 0; 
	my $tag_point;
	my $case;		
	### Checking for the right trapped exon... according to race type
	### If a trapped exon is not identified, a tag point is defined from the trapblock
	if ($race eq "5" && $gene_strand eq "1") {
		$exon_id = $blocks[$last_block_num]->{'region_id'};
		unless ($exon_id) {
			$insert {'insertion'}{'putative_insertion_start'} = $blocks[$last_block_num]->{'end'};
			$insert{'insertion'}{'trapblock_id'} = $blocks[$last_block_num]->{'trapblock_id'};

			$tag_point = $blocks[$last_block_num]->{'end'};
		}
		if ($gene_strand eq $trap_strand) {
			$case = 2;
		} elsif ($gene_strand ne $trap_strand) {
			### This is more or less the case in which the sequence-tag is not reversed
			### Note this should not happen for TIGEM at this point, since I have chosen 
			### all the trapped genes located on the opposite strand where race=5' and project is TIGEM
			$case = 3;
		}
	} elsif ($race eq "5" && $gene_strand eq "-1") {
		$exon_id = $blocks[0]->{'region_id'};
		unless ($exon_id) {
			$insert {'insertion'}{'putative_insertion_end'} = $blocks[0]->{'start'};
			$insert{'insertion'}{'trapblock_id'} = $blocks[0]->{'trapblock_id'};

			$tag_point = $blocks[0]->{'start'};
		}
		if ($gene_strand eq $trap_strand) {
			$case = 2;
		} elsif ($gene_strand ne $trap_strand) {
			### This is more or less the case in which the sequence-tag is not reversed
			### Note this should not happen for TIGEM at this point, since I have chosen 
			### all the trapped genes located on the opposite strand where race=5' and project is TIGEM
			$case = 3;
		}
	}  elsif ($race eq "3" && $gene_strand eq "-1") {
		$exon_id = $blocks[$last_block_num]->{'region_id'};
		unless ($exon_id) {
			$insert {'insertion'}{'putative_insertion_start'} = $blocks[$last_block_num]->{'end'};
			$insert{'insertion'}{'trapblock_id'} = $blocks[$last_block_num]->{'trapblock_id'};

			$tag_point = $blocks[$last_block_num]->{'end'};
		}
		$case = 1;
	}  elsif ($race eq "3" && $gene_strand eq "1") {
		$exon_id = $blocks[0]->{'region_id'};
		unless ($exon_id) {
			$insert {'insertion'}{'putative_insertion_end'} = $blocks[0]->{'start'};
			$insert{'insertion'}{'trapblock_id'} = $blocks[0]->{'trapblock_id'};
			
			$tag_point = $blocks[0]->{'start'};
		}
		$case = 1;
	}
		
		### Once the trapped exon or the tag_point is identified, we identify the flaking exon according to the case
	if ($case) {
		my $mycase;
		if ($case == 1){$mycase = $case1}
		if ($case == 2){$mycase = $case2}
		if ($case == 3){$mycase = $case3}
		$insert {'insertion'}{'insertion_case'} = $mycase;
		if ($exon_id) {
			#$insert {'insertion'}{'region_id'} = $exon_id;
			
			my (%exon_hash, $trapped_exon_num, $trapped_exon_obj);
			foreach my $exon (@exons) {
				my $i = $exon->{'rank'};
				my $this_exon_id = $exon->{'region_id'};
				my $exon_start = $exon->{'region_start'};
				my $exon_end = $exon->{'region_end'};
				if($exon_start > $exon_end){
					my $tag = $self->_check_region_range($exon, 'EXON',$gene_strand);
					unless ($tag){#print STDERR "EXON START $exon_start > END $exon_end\n"; die;
					}
				}
				if ($this_exon_id eq $exon_id) {
					$trapped_exon_num = $exon->{'rank'};
					$trapped_exon_obj = $exon;
					$insert {'insertion'}{'trapped_exon_rank'} = $trapped_exon_num;
				}
				$exon_hash {$i} {'exon'} = $this_exon_id;
				$exon_hash {$i} {'obj'} = $exon;
				
			}
			if ($trapped_exon_num == 1 && $case == 1) {
				if ($exons == 1) {
				}
				if ($gene_strand == '1') {
					$insert {'insertion'}{'putative_insertion_end'} = $blocks[0]->{'start'};
					$insert{'insertion'}{'trapblock_id'} = $blocks[0]->{'trapblock_id'};
				} elsif ($gene_strand == '-1') {
					$insert {'insertion'}{'putative_insertion_start'} = $blocks[$last_block_num]->{'end'};
					$insert{'insertion'}{'trapblock_id'} = $blocks[$last_block_num]->{'trapblock_id'};
				}
				$insert {'insertion'}{'insertion_ambiguous'} = 1;
			} elsif ($trapped_exon_num == $exons && ($case == 3 || $case == 2)) {
				if ($gene_strand == '1') {
					 $insert {'insertion'}{'putative_insertion_start'} = $blocks[$last_block_num]->{'end'};
					 $insert{'insertion'}{'trapblock_id'} = $blocks[$last_block_num]->{'trapblock_id'};
				} elsif ($gene_strand == '-1') {
					 $insert {'insertion'}{'putative_insertion_end'} = $blocks[0]->{'start'};
					 $insert{'insertion'}{'trapblock_id'} = $blocks[0]->{'trapblock_id'};
				}
				$insert {'insertion'}{'insertion_ambiguous'} = 1;
			} else {
				my ($flanking_exon_id, $flanking_exon_obj);
				# NB:
				# case 1: intron between trapped_exon and upstream_exon
				# case 2: intron between trapped_exon and downstream_exon
				# case 3: intron between trapped_exon and downstream_exon


				### Calculate the flanking_exon
				if ($case == 1) {
					 my $flanking_exon_num = $trapped_exon_num - 1;
					 $flanking_exon_id = $exon_hash {$flanking_exon_num} {'exon'};
					 die if ($exon_id - $flanking_exon_id) != 1;
					 $flanking_exon_obj = $exon_hash {$flanking_exon_num} {'obj'};
				} elsif ($case == 2 || $case == 3) {
					 my $flanking_exon_num = $trapped_exon_num + 1;
					 $flanking_exon_id = $exon_hash {$flanking_exon_num} {'exon'};
					 die if ($flanking_exon_id - $exon_id) != 1;
					 $flanking_exon_obj = $exon_hash {$flanking_exon_num} {'obj'};
				}


				### Handling the intron object
				if ($exon_id && $flanking_exon_id) {
					my $intron;
					if ($case == 1) {
						 $intron = new Bio::EnsEMBL::Intron ($flanking_exon_obj,$trapped_exon_obj);
					} elsif ($case == 2 || $case == 3) {
						 $intron = new Bio::EnsEMBL::Intron ($trapped_exon_obj,$flanking_exon_obj);
					}
					my $intron_start = $intron->start();
					my $intron_end = $intron->end();
					if ($intron_start < $intron_end) {
						 $insert {'insertion'}{'putative_insertion_start'} = $intron_start;
						 $insert {'insertion'}{'putative_insertion_end'} = $intron_end;
						 $insert {'insertion'}{'trapblock_id'} = $blocks[$last_block_num]->{'trapblock_id'};
					} else {
						 $insert {'insertion'}{'putative_insertion_start'} = $intron_end;
						 $insert {'insertion'}{'putative_insertion_end'} = $intron_start;
						 $insert {'insertion'}{'trapblock_id'} = $blocks[0]->{'trapblock_id'};
					}
					$insert {'insertion'}{'insertion_ambiguous'} = 0;
				}
			}
		} elsif ($tag_point) {
			my ($pre_exon_id, $pre_exon_obj, $next_exon_id, $next_exon_obj);
			foreach my $exon (@exons) {
				my $this_exon_id = $exon->{'region_id'};
				my $exon_start = $exon->{'region_start'};
				my $exon_end = $exon->{'region_end'};
				
				if ($gene_strand == '1') {
					if ($tag_point >= $exon_end) {
						$pre_exon_id = $this_exon_id;
						$pre_exon_obj = $exon;
						next;
					} else {
						$next_exon_id = $this_exon_id;
						$next_exon_obj = $exon;
						last;
					}
				} elsif ($gene_strand == '-1') {
					if ($tag_point <= $exon_start) {
						$pre_exon_id = $this_exon_id;
						$pre_exon_obj = $exon;
						next;
					} else {
						$next_exon_id = $this_exon_id;
						$next_exon_obj = $exon;
						last;
					}
				}
			}
			$insert {'insertion'}{'insertion_ambiguous'} = 1;
			if ($pre_exon_id && $next_exon_id) {
				if ($gene_strand == '1') {
					if ($case == '1') {
						$insert {'insertion'}{'putative_insertion_start'} = $pre_exon_obj->end;
						$insert{'insertion'}{'trapblock_id'} = $blocks[0]->{'trapblock_id'};
					} elsif ($case == '2' || $case == '3') {
						$insert {'putative_insertion_end'} = $next_exon_obj->start;
						$insert{'insertion'}{'trapblock_id'} = $blocks[0]->{'trapblock_id'};
					}
				} elsif ($gene_strand == '-1') {
					if ($case == '1') {
						$insert {'insertion'}{'putative_insertion_end'} = $pre_exon_obj->start;
						$insert{'insertion'}{'trapblock_id'} = $blocks[$last_block_num]->{'trapblock_id'};
					} elsif ($case == '2' || $case == '3') {
						$insert {'insertion'}{'putative_insertion_start'} = $next_exon_obj->end;
						$insert{'insertion'}{'trapblock_id'} = $blocks[$last_block_num]->{'trapblock_id'};
					}
				}
			}
			$insert {'insertion'}{'insertion_ambiguous'} = 1;
			$insert {'insertion'}{'label'} = 'novel_exon';


			if ($moltype ne 'mRNA'){
				if ($pre_exon_id) {$insert {'exon_id'} = $pre_exon_id;}			
				if ($next_exon_id) {$insert {'flanking_exon_id'} = $next_exon_id;}		
				if ($pre_exon_id && $next_exon_id) {
					my $intron;
					if ($gene_strand == '1') {$intron = new Bio::EnsEMBL::Intron ($pre_exon_obj,$next_exon_obj);
					}elsif ($gene_strand == '-1') {$intron = new Bio::EnsEMBL::Intron($next_exon_obj,$pre_exon_obj)}
					
					my $intron_start = $intron->start();
					my $intron_end = $intron->end();
					
					if ($intron_start < $intron_end) {
						 $insert {'putative_insertion_start'} = $intron_start;
						 $insert {'putative_insertion_end'} = $intron_end;
						 $insert{'insertion'}{'trapblock_id'} = $blocks[$last_block_num]->{'trapblock_id'};
					} else {
						 $insert {'putative_insertion_start'} = $intron_end;
						 $insert {'putative_insertion_end'} = $intron_start;
						 $insert{'insertion'}{'trapblock_id'} = $blocks[0]->{'trapblock_id'};
					}
					
					$insert {'insertion_ambiguous'} = 0;
					#$debug && print STDERR "The insertion is in the INTRON between $pre_exon_id & $next_exon_id - $intron_start..$intron_end\n";
				} else {
					### We are in the case in which we don't have both pre- and next- exons.
					$insert {'insertion_ambiguous'} = 1;
				}
			}
			
			
		} else {
			#$debug && print STDERR "No exon_id, no tag_point.\n";
			exit;
		}
	} else {
			#$debug && print STDERR "No case!!!\n";
			exit;
	}	
	#$debug && print STDERR "line 1019 ";
	#$debug && print STDERR ref $self;
	#$debug && print STDERR "->find_insertion_site region id ",$insert{'insertion'}{'trapblock_id'},"\n";	
	return (\%insert);
}

sub check_strand {
	my ($self, $trap, $tmf, $gene) = @_;
	my $amb = 0;
	
	my $seqdir = $trap->{'sequencing_direction'};
	my $moltype = $trap->{'mol_type'};
	$debug && print STDERR "$moltype, $seqdir,".$gene->strand.",".$tmf->strand."==>";
	if ($moltype eq "genomic DNA"){
		$amb = 1;
		if ($seqdir eq "F"){
			#true as it is only using mgi data
			if ($gene->strand eq "-1" and $tmf->strand eq "-1"){$amb = 0}
			elsif ($gene->strand eq "1" and $tmf->strand eq "1"){$amb = 0}
		}
		elsif ($seqdir eq "R"){
			if ($gene->strand eq "1" and $tmf->strand eq "1"){$amb = 0}
			elsif ($gene->strand eq "-1" and $tmf->strand eq "-1"){$amb = 0}
		}
		elsif ($seqdir eq "3SPK"){
		
			if ($gene->strand eq $tmf->strand ){$amb = 0}
		}
		elsif ($seqdir eq "5SPK"){
			#true as it is only using mgi data
			if ($gene->strand eq $tmf->strand ){$amb = 0}
		}
		else{
			if ($gene->strand eq "1" and $tmf->strand eq "1"){$amb = 0}
			elsif ($gene->strand eq "-1" and $tmf->strand eq "-1"){$amb = 0}
		}
	}
	elsif ($moltype eq "mRNA"){
		$amb = 1;
		if ($gene->strand eq "1" and $tmf->strand eq "1"){$amb = 0}
		elsif ($gene->strand eq "-1" and $tmf->strand eq "-1"){$amb = 0}
	}
	$debug && print STDERR "amb = $amb\n";
	return $amb;

}

sub define_up_downstream_trapped_region {
	my ($self,$hash,$gene,$tbf,$first,$last,$seqdir,$moltype) = @_;
	
	my $flanking_exon_id;
	my $region_id;
	my $stream_id;
	my $first_exon_id = $first->{'region_id'};
	my $last_exon_id = $last->{'region_id'};
	my $slice = $gene->slice;
	if ($gene->strand eq '1'){
		if ($tbf->end < $first->{'region_start'}){
			# insertion start somewhere in the 5' UPSTREAM 
			my $f = Bio::SeqFeature::Generic->new(	-display_name => "5' UPSTREAM",
								-start => $gene->slice->start + $gene->start,
								-end => $first->{'region_start'} - 1,
								-strand => $gene->strand	);
			$stream_id = $self->region->up_downstream($slice, $gene,$f); #load STREAM in region table
			if ($moltype eq "mRNA"){
				if ($seqdir eq "3") {
					#I cannot define the exon
					$region_id = $stream_id; #trapped region
					#the unitrap is the STREAM;
					$flanking_exon_id = $first_exon_id;
				}
				if($seqdir eq "5") {
					$region_id = $first_exon_id; #exon
					#the unitrap is the  STREAM;
					$flanking_exon_id = $first_exon_id;
				}
			}
			if ($moltype eq "genomic DNA"){
					
				$region_id = $stream_id; #trapped region
				#the unitrap is the STREAM;
				$flanking_exon_id = $first_exon_id;
				
			}
			
		}elsif ($tbf->start > $last->{'region_end'}){
			#insertion start somewhere in the 3' DOWNSTREAM
			my $f = Bio::SeqFeature::Generic->new(	-display_name => "3' DOWNSTREAM",
								-start => $last->{'region_end'} + 1,
								-end => $gene->slice->start + $gene->start + $gene->length -2,
								-strand => $gene->strand	);
			$stream_id = $self->region->up_downstream($slice, $gene,$f); #load STREAM in region table
			if ($moltype eq "mRNA"){
				if ($seqdir eq "5") {
					#I cannot define the exon
					$region_id = $stream_id; #trapped region
					#the unitrap is the STREAM;
					$flanking_exon_id = $last_exon_id;
				}
				if($seqdir eq "3") {
					$region_id = $last_exon_id; #exon
					#the unitrap is the  STREAM;
					$flanking_exon_id = $last_exon_id;
				}
				unless ($region_id || $flanking_exon_id){die "$!\n"}
				$debug && print STDERR "Region id $region_id\t Flanking $flanking_exon_id\n"; 
			}	
			if ($moltype eq "genomic DNA"){
						
				$region_id = $stream_id; #trapped region
				#the unitrap is the STREAM;
				$flanking_exon_id = $last_exon_id;
				
			}
		}
	}
	if ($gene->strand eq '-1'){#assuming that gene start is < gene end
		if ($tbf->end < $first->{'region_start'}){
			# insertion start somewhere in the 3' DOWNSTREAM
			my $f = Bio::SeqFeature::Generic->new(	-display_name => "3' DOWNSTREAM",
								-start => $gene->slice->start + $gene->start,
								-end => $first->{'region_start'} - 1,
								-strand => $gene->strand	);
			$stream_id = $self->region->up_downstream($slice, $gene,$f); #load STREAM in region table
			if ($moltype eq "mRNA"){
				if ($seqdir eq "5") {
					#I cannot define the exon
					$region_id = $stream_id; #trapped region
					#the unitrap is the STREAM;
					$flanking_exon_id = $first_exon_id;
				}
				if($seqdir eq "3") {
					$region_id = $first_exon_id; #exon
					#the unitrap is the  STREAM;
					$flanking_exon_id = $first_exon_id;
				}
				unless ($region_id || $flanking_exon_id){die "$!\n"}
				$debug && print STDERR "Region id $region_id\t Flanking $flanking_exon_id\n"; 
			}
			if ($moltype eq "genomic DNA"){
				$region_id = $stream_id; #trapped region
				#the unitrap is the STREAM;
				$flanking_exon_id = $first_exon_id;
			}
			
		}elsif ($tbf->start > $last->{'region_end'}){
			# insertion start somewhere in the 5' UPSTREAM
			my $f = Bio::SeqFeature::Generic->new(	-display_name => "5' UPSTREAM",
								-start => $last->{'region_end'} + 1,
								-end => $gene->slice->start + $gene->start + $gene->length - 2,
								-strand => $gene->strand	);
			$stream_id = $self->region->up_downstream($slice, $gene,$f); #load STREAM in region table
			if ($moltype eq "mRNA"){
				if ($seqdir eq "3") {
					$debug && print STDERR "I cannot define the exon\n";
					$region_id = $stream_id; #trapped region
					$debug && print STDERR "the 3'RACE unitrap is the STREAM\n";
					$flanking_exon_id = $last_exon_id;
				}
				if($seqdir eq "5") {
					$region_id = $last_exon_id; #exon
					$debug && print STDERR "the 5'RACE unitrap is the  STREAM\n";
					$flanking_exon_id = $last_exon_id;
				}
				unless ($region_id || $flanking_exon_id){die "$!\n"}
				$debug && print STDERR "Region id $region_id\t Flanking $flanking_exon_id\n"; 
			}
			elsif ($moltype eq "genomic DNA"){
				$region_id = $stream_id; #trapped region
				$debug && print STDERR "the genomic unitrap is the STREAM\n";
				$flanking_exon_id = $last_exon_id;				
			}
		}
	}
	
	$hash->{'insertion'}{'region_id'} = $region_id;
	$hash->{'insertion'}{'flanking_exon_id'} = $flanking_exon_id;
	return $hash;
	
}

sub define_internal_trapped_region {
	my ($self,$hash,$gene,$tbf,$first,$last,$prev,$post,$seqdir,$moltype) = @_;
	my $flanking_exon_id;
	my $region_id;
	my $case = 0;
	$debug && print STDERR ref $self;
	$debug && print STDERR "->define_internal_trapped_region\n";
	
	if ($gene->strand eq '1'){
			if ($moltype eq "mRNA"){
				if ($seqdir eq "3") {
					$region_id = $first->{'region_id'};
					$flanking_exon_id = $prev->{'region_id'};
					$case = 1;
				}
				if($seqdir eq "5") {
					$region_id = $last->{'region_id'};
					$flanking_exon_id = $post->{'region_id'};
					$case = 2;
				}
				unless ($region_id || $flanking_exon_id){die "$!\n"}
				$debug && print STDERR "CASE $case: Region id $region_id\t Flanking $flanking_exon_id\n"; 
			}
			if ($moltype eq "genomic DNA"){
				$region_id = $first->{'region_id'};
				$flanking_exon_id = $prev->{'region_id'};
			}
	}
	if ($gene->strand eq '-1'){#assuming that gene start is < gene end
		if ($moltype eq "mRNA"){
			if ($seqdir eq "5") {
				$region_id = $first->{'region_id'};
				$flanking_exon_id = $prev->{'region_id'};
				$case = 3;
			}
			if($seqdir eq "3") {
				$region_id = $last->{'region_id'};
				$flanking_exon_id = $post->{'region_id'};
				$case = 4;
			}
			unless ($region_id || $flanking_exon_id){die "$!\n"}
				$debug && print STDERR "CASE $case: Region id $region_id\t Flanking $flanking_exon_id\n"; 
		}
		if ($moltype eq "genomic DNA"){
			$region_id = $last->{'region_id'};
			$flanking_exon_id = $post->{'region_id'};
		}
	}
	
	$hash->{'insertion'}{'region_id'} = $region_id;
	$hash->{'insertion'}{'flanking_exon_id'} = $flanking_exon_id;
	return $hash;


}



1;