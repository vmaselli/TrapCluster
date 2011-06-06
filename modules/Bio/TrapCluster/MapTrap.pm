#!/usr/bin/perl -w

=head1 Bio::TrapCluster::MapTrap

=head2 Authors

=head3 Created by

             Vincenza Maselli
             v.maselli@ucl.ac.uk

=head3 Original script by

              Guglielmo Roma
              guglielmoroma@libero.it

=head2 Description
        
             This module run blast or read information from blastoutput and load the tables `trapmap` and `trapblock`	
             
=head2 Usage

	    my $obj = Bio::TrapCluster::MapTrap->new;
			my @blasts = @{$obj->run_blast};
			my @trap_maps = @{$obj->map_trap};
            $obj->load_db(\@trap_maps);
	   

	    }
            
=cut

package Bio::TrapCluster::MapTrap;

use strict;
use DBI;
use Carp;
use Data::Dumper;
use vars qw(@ISA);
use File::Spec;
use Bio::TrapCluster::Utils::Argument qw(rearrange);
use Bio::TrapCluster::Utils::Exception qw(throw warning deprecate);
use Bio::Tools::GFF;
@ISA = qw(Bio::Root::Root);

#setting global variables

require "$ENV{'TrapCluster'}/trapcluster_conf.pl";

my %conf =  %::conf;
my $debug = $conf{'global'}{'debug'};
my $debugSQL = $conf{'global'}{'debugSQL'};
my $mysql_path =$conf{'default'}{'mysql_path'};
my $tmpdir = $conf{'default'}{'tmp_dir'};
my $wrapdir = $conf{'application'}{'wrap_dir'};
my $swarm = $conf{'global'}{'swarm'};
my $run_blast = $conf{'application'}{'run_blast'};
my $annotation = $conf{'annotation'}{'do'};
my $mrna_blast = $conf{'application'}{'mrna_blast'};
my $exonerate_opt = $conf{'application'}{'exonerate'};
use Bio::TrapCluster::LoadTrap;
use Bio::TrapCluster::AnnotationTrap;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::TrapCluster::Utils::File;
use Bio::SeqFeature::Generic;
use Bio::TrapCluster::RetrieveTrap;
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
	$genomic_parameters,
	$mrna_parameters,
	$exonerate_parameters
	)
	= $self->_rearrange( [
	  'GENOMIC_PARAM',
	  'MRNA_PARAM',
	  'EXONERATE_PARAM'
	],
	@_
	);
	
	foreach my $key2 (keys %{$mrna_blast}){
		$mrna_parameters .= "-".$key2." ".$mrna_blast->{$key2}." ";
	}
	
	foreach my $key (keys %{$exonerate_opt}){
		next if $key eq 'ryo';
		$exonerate_parameters .= "--".$key." ".$exonerate_opt->{$key}." ";
	}	
	
	$self->mrna_param($mrna_parameters);
	$self->exonerate_param($exonerate_parameters);
	my $load = Bio::TrapCluster::LoadTrap->new;
	$self->load($load);
	my $ann = Bio::TrapCluster::AnnotationTrap->new;
	$self->annotation($ann);
	my $retrieve = Bio::TrapCluster::RetrieveTrap->new;
	$self->retrieve($retrieve);
	
	return $self;
}

# sub subname{
# 
#   my ($self, $value) = @_;
#   $self->{'subname'} = $value if defined $value;
#   
#   return $self->{'subname'};
# }

sub genomic_param{

  	my ($self, $value) = @_;
  	$self->{'genomic_param'} = $value if defined $value;
  
  	return $self->{'genomic_param'};
}

sub mrna_param{

  	my ($self, $value) = @_;
  	$self->{'mrna_param'} = $value if defined $value;
  
  	return $self->{'mrna_param'};
}

sub exonerate_param{

  	my ($self, $value) = @_;
  	$self->{'exonerate_param'} = $value if defined $value;
  
  	return $self->{'exonerate_param'};
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

=head2 retrieve

 Title    : retrieve
 Usage    : $obj->retrieve([$newval])
 Function : get/set method for attribute retrieve 
 Returns  : value of retrieve
 Args     : newval of retrieve (optional)

=cut 

sub retrieve{

  my ($self, $value) = @_;
  $self->{'retrieve'} = $value if defined $value;
  
  return $self->{'retrieve'};
}

=head2 annotation

 Title    : annotation
 Usage    : $obj->annotation([$newval])
 Function : get/set method for attribute annotation 
 Returns  : value of annotation
 Args     : newval of annotation (optional)

=cut 

sub annotation{

  my ($self, $value) = @_;
  $self->{'annotation'} = $value if defined $value;
  
  return $self->{'annotation'};
}

sub run_blast{
  	my ($self, $item,$type) = @_;
 	
 	my $fasta;
 	my $genome = File::Spec->catfile($ENV{'Blastdir'},$conf{'default'}{'blastdb'});
 	my $parameters = $self->mrna_param;
 	
 	if ($type eq 'trap'){
 		$fasta = Bio::TrapCluster::Utils::TextUtility->get_trap_by_id_fasta($item->{'trap_id'});	
 	}
 	elsif ($type eq 'fasta'){
 		$fasta = $item;
 	}
 	my $output = $fasta;
 	$output =~ s/\.fa/.bl/;
	unless ($output =~ /tmp/){$output = "/tmp/".$output}
	#$genome = "/lustre/scratch103/sanger/vvi/ncbi_m37_blast_index_no_NT/mus_musculus.ncbim37.v59.chromosomal.dna.fa";
	#$fasta = "$ENV{'HOME'}/workdir/11112010/244230.fa";
	#$output = $fasta;
	#$output =~ s/fa/bl/;
	my $timestmp = `date`;
	my $command = "/software/pubseq/bin/ncbi_blast_2.2.21/blastall -p blastn $parameters -d $genome -i $fasta -o $output";
	$debug && print STDERR "BLAST START AT $timestmp $command\n";
	
	`$command`;
	
	$timestmp = `date`;
	$command = "/software/pubseq/bin/ncbi_blast_2.2.21/blastall -p blastn $parameters -d $genome -i $fasta -o $output";
	$debug && print STDERR  "BLAST END AT $timestmp\n";

   	$self->{'run_blast'} = $output;
   	
   	return $self->{'run_blast'};

}

sub run_exonerate{
  	my ($self, $item,$type) = @_;
 	
 	my $fasta;
 	my $genome = File::Spec->catfile($ENV{'Blastdir'},$conf{'default'}{'blastdb'});
 	my $parameters = $self->exonerate_param;
 	
 	if ($type eq 'trap'){
 		$fasta = Bio::TrapCluster::Utils::TextUtility->get_trap_by_id_fasta($item->{'trap_id'});	
 	}
 	elsif ($type eq 'fasta'){
 		$fasta = $item;
 	}
 	my $results;
	my $ryo = $conf{'application'}{'exonerate'}{'ryo'};
	
	my $timestmp = `date`;
	my $command = "exonerate --ryo $ryo --query $fasta --target $genome $parameters";
	$debug && print STDERR  "EXONERATE START AT $timestmp $command\n";
	my $exon_fh;
	
	open ($exon_fh, $command."|") or die "cant read exonerate output\n";
	
	$timestmp = `date`;
	
	while(<$exon_fh>){
		chomp;
		next unless ($_ =~ /^>/);
		my ($q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end, $t_strand, $score, $q_length, $identical_bases, $percent_id, $score_2) = split/\s+/;
		my $hit = {
		   q_id => $q_id,
		   q_start => $q_start,
		   q_end => $q_end,
		   q_strand => $q_strand,
		   q_length => $q_length,
		   identical_bases => $identical_bases,
		   frac_identical => $identical_bases / $q_length,
		   t_id => $t_id,
		   t_start => $t_start,
		   t_end => $t_end,
		   t_strand => $t_strand,
		   score => $score,
		   identical_bases => $identical_bases,
		   percent_id => $percent_id
		};
		push @{$results->{$q_id}}, $hit;
	}
		
   	$self->{'run_exonerate'} = $results;
   	
   	return $self->{'run_exonerate'};

}

sub pars_exonerate_results {
	my ($self, $results) = @_;
	
	foreach my $query_name (keys %$results){
		$debug && print STDERR  "$query_name:\n";
		my $trap = $self->_get_trap($query_name);
		next unless $trap; #trap already mapped
        my $current = `date`;
		$debug && print STDERR  "$current -- trap not mapped yet \n";
		my $trap_id = $trap->{'trap_id'};
		my $genomic = 0;
		die unless $trap->{'mol_type'};
		if ($trap->{'mol_type'} eq 'genomic DNA'){$genomic = 1}
		
		my $chosen_hit;
		my $prev_id;
		my $qlen = $trap->{'seq_length'};
		unless($qlen){die;}
		my $min = $qlen * 0.7;
		$current = `date`;
		$debug && print STDERR  "$current -- start parsing hit \n";
		my $tot = scalar @{$results->{$query_name}};
		my $count = 0;
		foreach my $hit (@{$results->{$query_name}}){
			$count ++;
			
			$current = `date`;
			$debug && print STDERR  "$current -- looking at hit $count/$tot\n";
			
			my $t_start = $hit->{t_start};
			my $t_end = $hit->{t_end};
			my $hit_len = $t_end - $t_start + 1;
			
			$current = `date`;
			$debug && print STDERR  "$current -- start to laod hit\n";
			my $trapmap_id = $self->load->load_exonerate_trapmap($hit,$trap_id);
  			my $flag = 0;
  			
  			if ($hit->{frac_identical} < 0.90) {
  				next;
  			}
			if ($hit_len < $min){
				next;
			}
			
			$current = `date`;
			$debug && print STDERR  "$current -- calling _set_chosen_exonerate_hit\n";
			$flag = $self->_set_chosen_exonerate_hit($hit,$chosen_hit,$trapmap_id,$prev_id,$trap->{'mol_type'}, $genomic,$qlen);
			$current = `date`;
			$debug && print STDERR  "$current -- after _set_chosen_exonerate_hit\n";

			if ($flag == 1){last}
			$chosen_hit = $hit if $flag == 0;
			
			$hit->{frac_aligned_query} = $qlen/$hit_len;
			$current = `date`;
			$debug && print STDERR  "$current -- start to laod hps\n";
			my $trapblocks = $self->load->load_exonerate_trapblock($hit, $trapmap_id, $genomic,$qlen);
			
			$prev_id = $trapmap_id;
			
			#$debug && print STDERR  "-- $t_id, $t_start, $t_end, $t_strand, $identical_bases, $q_length, $identical_bases, $frac_identical, $percent_id, $score\n";
	    }
	    
	    if ($annotation){
			$current = `date`;
			$debug && print STDERR  "$current - calling anntrap\n";
			system("perl $ENV{'TrapCluster'}/scripts/anntrap.pl -id $trap_id -idcount 1");
			$debug && print STDERR  "$current - finished anntrap call\n";
		}
	}

}

sub parse_results {

	my ($self, $file, $format) = @_;		
	my $searchio;
	die unless defined $format;
	my $error;
	$debug && print STDERR  "check for $format results on $file\n";
	
	if ($format eq 'Blast') {
		$self->_parse_blast($file, $format);
	}
	elsif($format eq 'psl'){
		
		$self->_parse_blat($file, $format);
	}
	elsif($format =~ /gff/){
		
		$self->_parse_gff($file,$format);
	}
	
	
	return 1;
}

sub _parse_gff {
	my ($self, $file, $format) = @_;	
	my $version = substr($format,-1);
	# specify input via -fh or -file
	my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => $version);
	my $feature;
	# loop over the input stream
	while($feature = $gffio->next_feature()) {

		$self->parse_gff_feat($feature,$file);
	}
	
	$gffio->close();
	
}

sub parse_gff_feat{
	my ($self, $feature, $file) = @_;
	my $type = $feature->primary_tag;	#'Gene'
	my $source = $feature->source_tag;  #'Gene Trap'
	my $location = $feature->location;
	my $start = $location->start;
	my $end = $location->end;
	my $strand = $location->strand;
	my $chr = $feature->seq_id;
	my $score = $feature->score;
	
	my @tags = $feature->get_all_tags;
	my ($trap_name) = $feature->get_tag_values('Sequence_tag_ID'); 
	print "TRAP NAME in _parse_gff $trap_name\n";
	#$trap_name =~ s/\.//;
			
	my $trap = $self->load->fetch->get_trap_by_name($trap_name);
	my $trap_id = $trap->{'trap_id'};
	return 0 unless $trap_id; 
	foreach my $tag (@tags){
		#"Creator TIGM; Sequence_tag_ID IST11440A7BBF1; GenBank_ID EF833342; DBxref MGI:4005427; Type DNA";
		my $value = $feature->get_tag_values($tag);
		
	}
	$feature->add_tag_value('strand',$strand); 
	$feature->display_name($trap_id);
	$feature->add_tag_value('score',$score);
	$feature->add_tag_value('num_hsps',1);
	
	#print Dumper $feature;
	
	my $analysis;
	$analysis->{'analysis'}->{'name'} = "local alignment";
	$analysis->{'analysis'}->{'db'} = $feature->get_tag_values('DBxref') if $feature->has_tag('DBxref');
	$analysis->{'analysis'}->{'db_version'} = $conf{'default'}{'hit_db'};
	$analysis->{'analysis'}->{'db_file'} = $file;
	$analysis->{'analysis'}->{'program'} = 'BLAT';
	$analysis->{'analysis'}->{'program_version'} = 'MGI';
	$analysis->{'analysis'}->{'description'} = "generated by ".ref $self;
	
	my $trapmap_id = $self->load->load_trapmap($feature,$analysis);
	$self->load->update_trapcheck(1,1, $trap_id);
	$feature->add_tag_value('accepted',1);
	$feature->display_name($trapmap_id);
	$feature->add_tag_value('display_name',$trapmap_id);
	#print Dumper $feature;
	
	$self->load->load_trapblock($feature);
	my $check_annotation = $self->load->fetch->get_trapbmap_region_by_trapmap_id($trapmap_id);
	if ($annotation && !$check_annotation){
		my $current = `date`;
		$debug && print STDERR "$current - calling anntrap\n";
		$self->annotation->annotate($trap_id,1);
		$debug && print STDERR "$current - finished anntrap call\n";
	}
	
	return 1;
		
}

sub _parse_blat {
	my ($self, $file, $format) = @_;		
	
	my $searchio;
	die unless defined $format;
	my $error;
	
	eval{$searchio = Bio::SearchIO->new (-file => $file, -format => $format) || die $!;};	
	if ($@) {
		$debug && print STDERR "Could not create Bio::SearchIO object: $@\n";
		$error = $@;
	}
    
    my $current = `date`;
	$debug && print STDERR "$current -- opened blatfile\n";
	while (my $result = $searchio->next_result) {
		#some checks
		
		$current = `date`;
		
		$debug && print   STDERR "$current -- result: ".$result->query_name."\n";
		my $trap = $self->_get_trap($result->query_name);
		next unless $trap; 
		unless ($trap->{'mol_type'}){print STDERR "Molecular type error, exit program\n"; die;}
		my $qlen = $trap->{'seq_length'};
		
        $current = `date`;
		$debug && print STDERR "$current -- trap not mapped yet \n";
		my $trap_id = $trap->{'trap_id'};
		
		my $genomic = 0;
		if ($trap->{'mol_type'} eq 'genomic DNA'){$genomic = 1}
		
		my $chosen_hit;
		my $prev_id;
		
		while (my $hit = $result->next_hit) {
			$current = `date`;
			$debug && print STDERR "$current -- through hit \n"; 
			
            while (my $hsp = $hit->next_hsp){
            	my $strand = $searchio->strand; #modified psl.pm in order to obtain it
            	my @blocks = @{$hsp->gap_blocks("hit")};
				my $num_hsps = scalar @blocks;
				my $start = $hsp->feature2->location->start;
				my $end = $hsp->feature2->location->end;
				
				my $hsp_len = $end - $start + 1;
				
				$debug && print STDERR "NC = ".$hsp->num_conserved."\n";
				$debug && print STDERR "NI =".$hsp->num_identical."\n";
				
				my $mismatches = $hsp->mismatches;	
				my $frac_identical = ($hsp->num_conserved - $hsp->mismatches)/$result->query_length;
				my $gaps = $hsp->gaps('total');
				my $frac_aligned_query=$hsp->num_conserved/$result->query_length;
				
				$current = `date`;
				$debug && printf STDERR "$current -- hit: START %d END %d STRAND %s MISM %d GAP %d\n",$start,$end, $strand,$mismatches, $gaps;
				
				my $feature = Bio::SeqFeature::Generic->new( 
															-start        => $start, 
															-end          => $end,
															-strand       => $strand, 
															-primary_tag      => $hit->name, # -primary_tag is a synonym
															-seq_id			=> $hit->name,
															-source_tag   => 'blat',
															-display_name => $trap_id,
															-score        => $hsp->score,
															-tag          => {num_hsps => $num_hsps,
																			  frac_aligned_query => $frac_aligned_query,
																			  frac_identical => $frac_identical
																			 
																			 } );
						
  				my $trapmap_id = $self->load->load_trapmap($feature);
  				$self->_set_chosen_blocks($hsp, $strand,$trapmap_id,$genomic,$qlen);
  				
				my $flag = 0;
				$current = `date`;
				$debug && print STDERR "$current -- calling _set_chosen_hit\n";
				$flag = $self->_set_chosen_blat_hit($hsp,$chosen_hit,$trapmap_id,$prev_id,$trap->{'mol_type'}, $genomic,$qlen);
				$current = `date`;
				$debug && print STDERR "$current -- after _set_chosen_hit\n";
	
				if ($flag == 1){last}
				$chosen_hit = $hsp if $flag == 0;
				
				$prev_id = $trapmap_id;
				
			}
			
		}

  		if ($annotation){
            $current = `date`;
            $debug && print STDERR "$current - calling anntrap\n";
			$self->annotation->annotate($trap_id,1);
            $debug && print STDERR "$current - finished anntrap call\n";
 		}

	}

}


sub _parse_blast {
	my ($self, $blastfile, $format) = @_;		
	
	my $searchio;
	die unless defined $format;
	my $error;
	open (IN, $blastfile);
	while (<IN>){
			if ($_ =~ /No hit/){
				$debug && print STDERR  "No results, go to the next file\n";
				return 0;
			}
		}
	close(IN);
	
	eval{$searchio = Bio::SearchIO->new (-file => $blastfile, -format => $format) || die $!;};	
	if ($@) {
		$debug && print STDERR "Could not create Bio::SearchIO object: $@\n";
		$error = $@;
	}
    
    my $current = `date`;
	$debug && print STDERR "$current -- opened blastfile\n";
	while (my $result = $searchio->next_result) {
       
		
		#some checks
		$current = `date`;
		$debug && print STDERR "$current -- result: ".$result->query_name."\n";
		my $trap = $self->_get_trap($result->query_name);
		next unless $trap; #trap already mapped
		unless ($trap->{'mol_type'}){print STDERR "Molecular type error, exit program\n"; die;}
		
		my $qlen = $trap->{'seq_length'};
		unless($qlen){print STDERR "Sel length error, exit program\n"; die;}
        $current = `date`;
		$debug && print STDERR "$current -- trap not mapped yet \n";
		my $trap_id = $trap->{'trap_id'};
		
		my $genomic = 0;
		if ($trap->{'mol_type'} eq 'genomic DNA'){$genomic = 1}
		
		my $chosen_hit;
		my $prev_id;
		
		while (my $hit = $result->next_hit) {
            
            #check hsps
			next unless defined $hit->num_hsps;
			if (($hit->significance) && ($hit->significance > 1e-05)) { next;}
			next if $hit->num_hsps > 10;
			next if $hit->num_hsps == 0;
			
			$current = `date`;
			$debug && printf STDERR "$current -- hit: START %d END %d n hsp %d, frac ident %.2f, frac_aligned_query %.2f\n",$hit->start('hit'),$hit->end('hit'),$hit->num_hsps,$hit->frac_identical,$hit->frac_aligned_query;
			
	
			my $feature = Bio::SeqFeature::Generic->new( 
															-start        => $hit->start("hit"), 
															-end          => $hit->end("hit"),
															-strand       => $hit->strand, 
															-primary_tag      => $hit->name, # -primary_tag is a synonym
															-source_tag   => 'blast',
															-display_name => $trap_id,
															-score        => $hit->score,
															-tag          => { frac_aligned_query => $hit->frac_identical,
																			   significance => $hit->significance } );
						
  			my $trapmap_id = $self->load->load_trapmap($feature);
  			$self->_set_chosen_hsps($hit, $trapmap_id,$genomic,$qlen);
  			
  			my $flag = 0;
			$current = `date`;
			$debug && print STDERR "$current -- calling _set_chosen_hit\n";
			$flag = $self->_set_chosen_hit($hit,$chosen_hit,$trapmap_id,$prev_id,$trap->{'mol_type'}, $genomic,$qlen);
			$current = `date`;
			$debug && print STDERR "$current -- after _set_chosen_hit\n";

			if ($flag == 1){last}
			$chosen_hit = $hit if $flag == 0;
			
			$current = `date`;
			$debug && print STDERR "$current -- calling _load_trapblock\n";
			
			$current = `date`;
			
			$prev_id = $trapmap_id;
			
  		}
  		if ($annotation){
            $current = `date`;
            $debug && print STDERR "$current - calling anntrap\n";
			$self->annotation->annotate($trap_id,1);
            $debug && print STDERR "$current - finished anntrap call\n";
 		}

	}

}


sub _set_chosen_hsps{
	my ($self,$hit, $trapmap_id, $genomic,$qlen) = @_;
	unless ($trapmap_id){print STDERR $hit->name."not stored, exit program\n";die;}

	my @blocks;
	
  	my $hit_start = $hit->start("hit"); #start della hit 
  	my $hit_end = $hit->end("hit"); #end della hit
	$debug && print "INIT RANGE $hit_start $hit_end\n";
	my $new_end;
	my $new_start;
	my $nhsp = $hit->num_hsps; 
	my $genomic_start;
	my $genomic_end;
	my $max_percid = 0;
	my $count = 1;
	my $tag = 0;
	my $accepted = 1;
	if ($hit->frac_identical < 0.90) {
		$accepted = 0;
	}
	if ($hit->frac_aligned_query < 0.7 && $hit->num_hsps > 1){
		$accepted = 0;
	}
	if ($hit->frac_aligned_query < 0.5 ){
		$accepted = 0;
	} 
	
	my $prev_strand;
	my $skip = 0;
	my @hsps;

	while (my $hsp = $hit->next_hsp){
		my $start = $hsp->start("hsp");
		my $end = $hsp->end("hsp");
		my $feature = Bio::SeqFeature::Generic->new( 
													-start        => $start, 
													-end          => $end,
													-strand       => $hsp->strand, 
													-primary_tag      => $hit->name, # -primary_tag is a synonym
													-source_tag   => 'blast',
													-display_name => $trapmap_id,
													-score        => $hit->score,
													-tag          => { percent_identity => $hsp->percent_identity,
																	   accepted => $accepted } );
						
  			
	
		unless ($accepted){$self->load->load_trapblock($feature);next;}
	
		my $min = $qlen * 0.3;	
		my $perc_id = $hsp->percent_identity;
		if ($hsp->length < $min){
			$debug && print STDERR "hsp too short, not acceptded\n";
			$accepted = 0;
			$feature->remove_tag('accepted');
			$feature->add_tag_value('accepted',$accepted);
			$self->load->load_trapblock($feature);
			next;
		}
		elsif ($perc_id < 95){
			$debug && print STDERR "hsp has a percentage identity lt 95, not acceptded\n";
			$accepted = 0;
			$feature->remove_tag('accepted');
			$feature->add_tag_value('accepted',$accepted);
			$self->load->load_trapblock($feature);
			next;		
		}
		
		$self->load->load_trapblock($feature);
		
		unless($tag){
			$debug && print STDERR "initialize new range for the hit\n";
			$new_start = $start;
			$new_end = $end;
			$tag = 1;
		}
		else{
			if ($start < $new_start){$new_start = $start}
			elsif ($end > $new_end){$new_end = $end}
			$tag = 1;
		}
		$count ++;
		
		if ($genomic){		
			if ($perc_id > $max_percid){
			  $genomic_start =  $start;
			  $genomic_end =  $end;
			}  
		
			my $diff = ($end - $hit_start + 1) - ($hit_end - $hit_start +1); 
		  	if ($diff && $diff > 100){
		  		#new range defined";
				$new_end = $genomic_end;
				$new_start = $genomic_start;
		  	}
		  $debug && print STDERR "\n";
		}
		
		$hit_end = $new_end ; #definisco una nuova fine
		$hit_start = $new_start; #definisco un nuovo inizio
		
	}
	
	$debug && print STDERR "DATE RANGE $hit_start $hit_end $accepted\n";
	my $update = qq{UPDATE trapmap SET start = $hit_start, end = $hit_end WHERE trapmap_id = $trapmap_id};
	$self->fetch->update($update);
	return 1;
}

sub _set_chosen_blocks{
	my ($self,$hsp,$strand,$trapmap_id, $genomic,$qlen) = @_;
	
	my $accepted = 1;

	foreach my $block (@{$hsp->gap_blocks("hit")}) {
		my $start = $block->[0];
		my $end = $block->[0]+ $block->[1]-1;
		my $perc_id = "NULL";

		my $feature = Bio::SeqFeature::Generic->new( 
												-start        => $start, 
												-end          => $end,
												-strand       => $strand, 
												-primary_tag      => 'hsp', # -primary_tag is a synonym
												-source_tag   => 'blat',
												-display_name => $trapmap_id,
												-score        => $hsp->score,
												-tag          => { percent_identity => $perc_id,
																   accepted => $accepted } );  			

		# unless ($accepted){$self->load->load_trapblock($feature);next;}
# 	
# 		my $min = $qlen * 0.3;	
# 		if ($hsp->length < $min){
# 			$debug && print STDERR "hsp too short, not acceptded\n";
# 			$accepted = 0;
# 			$feature->remove_tag('accepted');
# 			$feature->add_tag_value('accepted',$accepted);
# 			$self->load->load_trapblock($feature);
# 			next;
# 		}
	
	
		$self->load->load_trapblock($feature);
	}

		
	return 1;
}

sub update_trapmap {
	my ($self, $chosen,$multiple,$chosen_filter,$checked, $trapmap_id) = @_;
	my $update = qq{UPDATE trapmap SET `chosen` = $chosen, `multiple` = $multiple, `chosen_filter` = \"$chosen_filter\", `checked` = $checked WHERE trapmap_id = $trapmap_id};
	$debug && print STDERR  "$update\n";
	$self->load->fetch->update($update);
	return $trapmap_id;

}

sub _set_chosen_hit {
	my ($self,$hit,$chosen_hit,$trapmap_id,$prev_id,$moltype, $genomic, $qlen) = @_;
	
	my $chosen = 0;
	my $multiple = 0;
	my $chosen_filter;
	my $checked = 0;
	my $logic_name;
	$debug && print STDERR  "Checking $hit,$chosen_hit,$trapmap_id,$prev_id,$moltype, $genomic, $qlen\n";
	if ($chosen_hit) {
		$debug && print STDERR  "CHOSEN HIT\n";
		if ($hit->frac_identical < $chosen_hit->frac_identical) {
			$logic_name = "following hit with less frac_identical";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			#the previous hit is better
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $trapmap_id);
			$self->load->load_trapblock($chosen_hit, $prev_id, $genomic, $qlen);
			return 1;
		} #discard the hit with less score -- last if wu-blast, next if blast
		elsif ($chosen_hit->num_hsps == 1 && $hit->num_hsps > 1) {
			$chosen_hit = $hit;
			$logic_name = "following hit with more hsps";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			
			if ($moltype eq 'mRNA'){
				$logic_name = "probable pseudogene";
				$debug && print STDERR $logic_name."\tHIT:".$hit->num_hsps."\tCHOSEN:".$chosen_hit->num_hsps."\n";
				my $new_chosen_filter = $self->load->load_logic_name($logic_name);
				$self->update_trapmap(0,$multiple,$new_chosen_filter,$checked, $prev_id);
				$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
				#the previous hit could be a pseudogene - not chosen
				$self->load->load_trapblock($hit, $trapmap_id, $genomic, $qlen);
				return 1;
			}
			
			if ($moltype ne 'mRNA'){
				$logic_name = "to check if pseudogene";
				$multiple = 1;
				$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
				#the previuos hit could be better, chosen both
				$self->load->load_trapblock($chosen_hit, $prev_id, $genomic, $qlen);
				$self->load->load_trapblock($hit, $trapmap_id, $genomic, $qlen);
				return 0;
			}
			
		} 
		elsif ($chosen_hit->num_hsps > 1 && ($hit->num_hsps == 1)) { 
			$multiple = 1;
			$logic_name = "following hit with one hsp";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			if ($moltype eq 'mRNA'){
				#the previuos hit is better
				$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $trapmap_id);
				$self->load->load_trapblock($chosen_hit, $prev_id, $genomic, $qlen);
			}
			if ($moltype ne 'mRNA'){
				$multiple = 1;
				$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
				$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $prev_id);
				#this hit is better -- previous discarded
				$self->load->load_trapblock($hit, $trapmap_id, $genomic, $qlen);
			}
			return 1;
		}
		elsif ($hit->significance > $chosen_hit->significance) {
			$logic_name = "following hit with less significance";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $trapmap_id);
			#the previuos hit is better
			
			$self->load->load_trapblock($chosen_hit, $prev_id, $genomic, $qlen);
			return 1;
		} 
		elsif ($hit->significance < $chosen_hit->significance) {
			$logic_name = "following hit with more significance";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			$chosen_hit = $hit;
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $prev_id);
			$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
			#this hit is better
			$self->load->load_trapblock($hit, $trapmap_id, $genomic, $qlen);
			return 1;
		} 
		elsif ($hit->significance == $chosen_hit->significance) {
			$logic_name = "following hit with equal significance";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			$multiple = 1;
			$chosen = 1;
			my $logic_name = "multiple mapping";
			my $chosen_filter = $self->load->load_logic_name($logic_name);
			$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
			#both are good
			$self->load->load_trapblock($chosen_hit, $prev_id, $genomic, $qlen);
			$self->load->load_trapblock($hit, $trapmap_id, $genomic, $qlen);
			return 0;
		}
		else{
			#$debug && print STDERR "SOMETHING ELSE...\n";
			die;
		}
	} 
	else {
		$debug && print STDERR  "NOT CHOSEN HIT\n";
		$chosen_hit = $hit;
		$chosen=1;
		$logic_name = "first choice";
		$chosen_filter = $self->load->load_logic_name($logic_name);
		$debug && print STDERR  "UPDATING TRAPMAP\n";
		$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked,$trapmap_id);
		$prev_id = $trapmap_id;
		return 0;
	}
	return 0;
}

sub _set_chosen_blat_hit {
	my ($self,$hit,$chosen_hit,$trapmap_id,$prev_id,$moltype, $genomic, $qlen) = @_;
	
	my $chosen = 0;
	my $multiple = 0;
	my $chosen_filter;
	my $checked = 0;
	my $logic_name;
	my @blocks = @{$hit->gap_blocks("hit")};
	my $num_hsps = scalar @blocks;
	
	my $frac_identical = ($hit->num_conserved - $hit->mismatches)/$qlen;
	my $frac_aligned_query=$hit->num_conserved/$qlen;
	if ($frac_identical < 0.90) {
		$logic_name = "frac_identical is too low";
		$chosen_filter = $self->load->load_logic_name($logic_name);
		if($frac_aligned_query < 0.7 && $num_hsps > 1){
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $trapmap_id);
			return 1;
		}
		elsif($frac_aligned_query < 0.5 && $num_hsps > 1){
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $trapmap_id);
			return 1;
		}
		else{
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $trapmap_id);
			return 1;
		}
	}
	
	if ($chosen_hit) {
		my @chosen_blocks = @{$chosen_hit->gap_blocks("hit")};
		my $chosen_hit_num_hsps = scalar @chosen_blocks;
		my $chosen_hit_frac_identical = ($chosen_hit->num_conserved - $chosen_hit->mismatches)/$qlen;
		my $chosen_hit_frac_aligned_query=$chosen_hit->num_conserved/$qlen;
		$debug && print STDERR  "CHOSEN HIT\n";
		if ($frac_identical < $chosen_hit_frac_identical) {
			$logic_name = "following hit with less frac_identical";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			#the previous hit is better
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $trapmap_id);
			return 1;
		} #discard the hit with less score -- last if wu-blast, next if blast
		elsif ($chosen_hit_num_hsps == 1 && $num_hsps > 1) {
			$chosen_hit = $hit;
			$logic_name = "following hit with more hsps";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			
			if ($moltype eq 'mRNA'){
				$logic_name = "probable pseudogene";
				$debug && print STDERR $logic_name."\tHIT:".$num_hsps."\tCHOSEN:".$chosen_hit_num_hsps."\n";
				my $new_chosen_filter = $self->load->load_logic_name($logic_name);
				$self->update_trapmap(0,$multiple,$new_chosen_filter,$checked, $prev_id);
				$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
				#the previous hit could be a pseudogene - not chosen
				return 1;
			}
			
			if ($moltype ne 'mRNA'){
				$logic_name = "to check if pseudogene";
				$multiple = 1;
				$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
				#the previuos hit could be better, chosen both
				return 0;
			}
			
		} 
		elsif ($chosen_hit_num_hsps > 1 && ($num_hsps == 1)) { 
			$multiple = 1;
			$logic_name = "following hit with one hsp";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			if ($moltype eq 'mRNA'){
				#the previuos hit is better
				$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $trapmap_id);
			}
			if ($moltype ne 'mRNA'){
				$multiple = 1;
				$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
				$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $prev_id);
				#this hit is better -- previous discarded
			}
			return 1;
		}
		elsif ($hit->score > $chosen_hit->score) {
			$logic_name = "following hit with less score";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $trapmap_id);
			#the previuos hit is better
			
			return 1;
		} 
		elsif ($hit->score < $chosen_hit->score) {
			$logic_name = "following hit with more score";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			$chosen_hit = $hit;
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $prev_id);
			$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
			#this hit is better
			return 1;
		} 
		elsif ($hit->score == $chosen_hit->score) {
			$logic_name = "following hit with equal score";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			$multiple = 1;
			$chosen = 1;
			my $logic_name = "multiple mapping";
			my $chosen_filter = $self->load->load_logic_name($logic_name);
			$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
			#both are good
			return 0;
		}
		else{
			#$debug && print STDERR "SOMETHING ELSE...\n";
			die;
		}
	} 
	else {
		$debug && print STDERR  "NOT CHOSEN HIT\n";
		$chosen_hit = $hit;
		$chosen=1;
		$logic_name = "first choice";
		$chosen_filter = $self->load->load_logic_name($logic_name);
		$debug && print STDERR  "UPDATING TRAPMAP\n";
		$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked,$trapmap_id);
		$prev_id = $trapmap_id;
		return 0;
	}
	return 0;
}

sub _set_chosen_exonerate_hit {
	my ($self,$hit,$chosen_hit,$trapmap_id,$prev_id,$moltype, $genomic, $qlen) = @_;
	
	my $chosen = 0;
	my $multiple = 0;
	my $chosen_filter;
	my $checked = 0;
	my $logic_name;
	 my $t_id = $hit->{t_id};
    my $t_start = $hit->{t_start};
    my $t_end = $hit->{t_end};
    my $t_strand = $hit->{t_strand};
    my $q_length = $hit->{q_length};
    my $identical_bases = $hit->{identical_bases};
    my $frac_identical = $hit->{frac_identical};
    my $percent_id = $hit->{percent_id};
    my $score = $hit->{score};
    
	if ($chosen_hit) {
		if ($hit->{frac_identical} < $chosen_hit->{frac_identical}) {
			$logic_name = "following hit with less frac_identical";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			#the previous hit is better
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $trapmap_id);
			$self->load->load_exonerate_trapblock($chosen_hit, $prev_id, $genomic, $qlen);
			return 1;
		} #discard the hit with less score -- last if wu-blast, next if blast
		
		elsif ($hit->{score} > $chosen_hit->{score}) {
			$logic_name = "following hit with less score";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $trapmap_id);
			#the previuos hit is better
			
			$self->load->load_exonerate_trapblock($chosen_hit, $prev_id, $genomic, $qlen);
			return 1;
		} 
		elsif ($hit->{score} < $chosen_hit->{score}) {
			$logic_name = "following hit with more score";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			$chosen_hit = $hit;
			$self->update_trapmap(0,$multiple,$chosen_filter,$checked, $prev_id);
			$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
			#this hit is better
			$self->load->load_exonerate_trapblock($hit, $trapmap_id, $genomic, $qlen);
			return 1;
		} 
		elsif ($hit->{score} == $chosen_hit->{score}) {
			$logic_name = "following hit with equal score";
			$chosen_filter = $self->load->load_logic_name($logic_name);
			$multiple = 1;
			$chosen = 1;
			my $logic_name = "multiple mapping";
			my $chosen_filter = $self->load->load_logic_name($logic_name);
			$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked, $trapmap_id);
			#both are good
			$self->load->load_exonerate_trapblock($chosen_hit, $prev_id, $genomic, $qlen);
			$self->load->load_exonerate_trapblock($hit, $trapmap_id, $genomic, $qlen);
			return 0;
		}
		else{
			#$debug && print STDERR "SOMETHING ELSE...\n";
			die;
		}
	} 
	else {
		$chosen_hit = $hit;
		$chosen=1;
		$logic_name = "first choice";
		$chosen_filter = Bio::TrapCluster::LoadTrap->load_logic_name($logic_name);
		$self->update_trapmap($chosen,$multiple,$chosen_filter,$checked,$trapmap_id);
		$prev_id = $trapmap_id;
		return 0;
	}
	return 0;
}

sub _get_trap{
	my ($self, $value) = @_;
	
	my ($id,$trapname);
        my $current = `date`;
	$debug && print STDERR "$current -- started get_trap\n";
	# ti|1|gi|291202702|gb|AB550827|0|Ayu21-W339||mRNA|||
	#ti|trap_id|gi|gb_id|gb|gb_locus|gss_id|trap_name|clone_id|mol_type|strain|cell_line|vector_name
	my @f = split /\|/, $value;
	if ($f[0] eq 'ti'){$id = $f[1]; $trapname = $f[7]}
	else{$trapname = $value}
	
	unless($trapname){die "no trapname $value\n";}
	$debug && print STDERR "TRAPNAME $trapname\n";
	my $trap = $self->load->fetch->get_trap_by_name($trapname);
    my $trap_id = $trap->{'trap_id'};
	unless($trap_id){
		$self->retrieve->fromquery(1);
		my $traps = $self->retrieve->get_from_query($trapname);
		return 0 unless $traps;
		$self->load->load_db($traps);
		$current = `date`;
		$debug && print STDERR "$current -- loaded trap $trapname\n";
		$trap = $self->load->fetch->get_trap_by_name($trapname);
		$trap_id = $trap->{'trap_id'};
	}
	unless ($trap_id){return 0;}
	$current = `date`;
	$debug && print STDERR "$current -- fetched trap $trap_id\n" if $trap_id;
	
	my $query = qq{SELECT trapmap_id FROM trapmap WHERE trap_id = $trap_id};
	my $trapmap_id = $self->load->fetch->get_trapmap_id_by_query($query);
    
    $current = `date`;
	if ($trapmap_id){
		$debug && print STDERR "$current -- trapcheck returns ".$trapmap_id."\n";
		$self->load->update_trapcheck(1,1, $trap_id);
		
	}
	return $trap;
}



1;
