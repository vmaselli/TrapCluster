#!/usr/bin/perl -w

=head1 $self->load

=head2 Authors

=head3 Created by

Vincenza Maselli v.maselli@ucl.ac.uk

=head2 Description

This module read information from file/web/local database and load the table
`project`, `product`, `genbank_files`, `trap`

=head2 Usage
	
	use $self->load;
	my $traps = $retrieve->get_from_file();
	$self->load->load_db($traps);


=cut

package Bio::TrapCluster::LoadTrap;
require "$ENV{'TrapCluster'}/trapcluster_conf.pl";

use strict; 
use DBI; 
use Carp;
use Data::Dumper; 
use vars qw(@ISA);

use Bio::TrapCluster::Utils::Argument qw(rearrange); 
use Bio::TrapCluster::Utils::Exception qw(throw warning deprecate);
use Bio::Search::HSP::HSPFactory;
@ISA = qw(Bio::Root::Root);

#setting global variables

my %conf =  %::conf;
#my $debug = $conf{'global'}{'debug'};
my $debug = 1;
my $debugSQL =
$conf{'global'}{'debugSQL'};
my $mysql_path =$conf{'default'}{'mysql_path'};
 
my $tmpdir = $conf{'default'}{'tmp_dir'};
my $base_keyword =$conf{'retrieve'}{'base_keyword'};
my $swarm = $conf{'global'}{'swarm'};
 
use Bio::TrapCluster::Utils::File;
use Bio::TrapCluster::Fetch;
use Bio::SeqIO;

=head2 new

  Arg [..]: Take a set of named argumnts from a config file
  Example: my $retrieve = Bio::TrapCluster::RetrieveTrap->new
  Description:
  Returntype: Bio::TrapCluster::
  Exceptions: source (fromfile, fromdb, fromquery) not defined;
  Caller:
  Status: 

=cut

sub new{
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

 
  my $fetch = Bio::TrapCluster::Fetch->new;  
  $self->fetch($fetch);
  return $self;
}
=head2 load_db

Arg : Take an array of trap hash 
Example: $self->load->load_db($traps) 
Description: Calls method to load trap preliminary information into the database 
Returntype: $self->load; 
Exceptions: Insert failure 
Caller: 
Status: stable

=cut

sub load_db { 
	my ($self, $traps) = @_;
	$debug && print STDERR  "USING DATABASE ".$self->fetch->dbname."\n";
	
 	foreach my $value (@{$traps}){ 			
 		
 		my $failure = "insert failure for ".$value->{'name'};
		$debug && print STDERR "Got name ".$value->{'name'}."\n";
 		if ($conf{'loading'}{'complete'}){
			my $trap = $self->load_project($value);
			unless ($trap->{'project_id'}){my $error_id =$self->load_error($failure,"project");warn("$failure")} 
					
			my $esclone_id = $self->load_esclone($trap);
			unless($esclone_id){my $error_id = $self->load_error($failure,"esclone");warn("$failure")} 
			my $analysis_id = $self->load_analysis($value);
			unless($analysis_id){my $error_id =$self->load_error($failure,"analysis");warn("$failure")}
			
			$value->{'esclone_id'} = $esclone_id;
 		}
 		my $trap_id = $self->load_trap($value);
 		unless($trap_id){my $error_id =$self->load_error($failure,"trap");warn("$failure")}
 		$value->{'trap_id'} = $trap_id;
		$self->load_trapadditional($value);
		if ($value->{'fromdb'}){
			my $db = $value->{'fromdb'}{'db'};
			my $oldtrap_id = $value->{'fromdb'}{'oldtrap_id'};
			$self->load->load_trapmap_from_db($db,$oldtrap_id, $trap_id);
		}
	}
	return $self;
}

=head2 load_analysis

Arg : Take a trap hash
Example: $self->load->load_analysis($trap) 
Description: loading the analysis table 
Returntype: analysis dbID; 
Exceptions: Insert failure 
Caller: 
Status: stable

=cut

sub load_analysis{ 
	my ($self, $trap) = @_;
	
	my $logic_name = $trap->{'analysis'}->{'name'};
 	my $logic_name_id =$self->load_logic_name($logic_name);
	my $db = $trap->{'analysis'}->{'db'};
	my $dbversion = $trap->{'analysis'}->{'db_version'};
	my $dbfile = 'NULL';
	$dbfile = $trap->{'analysis'}->{'db_file'} if $trap->{'analysis'}->{'db_file'};
	my $program = $trap->{'analysis'}->{'program'};
	my $program_version = $trap->{'analysis'}->{'program_version'};
	my $description = $trap->{'analysis'}->{'description'};
	my $insert = qq{INSERT INTO analysis SET logic_name_id = \"$logic_name_id\", db= \"$db\", db_version = \"$dbversion\", db_file = \"$dbfile\", program =\"$program\",program_version = \"$program_version\", description = \"description\"	};

	my $analysis_id = $self->fetch->store($insert);
 	return $analysis_id;
}

=head2 load_error

Arg : 'type of error', 'table name', 'table id'
Example: $self->load->load_error($error_note, $table_name, $table_id) 
Description: loading the trap_error table 
Returntype: error dbID; 
Exceptions:  
Caller: 
Status: stable

=cut

sub load_error{ 
	my ($self, $error_note, $table_name, $table_id) = @_;
 	$table_id = 0 unless defined $table_id;
	my $logic_name = $table_name."_error";
 	my $logic_name_id =$self->load_logic_name($logic_name);
 	my $query = qq{SELECT trap_error_id FROM trap_error WHERE logic_name_id = $logic_name_id AND table_id =$table_id AND table_name = \"$table_name\" };
 	my $error_id = $self->fetch->get_error_id_by_query($query);
 	unless ($error_id){
 		my $error_insert= qq{insert into trap_error set logic_name_id = $logic_name_id, table_id =$table_id, table_name = \"$table_name\", note = \"$error_note\"};
 		$error_id =$self->fetch->store($error_insert);
 	}
 	return $error_id;
}

=head2 load_logic_name

Arg : 'logic name'
Example: $self->load->load_logic_name($logic_name) 
Description: loading the logic_name table 
Returntype: logic_name dbID; 
Exceptions: 
Caller: 
Status: stable

=cut

sub load_logic_name{
	
	my ($self,$logic_name) = @_;

 	my $logic_name_id = $self->fetch->get_logic_name_id_by_name($logic_name);
 	unless ($logic_name_id){
 		my $logic_insert = qq{insert into logic_name set name = \"$logic_name\"};
 		$logic_name_id =$self->fetch->store($logic_insert);
 	}
 	return $logic_name_id;

}

=head2 fetch

 Title    : fetch
 Usage    : $obj->fetch([$newval])
 Function : get/set method for attribute fetch 
 Returns  : value of fetch
 Args     : newval of fetch (optional)

=cut 

sub fetch{
  my ($self, $value) = @_;
  $self->{'fetch'} = $value if defined $value;
  
  return $self->{'fetch'};
}

=head2 load_trap

Arg : trap hash
Example: $self->load->load_trap($trap) 
Description: loading the trap table 
Returntype: trap dbID; 
Exceptions: 
Caller: 
Status: stable

=cut

sub load_trap { 
	my ($self, $trap) = @_;
 	my $trap_name = $trap->{'name'};
 	my $trap_id = $self->fetch->get_trap_id_by_name($trap_name);
 	return $trap_id if $trap_id;
 	my $sequence = $trap->{'sequence'};
	my $gb_id = $trap->{'gb_id'};
	my $gb_locus = $trap->{'gb_locus'};
	my $gss_id = $trap->{'gss_id'};
	my $mol_type = $trap->{'mol_type'};
	my $seq_length = $trap->{'seq_length'};
	my $seq_length_not_N = $trap->{'seq_length_not_N'};
	my $max_frag_length_N_splitted = $trap->{'max_frag_length_N_splitted'};
	my $x_masked_seq = $trap->{'x_masked_seq'};
	$x_masked_seq |= "NA";
	my $nrepeat = $trap->{'nrepeat'};
	my $xrepeat = $trap->{'xrepeat'} ;
	my $x_percent_masked = $trap->{'x_percent_masked'} ;
	my $n_percent_masked = $trap->{'n_percent_masked'};
	my $sequencing_direction = $trap->{'sequencing_direction'};
	my $esclone_id = $trap->{'esclone_id'};	
	my $paired = 0;
	$paired = $trap->{'paired_id'} if $trap->{'paired_id'};
 	
	my $insert = qq{INSERT INTO trap SET  
	trap_name                  = \"$trap_name\",
	sequence                   = \"$sequence\",
	gb_id                      = $gb_id,
	gb_locus                   = \"$gb_locus\",
	gss_id                     = $gss_id,
	mol_type                   = \"$mol_type\",
	seq_length                 = $seq_length,
	seq_length_not_N           = $seq_length_not_N,
	max_frag_length_N_splitted = $max_frag_length_N_splitted,
	x_masked_seq         = \"$x_masked_seq\",
	nrepeat              = $nrepeat,
	xrepeat              = $xrepeat,
	x_percent_masked     = $x_percent_masked,
	n_percent_masked     = $n_percent_masked,
	sequencing_direction = \"$sequencing_direction\",
	esclone_id           = $esclone_id,
	paired_tag_id        = $paired
	};
	
	my $splk = 0;
	if ($trap->{'splk'}){$splk = 1}
	
	
	$trap_id = $self->fetch->store($insert);
	my $check_insert = qq{INSERT INTO trapcheck SET trap_id = $trap_id, splk = $splk};
	$self->fetch->store($check_insert);
	
	$debug && print STDERR "Loaded trap $trap_name\n";
	return $trap_id;
}

=head2 load_esclone

Arg : trap hash
Example: $self->load->load_esclone($trap) 
Description: loading the esclone table 
Returntype: esclone dbID; 
Exceptions: 
Caller: 
Status: stable

=cut

sub load_esclone { 
	my ($self, $trap) = @_;

	#$debug && print STDERR Dumper $trap;
 	my $clone = $trap->{'clone_id'};
	unless ($clone){$clone = $trap->{'name'}}
 	my $esclone_id = $self->fetch->get_esclone_id_by_clone($clone);

 	unless ($esclone_id){ 
		my $design_type =$trap->{'sequencing_direction'};
		if ($trap->{'splk'}){$design_type .="splk"}
		else{$design_type .= "\' RACE";} 
		
		if ($trap->{'sequencing_direction'} eq "I" || $trap->{'project_name'} eq 'tigm'){$design_type = "Inverse PCR"}
		elsif($trap->{'project_name'} eq 'escells'){
			if ($trap->{'sequencing_direction'} eq "F"){$design_type = "L"}
			if ($trap->{'sequencing_direction'} eq "R"){$design_type = "R"}
		}
		my $strain =$trap->{'strain'};
		my $clone_genbank_id = 0;
		$clone_genbank_id=$trap->{'clone_genbank_id'} if $trap->{'clone_genbank_id'};
		
		my $cell_line= $trap->{'cell_line'};
		my $vector_name = $trap->{'vector_name'};
		my $vector_type = $trap->{'vector_type'};
		my $insert = qq{INSERT INTO esclone SET clone_genbank_id = \"$clone_genbank_id\", cell_line = \"$cell_line\", strain = \"$strain\", vector_name = \"$vector_name\", design_type = \"$design_type\", clone = \"$clone\", project_id = $trap->{'project_id'}};
		$esclone_id =$self->fetch->store($insert);
	}

 	return $esclone_id;
}

=head2 load_project

Arg : trap hash
Example: $self->load->load_project($trap) 
Description: loading the project table 
Returntype: project dbID; 
Exceptions: 
Caller: 
Status: stable

=cut

sub load_project { 
	my ($self, $trap) = @_;

	my $project_id =$self->fetch->get_project_id_by_name(lc($trap->{'project_name'}));
	unless ($project_id){ 
		my $pname = lc($trap->{'project_name'});
		my $url =$trap->{'project_url'};
		my $wordkey = $trap->{'project_wordkey'};
		my $private = $trap->{'private'};
		my $insert = qq{insert into project set project_name = \"$pname\", url = \"$url\", wordkey = \"$wordkey\", private = $private};
		#$debug && print STDERR "Insert statement\n". $insert."\n";
		$project_id =$self->fetch->store($insert);
 	}
 	$trap->{'project_id'}=$project_id;
 	return $trap;
}



sub _get_from_genbank{ 
	my ($self, $wordkey) = @_;

	my $trap_dbID;
 	my @traps;
 	my $seqio = $self->genbank_exec($base_keyword.$wordkey, $self->mindate);
 	return $seqio;
}

=head2 load_trapmap_from_db

Arg : Database name ,old trap dbID, new trap dbID
Example: $self->load->load_trapmap_from_db('TrapDB', 109, 3) 
Description: loading the trapmap table using information imported from another database
Returntype: boolean; 
Exceptions: 
Caller: 
Status: not finished

=cut

sub load_trapmap_from_db {
	my ($self, $db,$old_trap_id, $new_trap_id) = @_;
	my $trap = $self->fetch->get_trap_by_id($new_trap_id);
	my $tmquery = $conf{'mapping'}{'tm.query'}." WHERE tm.trap_id = $old_trap_id";			
	my $tmres = $db->select_many_from_table($tmquery);
	use Bio::Search::Hit::GenericHit;
	use Bio::Factory::ObjectFactoryI;
	use Bio::Search::HSP::GenericHSP;
	my $hsps;
	foreach my $trapmap (@{$tmres}){
		my $old_trapmap_id = $trapmap->{'trapmap_id'};
		my $tbquery = $conf{'mapping'}{'tb.query'}." WHERE tb.trapmap_id = $old_trapmap_id AND tb.trap_id = $old_trap_id";
		my $tbres = $db->select_many_from_table($tbquery);
		foreach my $trapblock (@{$tbres}){
			my $hsp = Bio::Search::HSP::GenericHSP->new(
														-percent_identity => $trapblock->{'percent_identity'},   
														-hit_name    => $trapmap->{'hit_id'},
														-hit_start   => $trapmap->{'start'},
														-hit_end     => $trapmap->{'end'},
														-hit_length  => ($trapmap->{'end'} - $trapmap->{'start'} + 1),
														-query_length =>$trap->{'seq_length'}, 
														-query_name  => $trap->{'trap_name'},
														-query_start => 1,
														-query_end   => ($trap->{'seq_length'} + 1),
														-query_seq   => $trap->{'sequence'},
														-query_desc  => "trap",
														-start       => $trapblock->{'start'},
														-end         => $trapblock->{'end'},
														-strand      => $trapblock->{'strand'}
        												);
			
			push(@{$hsps}, $hsp);
		}
		my $hit = Bio::Search::Hit::GenericHit->new(
													-name         => $trapmap->{'hit_id'},
													-score        => $trapmap->{'score'},
													-significance => $trapmap->{'significance'},
													-hsps         => $hsps, 
													-hsp_factory  => "Bio::Factory::ObjectFactoryI",
													-lenght => $trapmap->{'frac_identical'},
													-frac_aligned_query => $trapmap->{'frac_aligned_query'}
													);
		
		my $trapmap_id = $self->load_trapmap($hit, $new_trap_id);
		$self->load_trapblock($hit,$trapmap_id);
	}
	return 1;
}

=head2 load_trapmap

Arg : Bio::Search::Hit::HitI 
Example: $self->load->load_trapmap($hit, 3) 
Description: loading the trapmap table
Returntype: trapmap dbID; 
Exceptions: 
Caller: 
Status: stable

=cut

sub load_trapmap {
	my ($self, $feat, $analysis) = @_;
	my $hit_id = $feat->seq_id;
	my $hit_db = $conf{'default'}{'hit_db'};
	my $start = $feat->start; #start della hit 
  	my $end = $feat->end; #end della hit
  	my $strand = $feat->strand;
  	$debug && "HIT START ".$feat->start." $start HIT END $end\n";
	my ($frac) = $feat->get_tag_values('frac_aligned_query') if $feat->has_tag ('frac_aligned_query');
	$frac ||= "NULL";
	my ($nhsps) = $feat->get_tag_values('num_hsps');
	my ($frac_id) = $feat->get_tag_values('frac_identical') if $feat->has_tag('frac_identical');
	$frac_id ||= "NULL";
	my $score = $feat->score;
	$score ||= 0;
	my ($signif) = $feat->get_tag_values('significance') if $feat->has_tag('significance');
	my $trap_id = $feat->display_name;
	$signif ||= "NULL";
	my $logic_name = "EnsEMBL";
	my $target_type = $self->load_logic_name($logic_name);
	unless(defined $analysis){
		$analysis->{'analysis'}->{'name'} = "local alignment";
		$analysis->{'analysis'}->{'db'} = $conf{'default'}{'db_name'};
		$analysis->{'analysis'}->{'db_version'} = $hit_db;
		$analysis->{'analysis'}->{'db_file'} = $conf{'default'}{'db_file'};
		$analysis->{'analysis'}->{'program'} = $conf{'default'}{'program'};
		$analysis->{'analysis'}->{'program_version'} = $conf{'default'}{'program_version'};
		$analysis->{'analysis'}->{'description'} = "generated by ".ref $self;
	}
	my $analysis_id = $self->load_analysis($analysis);
	
	my $query = qq{SELECT trapmap_id FROM trapmap WHERE `trap_id` = $trap_id AND `hit_id` = \"$hit_id\" AND `hit_db` = \"$hit_db\" AND `start` = $start AND `end` = $end  AND  `target_type` = $target_type};
	
	my $trapmap_id = $self->fetch->get_trapmap_id_by_query($query);
	unless($trapmap_id){
		my $insert = qq{INSERT INTO trapmap SET  `trap_id` = $trap_id, `hit_id` = \"$hit_id\", `hit_db` = \"$hit_db\", `start` = $start, `end` = $end, `strand` = \"$strand\",`frac_aligned_query` = $frac, `num_hsps` = $nhsps, `frac_identical` = $frac_id, `score` = $score, `significance` = $signif,  `target_type` = $target_type, `analysis_id` = $analysis_id};
		$trapmap_id = $self->fetch->store($insert);
		$self->update_trapcheck(1,0, $trap_id);
	}
	return $trapmap_id;
	
}

sub load_maxicluster{
	my ($self, $hash, $analysis) = @_;
	
	my $maxicluster_accession = $hash->{'accession'};
	
	my $query = qq{SELECT maxicluster_id FROM maxicluster WHERE  `accession` = \"$maxicluster_accession\"  };
	my $maxicluster = $self->fetch->select_from_table($query);
	
	my $maxicluster_id = $maxicluster->{'maxicluster_id'} if defined $maxicluster;
	unless($maxicluster_id){
		my $insert = qq{INSERT INTO maxicluster SET  `accession` = \"$maxicluster_accession\"};
		$maxicluster_id = $self->fetch->store($insert);
	}
	
	return $maxicluster_id;
	
}


sub load_maxiclustermap {
	my ($self, $hash, $analysis) = @_;
	my $hit_id = $hash->{'hit_id'};
	my $hit_db = $conf{'default'}{'hit_db'};
	my $start = $hash->{'start'}; #start della hit 
  	my $end = $hash->{'end'}; #end della hit
  	my $strand = $hash->{'strand'};
	my $maxicluster_id = $hash->{'maxicluster_id'};
	my $query = qq{SELECT maxiclustermap_id FROM maxiclustermap WHERE  `hit_id` = \"$hit_id\" AND `hit_db` = \"$hit_db\" AND `start` = $start AND `end` = $end  };
	my $maxiclustermap = $self->fetch->select_from_table($query);
	
	my $maxiclustermap_id = $maxiclustermap->{'maxiclustermap_id'} if defined $maxiclustermap;
	
	unless($maxiclustermap_id){
		my $insert = qq{INSERT INTO maxiclustermap SET  `hit_id` = \"$hit_id\", `hit_db` = \"$hit_db\", `start` = $start, `end` = $end, `strand` = \"$strand\",  `maxicluster_id` = $maxicluster_id};
		$maxiclustermap_id = $self->fetch->store($insert);
	}
	
	return $maxiclustermap_id;
	
}


sub load_maxiclusterblock {
	my ($self, $hash) = @_;
	my $maxiclustermap_id = $hash->{'maxiclustermap_id'};
	my $start = $hash->{'start'}; #start della hit 
  	my $end = $hash->{'end'}; #end della hit
  	my $sequence = $hash->{'sequence'};
	my $strand = $hash->{'strand'};
	my $query = qq{SELECT maxiclusterblock_id FROM maxiclusterblock WHERE `start` = $start AND `end` = $end AND  maxiclustermap_id = $maxiclustermap_id};
	my $maxiclusterblock = $self->fetch->select_from_table($query);
	
	my $maxiclusterblock_id = $maxiclusterblock->{'maxiclusterblock_id'} if defined $maxiclusterblock;
	
	unless($maxiclusterblock_id){
		my $insert = qq{INSERT INTO maxiclusterblock SET  `start` = $start, `end` = $end, `strand` = \"$strand\", `maxiclustermap_id` = $maxiclustermap_id, `sequence` = \"$sequence\"};
		$maxiclusterblock_id = $self->fetch->store($insert);
	}
	
	return $maxiclusterblock_id;
	
}


sub load_trap_maxicluster {
	my ($self, $hash) = @_;
	my $trap_id = $hash->{'trap_id'};
	my $trapmap_id = $hash->{'trapmap_id'};
	my $maxicluster_id = $hash->{'maxicluster_id'};
	my $query = qq{SELECT trap_maxicluster_id FROM trap_maxicluster WHERE trap_id = $trap_id AND trapmap_id = $trapmap_id};
	my $trap_maxicluster = $self->fetch->select_from_table($query);
	
	my $trap_maxicluster_id = $trap_maxicluster->{'trap_maxicluster_id'} if defined $trap_maxicluster;
	
	unless($trap_maxicluster_id){
		my $insert = qq{INSERT INTO trap_maxicluster SET trap_id = $trap_id, maxicluster_id = $maxicluster_id, trapmap_id = $trapmap_id };
		$trap_maxicluster_id = $self->fetch->store($insert);
	}

	return $trap_maxicluster_id;
	
}

sub load_trapclustermap {
	my ($self, $hash) = @_;
	my $hit_id = $hash->{'hit_id'};
	my $hit_db = $conf{'default'}{'hit_db'};
	my $start = $hash->{'start'}; #start della hit 
  	my $end = $hash->{'end'}; #end della hit
  	my $strand = $hash->{'strand'};
  	my $trapcluster_id = $hash->{'trapcluster_id'};

	my $query = qq{SELECT trapclustermap_id FROM trapclustermap WHERE `trapcluster_id` = $trapcluster_id AND `hit_id` = \"$hit_id\" AND `hit_db` = \"$hit_db\" AND `start` = $start AND `end` = $end  AND  `strand` = \"$strand\"};
	
	my $trapclustermap = $self->fetch->select_from_table($query);
	my $trapclustermap_id = $trapclustermap->{'trapclustermap_id'} if defined $trapclustermap;
	
	unless($trapclustermap_id){
		my $insert = qq{INSERT INTO trapclustermap SET  `trapcluster_id` = $trapcluster_id, `hit_id` = \"$hit_id\", `hit_db` = \"$hit_db\", `start` = $start, `end` = $end, `strand` = \"$strand\"};
		$trapclustermap_id = $self->fetch->store($insert);
		
	}
	return $trapclustermap_id;
	
}

sub load_trapclusterblock {
	my ($self, $hash) = @_;
	my $trapclustermap_id = $hash->{'trapclustermap_id'};
	my $start = $hash->{'start'}; #start della hit 
  	my $end = $hash->{'end'}; #end della hit
  	my $sequence = $hash->{'sequence'};
	my $strand = $hash->{'strand'};
	
	my $query = qq{SELECT trapclusterblock_id FROM trapclusterblock WHERE `start` = $start AND `end` = $end AND  trapclustermap_id = $trapclustermap_id};
	
	my $trapclusterblock = $self->fetch->select_from_table($query);
	my $trapclusterblock_id = $trapclusterblock->{'trapclusterblock_id'} if defined $trapclusterblock;
	
	unless($trapclusterblock_id){
		my $insert = qq{INSERT INTO trapclusterblock SET  `start` = $start, `end` = $end, `strand` = \"$strand\", `trapclustermap_id` = $trapclustermap_id, `sequence` = \"$sequence\"};
		$trapclusterblock_id = $self->fetch->store($insert);
	}
	
	return $trapclusterblock_id;
	
}


sub load_trap_trapcluster {
	my ($self, $hash) = @_;
	my $trap_id = $hash->{'trap_id'};
	my $trapmap_id = $hash->{'trapmap_id'};
	my $trapcluster_id = $hash->{'trapcluster_id'};
	my $query = qq{SELECT trap_trapcluster_id FROM trap_trapcluster WHERE trap_id = $trap_id AND trapmap_id = $trapmap_id};
	
	my $trap_trapcluster = $self->fetch->select_from_table($query);
	my $trap_trapcluster_id = $trap_trapcluster->{'trap_trapcluster_id'} if defined $trap_trapcluster;
	
	
	unless($trap_trapcluster_id){
		my $insert = qq{INSERT INTO trap_trapcluster SET trap_id = $trap_id, trapcluster_id = $trapcluster_id, trapmap_id = $trapmap_id };
		$trap_trapcluster_id = $self->fetch->store($insert);
	}

	return $trap_trapcluster_id;	
}



sub load_trapcluster{
	my ($self, $hash) = @_;
	
	my $trapcluster_accession = $hash->{'accession'};
	my $maxicluster_id = $hash->{'maxicluster_id'};
	my $ens = $hash->{'link_to_ensembl'};
	my $ucsc = $hash->{'link_to_ucsc'};
	
	my $query = qq{SELECT trapcluster_id FROM trapcluster WHERE  `accession` = \"$trapcluster_accession\"  };
		
	my $trapcluster = $self->fetch->select_from_table($query);
	print Dumper $trapcluster;
	my $trapcluster_id = $trapcluster->{'trapcluster_id'} if defined $trapcluster;
	unless($trapcluster_id){
		my $insert = qq{INSERT INTO trapcluster SET  `accession` = \"$trapcluster_accession\", maxicluster_id = $maxicluster_id, link_to_ensembl = \"$ens\", link_to_ucsc = \"$ucsc\"};
		$trapcluster_id = $self->fetch->store($insert);
	}
	$debug && print STDOUT "$trapcluster_accession has ID $trapcluster_id\n";
	return $trapcluster_id;
	
}



sub load_exonerate_trapmap {
	my ($self, $hit, $trap_id) = @_;
	
	my $hit_id = $hit->{t_id};
	my $hit_db = $conf{'default'}{'hit_db'};
    my $start = $hit->{t_start};
    my $end = $hit->{t_end};
    my $frac = $hit->{frac_aligned_query};
    my $nhsps = 1;
    my $frac_id = $hit->{frac_identical};
    my $score = $hit->{score};
    my $signif = $score;
    
	my $logic_name = "EnsEMBL";
	my $target_type = $self->load_logic_name($logic_name);
	my $analysis;
	$analysis->{'analysis'}->{'name'} = "alignment";
	$analysis->{'analysis'}->{'db'} = $conf{'default'}{'blastdb'};
	$analysis->{'analysis'}->{'db_version'} = $hit_db;
	$analysis->{'analysis'}->{'db_file'} = $conf{'default'}{'blastdb'};
	$analysis->{'analysis'}->{'program'} = "exonerate";
	$analysis->{'analysis'}->{'program_version'} = "2.2.0";
	$analysis->{'analysis'}->{'description'} = "generated by ".ref $self;
	
	my $analysis_id = $self->load_analysis($analysis);
	
	my $query = qq{SELECT trapmap_id FROM trapmap WHERE `trap_id` = $trap_id AND `hit_id` = \"$hit_id\" AND `hit_db` = \"$hit_db\" AND `start` = $start AND `end` = $end AND `num_hsps` = $nhsps AND `score` = $score AND  `target_type` = $target_type};
	
	my $trapmap_id = $self->fetch->get_trapmap_id_by_query($query);
	unless($trapmap_id){
		my $insert = qq{INSERT INTO trapmap SET  `trap_id` = $trap_id, `hit_id` = \"$hit_id\", `hit_db` = \"$hit_db\", `start` = $start, `end` = $end, `frac_aligned_query` = $frac, `num_hsps` = $nhsps, `frac_identical` = $frac_id, `score` = $score, `significance` = $signif,  `target_type` = $target_type, `analysis_id` = $analysis_id};
		$trapmap_id = $self->fetch->store($insert);
		$self->_update_trapcheck(1,0, $trap_id);
	}
	return $trapmap_id;	
}


sub load_trapclustermap_region{
	my ($self, $hash) = @_;
	my $region_id = $hash->{'region_id'};
	my $trapclustermap_id =  $hash->{'trapclustermap_id'};
	my $annotation_ambigous = $hash->{'annotation_ambiguous'};	
	my $number_trapclusterblocks = $hash->{'number_trapclusterblocks'};
	$number_trapclusterblocks ||=0;
	my $number_annotated_trapclusterblocks = $hash->{'number_annotated_trapclusterblocks'};
	unless ($number_annotated_trapclusterblocks){$number_annotated_trapclusterblocks = 0}
	my $exons = $hash->{'total_number_exons'};
	$exons ||= 0;
	my $overlap = $hash->{'overlap'};
	unless ($overlap){$overlap = 0}
	my $type = $hash->{'type'};
	my $query = qq{SELECT trapclustermap_region_id FROM `trapclustermap_region` WHERE `region_id` = $region_id AND `type` = \"$type\"  AND `trapclustermap_id` = $trapclustermap_id};
	my $trapclustermap_region_id = $self->fetch->get_trapclustermap_region_id_by_query($query);
	unless ($trapclustermap_region_id){
		my $insert = qq{INSERT INTO trapclustermap_region SET `region_id` = $region_id, `type` = \"$type\", `annotation_ambiguous` = $annotation_ambigous, `total_number_exons` = $exons,`number_trapclusterblocks` = $number_trapclusterblocks, `number_annotated_trapclusterblocks` = $number_annotated_trapclusterblocks, `overlap` = $overlap, `trapclustermap_id` = $trapclustermap_id};
		$trapclustermap_region_id = $self->fetch->store($insert);
		#$debug && print STDERR "$insert\n";
	}
	else{$debug && print STDERR "\t\tLoadtrapcluster->load_trapclustermap_region $type ALREADY MAPPED\n";}
	return $trapclustermap_region_id;

}

sub load_trapcluster_annotation{
	my ($self, $hash) = @_;
	my $failure;
	my $trapblock_annotation_id = $self->load_trapclusterblock_annotation($hash);
		$failure = "insert failure for ".$hash->{'display_name'}." - ";
		$failure .= $hash->{'region_id'} if $hash->{'region_id'};
		unless($trapblock_annotation_id){my $error_id =$self->load_error($failure,"trapblock_annotation");warn("$failure")}
		
	return $trapblock_annotation_id;
}

sub load_trapclusterblock_annotation{
	my ($self, $hash) = @_;
	$debug && print STDERR "line 365 ";
	$debug && print STDERR ref $self;
	$debug && print STDERR "->load_trapclusterblock_annotation\n";
	my $trapclusterblock_id = $hash->{'trapclusterblock_id'};
	my $rank = 0;
	$rank = $hash->{'rank'} if defined $hash->{'rank'};
	my $display_name = $hash->{'display_name'};
	my $dS_s = 0;
	$dS_s =  $hash->{'dS_s'} if defined $hash->{'dS_s'};
	my $dE_e = 0;
	$dE_e = $hash->{'dE_e'} if defined $hash->{'dE_e'};
	my $logic_name = $hash->{'logic_name'};
	my $source_logic_name_id = $self->load_logic_name($logic_name);
	my $comment = $hash->{'comment'};
	my $label = $hash->{'label'};
	my $label_logic_name_id = $self->load_logic_name($label);
	my $coverage = $hash->{'coverage'};
	my $region_id = 0;
	$region_id = $hash->{'region_id'} if defined $hash->{'region_id'};
	my $flanking_exon_id = 0;
	$flanking_exon_id = $hash->{'flanking_exon_id'} if defined $hash->{'region_id'};
	
	my $type = $self->load_logic_name($hash->{'type'});
	
	my $query = qq{SELECT trapclusterblock_annotation_id FROM `trapclusterblock_annotation` WHERE `trapclusterblock_id` = $trapclusterblock_id  AND `display_name` = \"$display_name\" };
	my $annotation_id = $self->fetch->get_trapclusterblock_annotation_id_by_query($query);
	unless($annotation_id){
	
		my $insert = qq{ INSERT INTO `trapclusterblock_annotation` SET `source_logic_name_id` = $source_logic_name_id, `trapclusterblock_id` = $trapclusterblock_id ,`rank` = $rank , `display_name` = \"$display_name\" , `dS_s` = $dS_s , `dE_e` = $dE_e , `label_logic_name_id` = $label_logic_name_id , `comment` = \"$comment\" , `coverage` = $coverage , `region_id` = $region_id, `flanking_exon_id` = $flanking_exon_id, `type` = $type};
		$annotation_id = $self->fetch->store($insert);
	}
	return $annotation_id;
}
sub load_exonerate_trapblock{
	my ($self,$hit, $trapmap_id, $genomic,$qlen) = @_;
	
    my $start = $hit->{t_start};
    my $end = $hit->{t_end};
    
    my $strand = $hit->{t_strand};
    my $q_length = $hit->{q_length};
    my $identical_bases = $hit->{identical_bases};
    my $frac_identical = $hit->{frac_identical};
    my $perc_id = $hit->{percent_id};
    my $score = $hit->{score};
	my $frac = $hit->{frac_aligned_query};
	my @blocks;

	my $accepted = 0;
	my $min = $qlen * 0.3;
	
	$accepted = 1 if $perc_id >= 95; 
				
	my $query = qq{SELECT trapblock_id FROM trapblock WHERE `trapmap_id`=$trapmap_id AND `start` = $start AND `end` = $end AND `strand` = $strand};
	my $trapblock_id = $self->fetch->get_trapblock_id_by_query($query);
	unless ($trapblock_id){
		my $insert = qq{INSERT INTO trapblock SET `trapmap_id`=$trapmap_id, `start` = $start,`end` = $end, `strand` = $strand, `perc_ident` = $perc_id, `accepted` = $accepted};
		$trapblock_id = $self->fetch->store($insert);
	}
	push (@blocks,$trapblock_id);	
	return \@blocks;
}



sub update_trapcheck {
	my ($self, $mapped, $checked, $trap_id) = @_;

	my $update = qq{UPDATE trapcheck SET mapped = $mapped, checked = $checked WHERE trap_id = $trap_id };
	my $trapcheck_id = $self->fetch->update($update);
	return $trapcheck_id;
}




=head2 load_region

Arg : region hash 
Example: $self->load_region($hash) 
Description: loading region table
Returntype: region dbID; 
Exceptions: 
Caller: 
Status: stable

=cut

sub load_region{
	my ($self, $hash) = @_;
	
	my $region_name = $hash->{'name'};
	unless (defined $hash->{'seq_id'}){$hash->{'seq_id'} = "Exon.".$hash->{'parent_id'}.".".$hash->{'rank'}}
	my $seq_id = $hash->{'seq_id'};
	my $strand = $hash->{'strand'};
	my $start = $hash->{'start'};
	my $end = $hash->{'end'};
	my $refseq = $hash->{'refseq'};
	unless (defined $refseq){$refseq = "NULL";}
	my $parent_id = $hash->{'parent_id'};
	my $rank = $hash->{'rank'};
	my $description = $hash->{'description'};
	my $external_name = "NULL";
	$external_name = $hash->{'external_name'} if $hash->{'external_name'};
	my $biotype = $hash->{'biotype'};
	unless (defined $biotype){$biotype = "NULL";}
	my $logic_name = $biotype if $biotype ne "NULL";
	my $logic_name_id = 0;
	$logic_name_id = $self->load_logic_name($logic_name) if $logic_name;
	
	my $query;
	
	if ($description eq 'genomic'){
		$query = qq{SELECT region_id FROM region WHERE `region_name` = \"$region_name\"  AND `region_strand` = \"$strand\" AND `region_start` = $start AND`region_end` = $end AND `region_strand` = \"$strand\" };
	}
	else{
		$query = qq{SELECT region_id FROM region WHERE `seq_id` = \"$seq_id\" and parent_id = $parent_id};
	}
# 	$debug && print STDERR ref $self;
# 	$debug && print STDERR "->load_region\n";
# 	$debug && print STDERR "REGION $seq_id\n$query ";
	my $region_id = $self->fetch->get_region_id_by_query($query);
# 	$debug && print STDERR "region_id $region_id " if $region_id;
# 	$debug && print STDERR "\n";
	
	my $insert =qq{INSERT INTO region SET `region_name` = \"$region_name\",`seq_id` = \"$seq_id\", `region_strand` = \"$strand\", `region_start` = $start,`region_end` = $end, `refseq_id` = \"$refseq\",`description` = \"$description\", `parent_id` = \"$parent_id\",`rank` = $rank, external_name = \"$external_name\", biotype = $logic_name_id};
	
	unless($region_id){
		#$debug && print STDERR "-- at this time ". `date` . " no region_id exists for $seq_id\n";
		#$debug && print STDERR "try to $insert\n";
		eval{$region_id = $self->fetch->store($insert)};
		unless ($region_id){
			warn $@;
			$debug && print STDERR "-- at this time ". `date` . "Found a duplicate, try to get region_id by $query\n"; 
			$region_id = $self->fetch->get_region_id_by_query($query)
		}
		#else{print STDERR "GOT THIS ID: $region_id\n";}
	}
	
	return $region_id;
}

sub load_trapblock{
	my ($self,$feat) = @_;
	
	my $start = $feat->start; 
	my $end = $feat->end;
	my $strand = $feat->strand;
	my $trapmap_id = $feat->display_name;
	my ($perc_id) = $feat->get_tag_values('percent_identity') if $feat->has_tag('percent_identity');
	$perc_id ||='NULL';
	my ($accepted) = $feat->get_tag_values('accepted');
						
			
	my $query = qq{SELECT trapblock_id FROM trapblock WHERE `trapmap_id`=$trapmap_id AND `start` = $start AND `end` = $end AND `strand` = $strand};
	my $trapblock_id = $self->fetch->get_trapblock_id_by_query($query);
	unless ($trapblock_id){
		my $insert = qq{INSERT INTO trapblock SET `trapmap_id`=$trapmap_id, `start` = $start,`end` = $end, `strand` = \"$strand\", `perc_ident` = $perc_id, `accepted` = $accepted};

		$trapblock_id = $self->fetch->store($insert);
	}
	return $trapblock_id;
}

=head2 load_trapblock_annotation

Arg : annotation hash 
Example: $self->_load_trapblock_annotation($hash) 
Description: loading the trapblock_annotation table
Returntype: trapblock_annotation dbID; 
Exceptions: 
Caller: 
Status: stable

=cut

sub load_annotation{
	my ($self, $hash) = @_;
	my $failure;
	my $trapblock_annotation_id = $self->load_trapblock_annotation($hash);
		$failure = "insert failure for ".$hash->{'display_name'}." - ";
		$failure .= $hash->{'region_id'} if $hash->{'region_id'};
		unless($trapblock_annotation_id){my $error_id =$self->load_error($failure,"trapblock_annotation");warn("$failure")}
		
	return $trapblock_annotation_id;
}


sub load_trapblock_annotation{
	my ($self, $hash) = @_;
	$debug && print STDERR "line 365 ";
	$debug && print STDERR ref $self;
	$debug && print STDERR "->load_trapblock_annotation\n";
	my $trapblock_id = $hash->{'trapblock_id'};
	my $rank = 0;
	$rank = $hash->{'rank'} if defined $hash->{'rank'};
	my $display_name = $hash->{'display_name'};
	my $dS_s = 0;
	$dS_s =  $hash->{'dS_s'} if defined $hash->{'dS_s'};
	my $dE_e = 0;
	$dE_e = $hash->{'dE_e'} if defined $hash->{'dE_e'};
	my $logic_name = $hash->{'logic_name'};
	my $source_logic_name_id = $self->load_logic_name($logic_name);
	my $comment = $hash->{'comment'};
	my $label = $hash->{'label'};
	my $label_logic_name_id = $self->load_logic_name($label);
	my $coverage = $hash->{'coverage'};
	my $region_id = 0;
	$region_id = $hash->{'region_id'} if defined $hash->{'region_id'};
	my $flanking_exon_id = 0;
	$flanking_exon_id = $hash->{'flanking_exon_id'} if defined $hash->{'region_id'};
	
	my $type = $self->load_logic_name($hash->{'type'});
	
	my $query = qq{SELECT trapblock_annotation_id FROM `trapblock_annotation` WHERE `trapblock_id` = $trapblock_id  AND `display_name` = \"$display_name\" };
	my $annotation_id = $self->fetch->get_trapblock_annotation_id_by_query($query);
	unless($annotation_id){
	
		my $insert = qq{ INSERT INTO `trapblock_annotation` SET `source_logic_name_id` = $source_logic_name_id, `trapblock_id` = $trapblock_id ,`rank` = $rank , `display_name` = \"$display_name\" , `dS_s` = $dS_s , `dE_e` = $dE_e , `label_logic_name_id` = $label_logic_name_id , `comment` = \"$comment\" , `coverage` = $coverage , `region_id` = $region_id, `flanking_exon_id` = $flanking_exon_id, `type` = $type};
		$annotation_id = $self->fetch->store($insert);
	}
	return $annotation_id;
}

=head2 load_trapmap_region

Arg : trapmap_region hash 
Example: $self->load->_load_trapmap_region($hash) 
Description: loading the trapmap_region table
Returntype: trapmap_region dbID; 
Exceptions: 
Caller: 
Status: stable

=cut

sub load_trapmap_region{
	my ($self, $hash) = @_;
	my $region_id = $hash->{'region_id'};
	my $trapmap_id =  $hash->{'trapmap_id'};
	my $annotation_ambigous = $hash->{'annotation_ambiguous'};	
	my $number_trapblocks = $hash->{'number_trapblocks'};
	$number_trapblocks ||=0;
	my $number_annotated_trapblocks = $hash->{'number_annotated_trapblocks'};
	unless ($number_annotated_trapblocks){$number_annotated_trapblocks = 0}
	my $exons = $hash->{'total_number_exons'};
	$exons ||= 0;
	my $overlap = $hash->{'overlap'};
	unless ($overlap){$overlap = 0}
	my $type = $hash->{'type'};
	my $query = qq{SELECT trapmap_region_id FROM `trapmap_region` WHERE `region_id` = $region_id AND `type` = \"$type\"  AND `trapmap_id` = $trapmap_id};
	my $trapmap_region_id = $self->fetch->get_trapmap_region_id_by_query($query);
	unless ($trapmap_region_id){
		my $insert = qq{INSERT INTO trapmap_region SET `region_id` = $region_id, `type` = \"$type\", `annotation_ambiguous` = $annotation_ambigous, `total_number_exons` = $exons,`number_trapblocks` = $number_trapblocks, `number_annotated_trapblocks` = $number_annotated_trapblocks, `overlap` = $overlap, `trapmap_id` = $trapmap_id};
		$trapmap_region_id = $self->fetch->store($insert);
		#$debug && print STDERR "$insert\n";
	}
	else{$debug && print STDERR "\t\tLoadTrap->load_trapmap_region $type ALREADY MAPPED\n";}
	return $trapmap_region_id;

}

=head2 load_insertion

Arg : insertion hash 
Example: $self->load_insertion($hash) 
Description: loading the isertion table
Returntype: insertion dbID; 
Exceptions: 
Caller: 
Status: stable

=cut

sub load_insertion{
	my ($self, $hash) = @_;
	
	my $trapmap_id =  $hash->{'trapmap_id'};
	my $trapblock_id =  $hash->{'trapblock_id'};
	my $trap_id =  $hash->{'trap_id'};
	my $gene_id = $hash->{'gene_id'};
	my $putative_insertion_start = $hash->{'putative_insertion_start'};
	my $putative_insertion_end = $hash->{'putative_insertion_end'};
	my $insertion_case = $hash->{'insertion_case'};
	my $insertion_case_id = $self->load_logic_name($insertion_case);
	
	my $insertion_ambigous = $hash->{'insertion_ambiguous'};
	unless (defined $insertion_ambigous){$insertion_ambigous = 0}
	my $insertion_region_start = $hash->{'region_start'};
	my $insertion_region_end = $hash->{'region_end'};
	
	my $query = qq{SELECT insertion_id FROM `insertion` WHERE `trap_id` = $trap_id AND `trapblock_id` = \"$trapblock_id\" AND `trapmap_id` = $trapmap_id AND `gene_id` = $gene_id };
	my $insertion_id = $self->fetch->get_insertion_id_by_query($query);
	unless ($insertion_id){
		$debug && print STDERR "$query\n";
		my $insert = qq{INSERT INTO insertion SET `trap_id` = $trap_id , `trapblock_id` = $trapblock_id, `putative_insertion_start` = $putative_insertion_start, `putative_insertion_end` = $putative_insertion_end, `insertion_case` = $insertion_case_id, `insertion_ambiguous` = $insertion_ambigous, `trapmap_id` = $trapmap_id,`gene_id` = $gene_id, `region_start` = $insertion_region_start,`region_end` = $insertion_region_end};
		#$debug && print STDERR "$insert\n";
		$insertion_id = $self->fetch->store($insert);
	}
	return $insertion_id;

}

=head2 load_trapadditional

Arg : additional hash 
Example: $self->load_trapadditional($hash) 
Description: loading the trapadditional table
Returntype: trapadditional dbID; 
Exceptions: 
Caller: 
Status: stable

=cut

sub load_trapadditional{
	my ($self, $hash) = @_;
		
	my $processing = $hash->{'processing'};
	unless (defined $processing){$processing = "NA"}
	my $quality = $hash->{'quality'};
	unless (defined $quality){$quality = "NA"}
	my $comment = $hash->{'comment'};
	unless (defined $comment){$comment = "NA"}
	my $verified = $hash->{'verified'};
	unless (defined $verified){$verified = 0}
	my $label = $hash->{'label'};
	unless (defined $label){$label = "NA"}
	my $tigem_exref = $hash->{'tigem_exref'};
	unless (defined $tigem_exref){$tigem_exref = "NA"}
	my $user = $hash->{'user'};
	unless (defined $user){$user = "NA"}
	my $note = $hash->{'note'};
	unless (defined $note){$note = "NA"}
	
	my $trap_id = $hash->{'trap_id'};
	my $sql = qq{SELECT trapadditional_id FROM trapadditional WHERE trap_id = $trap_id and label = \"$label\"};
	my $trapadd = $self->fetch->select_from_table($sql);
	my $trapadditional_id = $trapadd->{'trapadditional_id'};
	unless($trapadditional_id){
		my $insert = qq{INSERT INTO `trapadditional` SET `processing` = \"$processing\" ,`quality` = \"$quality\" , `comment` = \"$comment\" , `verified` = \"$verified\" , `label` = \"$label\" , `tigem_ex_ref` = \"$tigem_exref\" , `trap_id` = $trap_id , `user` = \"$user\" ,`note` = \"note\"};
		$trapadditional_id = $self->fetch->store($insert);
	}
	return $trapadditional_id;
}

sub load_trapclusteradditional{
	my ($self, $hash) = @_;
		
	my $processing = $hash->{'processing'};
	unless (defined $processing){$processing = "NA"}
	my $quality = $hash->{'quality'};
	unless (defined $quality){$quality = "NA"}
	my $comment = $hash->{'comment'};
	unless (defined $comment){$comment = "NA"}
	my $verified = $hash->{'verified'};
	unless (defined $verified){$verified = 0}
	my $label = $hash->{'label'};
	unless (defined $label){$label = "NA"}
	my $tigem_exref = $hash->{'tigem_exref'};
	unless (defined $tigem_exref){$tigem_exref = "NA"}
	my $user = $hash->{'user'};
	unless (defined $user){$user = "NA"}
	my $note = $hash->{'note'};
	unless (defined $note){$note = "NA"}
	
	my $trapcluster_id = $hash->{'trapcluster_id'};
	my $sql = qq{SELECT trapclusteradditional_id FROM trapclusteradditional WHERE trapcluster_id = $trapcluster_id and label = \"$label\"};
	my $trapclusteradd = $self->fetch->select_from_table($sql);
	my $trapclusteradditional_id = $trapclusteradd->{'trapclusteradditional_id'};
	unless($trapclusteradditional_id){
		my $insert = qq{INSERT INTO `trapclusteradditional` SET `processing` = \"$processing\" ,`quality` = \"$quality\" , `comment` = \"$comment\" , `verified` = \"$verified\" , `label` = \"$label\" , `tigem_ex_ref` = \"$tigem_exref\" , `trapcluster_id` = $trapcluster_id , `user` = \"$user\" ,`note` = \"note\"};
		$trapclusteradditional_id = $self->fetch->store($insert);
	}
	return $trapclusteradditional_id;
}


=head2 load_unitrap

Arg : unitrap hash 
Example: $self->load_unitrap($hash) 
Description: loading the unitrap table
Returntype: unitrap dbID; 
Exceptions: 
Caller: 
Status: stable

=cut

sub load_unitrap{
	my ($self,$toinsert)=@_;

	my $accession = $toinsert->{'accession'};
	my $chr = $toinsert->{'chr'};
	my $start = $toinsert->{'start'};
	my $end = $toinsert->{'end'};
	my $hit_db = $toinsert->{'hit_db'};	
	my $gene_id = $toinsert->{'gene_id'};
	my $ambiguous = 0;
	
	my $query = qq{SELECT unitrap_id FROM unitrap WHERE chr = \"$chr\" AND hit_db =\"$hit_db\" AND region_id = $gene_id AND start = $start AND end = $end};
	
	my $unitrap_id = $self->fetch->get_unitrap_id_by_query($query);
	
	unless ($unitrap_id){
		$debug && print STDERR "$query\n";
		my $insert = qq{INSERT INTO unitrap SET accession = \"$accession\", chr = \"$chr\", start = $start, end = $end, hit_db = \"$hit_db\", region_id = $gene_id};
		$unitrap_id = $self->fetch->store($insert);
	}
	return $unitrap_id;
}

sub load_unitrap_history{
	my ($self,$newID,$oldID) = @_;
	
	my $query = qq{SELECT history_id FROM unitrap_history WHERE new_unitrap_id = \"$newID\" AND old_unitrap_id = \"$oldID\"};
	
	my $history_id = $self->fetch->select_from_table($query);
	
	unless ($history_id){
		my $insert = qq{INSERT INTO unitrap_history SET new_unitrap_id = \"$newID\", old_unitrap_id = \"$oldID\"};
		$history_id = $self->fetch->store($insert);
	}
	
	return $history_id;	
}


1;
