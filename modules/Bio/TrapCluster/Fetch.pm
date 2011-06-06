=head1 Bio::TrapCluster::Fetch

=head2 Author

             Vincenza Maselli
             v.maselli@ucl.ac.uk

=head2 Description
        
             This module have method to fetch the database
             
=head2 Example

            
=cut


package Bio::TrapCluster::Fetch;

use strict;
use DBI;
use Carp;
use Data::Dumper;
use File::Spec;
use vars qw(@ISA);
require "$ENV{'TrapCluster'}/trapcluster_conf.pl";
use Bio::TrapCluster::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::Root::Root; 
use Bio::TrapCluster::Utils::DBSQL::DB;
@ISA = qw(Bio::Root::Root Bio::TrapCluster::Utils::DBSQL::DB);



my %conf =  %::conf;
my $debug = $conf{'global'}{'debug'};
my $debugSQL = $conf{'global'}{'debugSQL'};
my $mysql_path =$conf{'default'}{'mysql_path'};
my $tmpdir = $conf{'default'}{'tmp_dir'};

sub new {
  my ($caller) =shift @_;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

my (
    
    $host,
    $db,
    $user,
    $pass,
    $port,
    $sth
    )
    = $self->_rearrange( [
    'HOST', 
    'DBNAME',
    'USER',
    'PASS',
    'PORT',
    'STH'
    ],
    @_
    );
  
  unless($host){$host = $conf{'dbaccess'}{'host'}};
  unless($db){$db= $conf{'dbaccess'}{'dbname'}};
  unless($user){$user= $conf{'dbaccess'}{'user'}};
  unless($pass){$pass= $conf{'dbaccess'}{'pass'}};
  unless($port){$port= $conf{'dbaccess'}{'port'}};


  $host && $self->host($host);
  $user && $self->user($user);
  $pass && $self->pass($pass);
  $db && $self->dbname($db);
  $port && $self->port($port);

  my $dbh = DBI->connect("DBI:mysql:database=$db;host=$host;port=$port", $user, $pass);
  unless ($dbh){warn("Can't connect $db; I'll try to create\n");$dbh = $self->create_db} ;
  unless ($dbh){die "Can't create $db: ", $DBI::errst}
  $dbh && $self->db_connection($dbh);


 return $self;
}

sub get_project_id_by_name{
  my ($self, $name) = @_;
  
  my $query = qq{SELECT  project_id FROM project WHERE project_name = \"$name\"};
  my $project = $self->select_from_table($query);
    
  return $project->{'project_id'};
}

sub get_project_by_id{
  my ($self, $id) = @_;
  
  my $query = qq{SELECT * FROM project WHERE project_id = $id};
  my $project = $self->select_from_table($query);
    
  return $project;
}

sub get_logic_name_id_by_name{
  my ($self, $name) = @_;
  
  my $query = qq{SELECT logic_name_id FROM logic_name WHERE name = \"$name\"};
  my $logic_name = $self->select_from_table($query);
    
  return $logic_name->{'logic_name_id'};
}

sub get_esclone_id_by_clone{
  my ($self, $name) = @_;
  
  my $query = qq{SELECT  esclone_id FROM esclone WHERE clone = \"$name\"};
  print "$query\n";
  my $esclone = $self->select_from_table($query);
    
  return $esclone->{'esclone_id'};

}

sub get_esclone_by_id{
  my ($self, $id) = @_;
  
  my $query = qq{SELECT  * FROM esclone WHERE esclone_id = $id};
  my $esclone = $self->select_from_table($query);
    
  return $esclone;

}

sub get_all_trap{
	my ($self) = @_;
	my $query = qq{select t.*, e.*, p.* from trap t, esclone e, project p where e.esclone_id = t.esclone_id and p.project_id = e.project_id limit 3};
	my $trap = $self->select_many_from_table($query);
  	return $trap;
}

sub get_trap_by_id{
	my ($self, $id) = @_;
	my $query = qq{select t.*, e.*, p.* from trap t, esclone e, project p where e.esclone_id = t.esclone_id and p.project_id = e.project_id and trap_id = $id };
	my $trap = $self->select_from_table($query);
  	return $trap;
}

sub get_trapcluster_by_id{
	my ($self, $id) = @_;
	my $query = qq{select t.* from trapcluster where  t.trapcluster_id = $id };
	my $trapcluster = $self->select_from_table($query);
  	return $trapcluster;
}

sub get_trap_id_by_name{
	my ($self, $name) = @_;
	my $trap = $self->get_trap_by_name($name);
  	return $trap->{'trap_id'};
}

sub test_trap_by_id{
	my ($self, $id) = @_;
	my $query = qq{select trap_id from trap where trap_id = $id };
	my $trap = $self->select_from_table($query);
  	return $trap->{'trap_id'};
}

sub get_trap_by_name{
	my ($self, $name) = @_;
	my $query = qq{select t.*, e.*, p.* from trap t, esclone e, project p where e.esclone_id = t.esclone_id and p.project_id = e.project_id and trap_name = \"$name\" };
	my $trap = $self->select_from_table($query);
  	return $trap;
}

sub get_trap_by_gbID{
	my ($self, $id) = @_;
	my $query = qq{select trap_id from trap where gb_id = $id };
	my $trap = $self->select_from_table($query);
  	return $trap;
}

sub get_trap_by_cloneID{
	my ($self, $name) = @_;
	my $query = qq{select t.trap_id from trap t, esclone e where t.esclone_id = e.esclone_id AND e.clone = \"$name\" };
	my $trap = $self->select_from_table($query);
  	return $trap;
}

sub get_trapmap_by_trap_id{
	my ($self, $id) = @_;
	my $query = qq{select * from trapmap where trap_id = $id};
	my $trapmap = $self->select_many_from_table($query);
  	return $trapmap;
}

sub get_trap_id_by_trapmap_id{
	my ($self, $id) = @_;
	my $query = qq{select trap_id from trapmap where trapmap_id = $id};
	my $trap = $self->select_from_table($query);
  	return $trap->{'trap_id'};
}

sub get_chosen_trapmap_by_trap_id{
	my ($self, $id) = @_;
	my $query = qq{select tm.* from trapmap tm  WHERE  trap_id = $id AND chosen = 1};
	my $trapmap = $self->select_many_from_table($query);
  	return $trapmap;
}

sub get_chosen_trapmap_region_by_trap_id{
	my ($self, $id) = @_;
	my $query = qq{select tm.*, tmr.* from trapmap tm, trapmap_region tmr  WHERE  tm.trapmap_id = tmr.trapmap_id AND trap_id = $id AND chosen = 1};
	my $trapmap = $self->select_many_from_table($query);
  	return $trapmap;
}

sub get_trapadditional_by_trap_id{
	my ($self, $id) = @_;
	my $query = qq{select * from trapadditional where trap_id = $id};
	my $trapadditional = $self->select_from_table($query);
  	return $trapadditional;
}

sub get_trapcheck_by_trap_id{
	my ($self, $id) = @_;
	my $query = qq{select * from trapcheck where trap_id = $id};
	my $trapcheck = $self->select_from_table($query);
  	return $trapcheck;
}

sub get_trap_id_by_chosen_trapmap{
	my ($self) = @_;
	my $query = qq{select distinct(trap_id) from trapmap where  chosen = 1};
	my $trapids = $self->select_many_from_table($query);
  	return $trapids;

}

sub get_mapped_but_not_annotated_trap_id{
	my ($self,$start,$limit) = @_;
	my $query = qq{select distinct(trap_id) from trapcheck where  annotated = 0 and checked = 0 and mapped = 1};
	if ($limit){
		$query = qq{select distinct(trap_id) from trapcheck where  annotated = 0 and checked = 0 and mapped = 1 and trap_id >= $start and trap_id < $limit}
	}
	my $trapids = $self->select_many_from_table($query);
  	return $trapids;

}

sub get_mapped_trap_id{
	my ($self) = @_;
	my $query = qq{select distinct(trap_id) from trapcheck where mapped = 1 and checked = 1};
	
	my $trapids = $self->select_many_from_table($query);
  	return $trapids;

}

sub get_not_mapped_trap_id{
	my ($self) = @_;
	my $query = qq{select distinct(trap_id) from trapcheck where mapped = 0};
	
	my $trapids = $self->select_many_from_table($query);
  	return $trapids;

}

sub get_annotated_trap_id{
	my ($self) = @_;
	my $query = qq{select distinct(trap_id) from trapcheck where annotated = 1};
	
	my $trapids = $self->select_many_from_table($query);
  	return $trapids;

}

sub get_trapblock_by_trapmap_id{
	my ($self, $id) = @_;
	
	my $query = qq{select * from trapblock   where   trapmap_id= $id};
# 	$debug && print STDERR ref $self;
# 	$debug && print STDERR "->get_trapblock_by_trapmap_id($id)\n";
# 	$debug && print STDERR "$query\n";
	my $trapblock = $self->select_many_from_table($query);
  	return $trapblock;
}

sub get_trapblock_and_annotation_by_id{
	my ($self, $id) = @_;
	
	my $query = qq{select tb.*, tba.* from trapblock tb, trapblock_annotation tba  where tb.trapblock_id = tba.trapblock_id  AND tb.trapblock_id= $id};
	my $trapblock = $self->select_from_table($query);
  	return $trapblock;
}

sub get_trapbmap_region_by_trapmap_id{
	my ($self, $id) = @_;
	my $trapmapregion_id;
	my $query = qq{select trapmap_region_id from trapmap_region where trapmap_id= $id};
	my $trapmapregion = $self->select_from_table($query);
	$trapmapregion_id = $trapmapregion->{'trapmap_region_id'} if defined $trapmapregion;
  	return $trapmapregion_id;
}

sub get_selected_trapblock_by_trapmap_id{
	my ($self,$trapmap_id) = @_;
	my @blocks;
	my %block_hash;
	my $max_length = 0;
	my $i = 0;
	my $prev_end = 0;

	
	my @tbs = @{$self->get_trapblock_by_trapmap_id($trapmap_id)};

	my $retry = 1;
	if (scalar @tbs == 0 && $retry){
		$retry = 0;
		@tbs = @{$self->get_trapblock_by_trapmap_id($trapmap_id)};
	}
	if (scalar @tbs == 0 && $retry == 0){return 0}
	if (scalar @tbs == 1){return \@tbs}
	my @sorted = sort  {$a->{'start'} <=> $b->{'start'}} @tbs;
	foreach my $block (@sorted){
		my $length = $block->{'end'} - $block->{'start'} + 1;
		if ($block->{'start'} > $prev_end){$i ++;}
		if ($length >= $max_length){
			$max_length = $length;
			$block_hash{$i} = $block;
		}
		$prev_end = $block->{'end'};
	}
	foreach my $k (sort keys %block_hash){
		push(@blocks,$block_hash{$k});
	}
	return \@blocks;
}

sub get_trapblock_and_annotation_by_trapmap_id{
	my ($self, $id) = @_;
	
	my $query = qq{select tb.*, tba.* from trapblock tb, trapblock_annotation tba where tb.trapblock_id = tba.trapblock_id  AND tb.trapmap_id= $id};
	my $trapblock = $self->select_many_from_table($query);
  	return $trapblock;
}

sub get_trapblock_by_trapblock_annotation_id{
	my ($self, $id) = @_;
	
	my $query = qq{select tb.*, tba.* from trapblock tb, trapblock_annotation tba where tb.trapblock_id = tba.trapblock_id AND tba.trapblock_annotation_id = $id};
	my $trapblock = $self->select_many_from_table($query);
  	return $trapblock;
}

sub get_trapmap_by_trapblock_id{
	my ($self, $id) = @_;
	my $query = qq{select tm.*, tmr.* from trapmap tm, trapblock tb, trapmap_region tmr  WHERE tm.trapmap_id = tb.trapmap_id AND tm.trapmap_id = tmr.trapmap_id AND tb.trapblock_id = $id};
	my $trapmap = $self->select_from_table($query);
  	return $trapmap;
}

sub get_trapmap_by_id{
	my ($self, $id) = @_;
	my $query = qq{select tm.*, tmr.* from trapmap tm, trapmap_region tmr  WHERE  tm.trapmap_id = tmr.trapmap_id AND tm.trapmap_id = $id};
	my $trapmap = $self->select_from_table($query);
  	return $trapmap;
}

sub test_trapmap_by_id{
	my ($self, $id) = @_;
	my $query = qq{select trapmap_id from trapmap  WHERE trapmap_id = $id};
	my $trapmap = $self->select_from_table($query);
  	return $trapmap->{'trapmap_id'};
}

sub test_ambiguous_by_id{
	my ($self, $id) = @_;
	my $query = qq{select tmr1.trapmap_id from trapmap_region tmr1, trapmap_region tmr2 where tmr1.trapmap_id = tmr2.trapmap_id and tmr1.annotation_ambiguous = 1 and tmr2.annotation_ambiguous = 0 and tmr1.type = tmr2.type and tmr1.type != 'genomic' and tmr2.trapmap_id = $id};
	my ($trapmap) = @{$self->select_many_from_table($query)};
  	return $trapmap->{'trapmap_id'};
}

sub get_trap_from_project{
	my ($self, $id) = @_;
	my $query = qq{select t.*, e.*, p.* from trap t, esclone e, project p where e.project_id = p.project_id and t.esclone_id = e.esclone_id and e.project_id = $id order by t.trap_name};
	my $trap = $self->select_many_from_table($query);
  	return $trap;
}

sub get_unitrap_by_id{
	my ($self, $id) = @_;
	my $query = qq{select * from unitrap where unitrap_id = $id};
	my $unitrap = $self->select_from_table($query);
  	return $unitrap;
}

sub get_unitrap_acc_by_trap{
	my ($self, $name) = @_;
	my $query = qq{select u.accession from unitrap u, trap_unitrap tu, trap t where tu.unitrap_id = u.unitrap_id AND tu.trap_id = t.trap_id AND t.trap_name = \"$name\"};
	my $unitrap = $self->select_from_table($query);
  	return $unitrap->{'accession'};
}

sub get_last_unitrap_id{
	my ($self, $id) = @_;
	my $query = qq{select * from unitrap order by unitrap_id desc limit 1};
	my $unitrap = $self->select_from_table($query);
  	return $unitrap->{'unitrap_id'};
}

sub get_region_id_by_query{
	my ($self, $query) =@_;
	my $region = $self->select_from_table($query);
	return $region->{'region_id'};
}

sub get_region_by_rank_parent_id{
	my ($self, $rank,$id) = @_;
	my $query = qq{select * from region where rank = $rank and parent_id = $id};
	my $region = $self->select_from_table($query);
	return $region;
}

sub get_intron_by_rank_parent_id{
	my ($self, $rank,$id) = @_;
	my $query = qq{select * from region where description  like "%intron" and rank = $rank and parent_id = $id};
	my $region = $self->select_from_table($query);
	return $region;
}
sub get_all_exons_by_gene_id{
	my ($self, $id) = @_;
	my $query = qq{select * from region where description like \"%exon%\" AND parent_id = $id};
	my $exons = $self->select_many_from_table($query);
	return $exons;
}

sub get_all_transcript_id{
	my ($self) = @_;
	my $query = qq{select region_id from region where description like \"%transcript\"};
	my $transcripts = $self->select_many_from_table($query);
	return $transcripts;
}


sub get_all_unitrap{
	my ($self) = @_;
	my $query = qq{SELECT * FROM `unitrap` };
	my $unitrap = $self->select_many_from_table($query);
	return $unitrap;
}

sub get_all_insertions{
	my ($self) = @_;
	my $query = qq{select i.*, r.seq_id, tm.hit_id, tm.hit_db from insertion i, region r, trapmap tm, trapmap_region tmr where r.region_id = i.gene_id and tm.trapmap_id = i.trapmap_id and tm.trapmap_id = tmr.trapmap_id and tm.chosen = 1 and tmr.region_id = i.gene_id and tmr.annotation_ambiguous = 0  order by i.gene_id, i.putative_insertion_start};
	my $insert = $self->select_many_from_table($query);
	return $insert;
}

sub get_insert_from_range{
	my ($self, $region_id, $start, $end) = @_;
	my $query = qq{select * from `insertion` where transcript_id = $region_id AND putative_insertion_start = $start AND putative_insertion_end = $end};
	my $insert = $self->select_many_from_table($query);
	return $insert;
}

sub get_insert_from_trapblock_id{
	my ($self, $id) = @_;
	my $query = qq{select * from `insertion` where trapblock_id = $id };
	my $insert = $self->select_many_from_table($query);
	return $insert;
}

sub get_insertion_by_region_id{
	my ($self, $id) = @_;
	my $query = qq{select * from `insertion` where gene_id = $id };
	my $insert = $self->select_many_from_table($query);
	return $insert;
}

sub get_region_by_id{
	my ($self, $id) = @_;
	my $query = qq{select * from region where region_id = $id};
	my $region = $self->select_from_table($query);
	return $region;
}

sub get_region_by_parent_id{
	my ($self, $id) = @_;
	my $query = qq{select * from region where parent_id = $id};
	my $region = $self->select_many_from_table($query);
	return $region;
}

sub get_region_by_seq_id{
	my ($self, $id) = @_;
	my $query = qq{select * from region where seq_id = \"$id\"};
	my $region = $self->select_from_table($query);
	return $region;
}

sub get_intron_by_seq_id{
	my ($self, $id) = @_;
	my $query = qq{select * from region where seq_id like \"$id%\"};
	my $region = $self->select_from_table($query);
	return $region;
}


sub get_region_by_external_name{
	my ($self, $id) = @_;
	my $query = qq{select * from region where external_name = \"$id\"};
	my $region = $self->select_from_table($query);
	return $region;
}

sub get_previuos_exon_by_intron_id{
	my ($self, $id, $parent_id) = @_;
	$id --;
	my $query = qq{select * from region where region_id = $id and parent_id = $parent_id};
	my $region = $self->select_from_table($query);
	return $region;
}

sub get_next_exon_by_intron_id{
	my ($self, $id, $parent_id) = @_;
	$id ++;
	my $query = qq{select * from region where region_id = $id and parent_id = $parent_id};
	my $region = $self->select_from_table($query);
	return $region;
}

sub get_all_region_in_trapblock_annotation{
	my ($self) = @_;
	my $query = qq{select distinct(region_id) from trapblock_annotation where source_logic_name_id != 22};
	my $region = $self->select_many_from_table($query);
	return $region;
}

sub get_all_annotated_trapblock{
	my ($self) = @_;
	my $query = qq{select tba.*, tb.* from trapblock_annotation tba, trapblock tb where tba.trapblock_id = tb.trapblock_id AND tba.label_logic_name_id != 22  ORDER BY tba.region_id ASC};
	my $region = $self->select_many_from_table($query);
	return $region;
}

sub get_traps_by_unitrap_id{
	my ($self,$id) = @_;
	my $query = qq{select distinct(tm.trapmap_id),tm.start, tm.end, tm.strand,u.chr, concat(t.sequencing_direction,"_",t.trap_name) as trap_name 
FROM trap_unitrap tu, trapmap tm , trap t, insertion i, unitrap u 
WHERE t.trap_id = tu.trap_id 
AND tu.unitrap_id = u.unitrap_id 
AND tu.trap_id = i.trap_id 
AND i.trapmap_id = tm.trapmap_id 
AND tm.trap_id = i.trap_id
AND tm.chosen = 1 
AND u.chr = tm.hit_id
AND tu.unitrap_id = $id};
	my $trap_unitrap = $self->get_trap_unitrap_by_query($query);
	return $trap_unitrap;
}

sub get_trap_unitrap_by_query{
	my ($self, $query) = @_;
	my $trap_unitrap = $self->select_many_from_table($query);
	return $trap_unitrap;
}

sub get_trapmap_id_by_query{
	my ($self, $query) = @_;
	my $trapmap = $self->select_from_table($query);
	return $trapmap->{'trapmap_id'}
}

sub get_error_id_by_query{
	my ($self, $query) = @_;
	my $error = $self->select_from_table($query);
	return $error->{'trap_error_id'}
}

sub get_trapblock_id_by_query{
	my ($self, $query) = @_;
	my $annotation = $self->select_from_table($query);
	return $annotation->{'trapblock_id'};
}

sub get_trapblock_annotation_id_by_query{
	my ($self, $query) = @_;
	my $annotation = $self->select_from_table($query);
	return $annotation->{'trapblock_annotation_id'};
}

sub get_trapmap_region_id_by_query{
	my ($self, $query) = @_;
	my $trapmap_region = $self->select_from_table($query);
	return $trapmap_region->{'trapmap_region_id'}
}

sub get_insertion_id_by_query{
	my ($self, $query) = @_;
	my $insertion = $self->select_from_table($query);
	return $insertion->{'insertion_id'}
}

sub get_trapcheck_by_query{
	my ($self, $query) = @_;
	my $check = $self->select_from_table($query);
	return $check;
}


sub get_trapadditional_id_by_query{
	my ($self, $query) = @_;
	my $trapmap_region = $self->select_from_table($query);
	return $trapmap_region->{'trapadditional_id'}
}

sub get_unitrap_id_by_query{
	my ($self, $query) = @_;
	my $unitrap = $self->select_from_table($query);
	return $unitrap->{'unitrap_id'}
}


sub store{
  my ($self, $insert) = @_;
  my $dbID = $self->insert_set($insert);
  return $dbID;

}

sub update{
	my ($self, $update) = @_;
	my $dbID = $self->update_set($update);
  	return $dbID;
}

sub delete{
	my ($self, $delete) = @_;
	$self->delete_set($delete);
  	return 1;
}
