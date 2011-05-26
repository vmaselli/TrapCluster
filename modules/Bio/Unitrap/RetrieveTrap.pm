#!/usr/bin/perl -w

=head1 Bio::Unitrap::RetrieveTrap

=head2 Authors

=head3 Created by

             Vincenza Maselli
             v.maselli@ucl.ac.uk

=head3 Modified from the original script by

              Guglielmo Roma
              guglielmoroma@libero.it

=head2 Description
        	
        	This module allow to retrieve trap information from different source. Detailed descriptions for each method are available. 
             
=head2 Usage
	
	my $retrieve = Bio::Unitrap::RetrieveTrap->new;
	if ($retrieve->fromfile && $retrieve->todb){
		my $traps = $retrieve->get_from_file();
		Bio::Unitrap::LoadTrap->load_db($traps);
	}
	if ($retrieve->fromdb && $retrieve->todb){
		print STDOUT "loading trap from db\n";
		my $stored_traps = $retrieve->get_from_db("UniTrap");
	}
	    
            
=cut

package Bio::Unitrap::RetrieveTrap;

use strict;
use DBI;
use Carp;
use Data::Dumper;
use vars qw(@ISA);
use Bio::DB::GenBank;
use Bio::Unitrap::Utils::Argument qw(rearrange);
use Bio::Unitrap::Utils::Exception qw(throw warning deprecate);
@ISA = qw(Bio::Root::Root);

#setting global variables

require "$ENV{'Unitrap'}/unitrap_conf.pl";

my %conf =  %::conf;
my $debug = $conf{'global'}{'debug'};
my $debugSQL = $conf{'global'}{'debugSQL'};
my $mysql_path =$conf{'default'}{'mysql_path'};
my $tmpdir = $conf{'default'}{'tmp_dir'};
my $base_keyword = $conf{'retrieve'}{'base_keyword'};
my $swarm = $conf{'global'}{'swarm'};
use Bio::Unitrap::Fetch;
use Bio::SeqIO;
use Bio::Unitrap::Utils::File;


=head2 new

  Arg [..]: Take a set of named argumnts from a config file
  Example: my $retrieve = Bio::Unitrap::RetrieveTrap->new
  Description:
  Returntype: Bio::Unitrap::
  Exceptions: source (fromfile, fromdb, fromquery) not defined;
  Caller:
  Status: 

=cut


sub new{
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my (
    $update,
    $mindate,
    $dir,
    $format,
    $fromfile,
    $fromquery,
    $fromdb,
    $todb,
    $togff,
    $tofasta
    )
    = $self->_rearrange( [
      'UPDATE',
      'MINDATE',
      'DIR',
      'FORMAT',
      'FROMFILE',
      'FROMQUERY',
      'FROMDB',
      'TODB',
      'TOGFF',
      'TOFASTA'
    ],
    @_
    );
	
  unless($update){$update = $conf{'retrieve'}{'update'};}
  unless($mindate){$mindate = $conf{'retrieve'}{'mindate'};}
  unless($dir){$dir = $conf{'retrieve'}{'dir'};}
  unless($format){$format = $conf{'retrieve'}{'format'};}
  unless($fromfile){$fromfile = $conf{'retrieve'}{'fromfile'};}
  unless($fromquery){$fromquery = $conf{'retrieve'}{'fromquery'};}
  unless($fromdb){$fromdb  = $conf{'retrieve'}{'fromdb'};}
  unless($todb){$todb  = $conf{'loading'}{'todb'};}
  unless($togff){$togff  = $conf{'loading'}{'togff'};}
  unless($tofasta){$tofasta  = $conf{'loading'}{'tofasta'};}
  
  if ($fromdb == 0 && $fromfile == 0 && $fromquery == 0){
  	$debug && print STDERR "ERROR: any source defined. Please check the configuration.";
  	exit;
  }
  
  $self->update($update) if defined $update;
  $self->mindate($mindate) if defined $mindate;
  if ($self->update){$dir .= "/Update";}

  $dir && $self->dir($dir);
  
  $format && $self->formatfile($format);
  $self->swarm($swarm)if defined ($swarm);
  $self->fromfile($fromfile) if defined $fromfile;
  $self->fromquery($fromquery)if defined $fromquery;
  $self->fromdb($fromdb) if defined $fromdb;
  
  $self->todb($todb) if defined $todb;
  $self->tofasta($tofasta) if defined $tofasta;
  $self->togff($togff) if defined $togff;
  my $fetch = Bio::Unitrap::Fetch->new;
  $self->fetch($fetch);
  my $utilityfile = Bio::Unitrap::Utils::File->new;
  $self->utilityfile($utilityfile);
  return $self;
}

=head2 update

 Title    : update
 Usage    : $obj->update([$newval])
 Function : get/set method for attribute update 
 Returns  : value of update
 Args     : newval of update (optional)

=cut


sub update{

  my ($self, $value) = @_;
  $self->{'update'} = $value if defined $value;
  return $self->{'update'};
}

=head2 mindate

 Title    : mindate
 Usage    : $obj->mindate([$newval])
 Function : get/set method for attribute mindate 
 Returns  : value of mindate
 Args     : newval of mindate (optional)

=cut 

sub mindate{

  my ($self, $value) = @_;
  $self->{'mindate'} = $value if defined $value;
  
  return $self->{'mindate'};
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

=head2 utilityfile

 Title    : utilityfile
 Usage    : $obj->utilityfile([$newval])
 Function : get/set method for attribute utilityfile 
 Returns  : value of utilityfile
 Args     : newval of utilityfile (optional)

=cut 

sub utilityfile{

  my ($self, $value) = @_;
  $self->{'utilityfile'} = $value if defined $value;
  
  return $self->{'utilityfile'};
}

=head2 format

 Title    : format
 Usage    : $obj->format([$newval])
 Function : get/set method for attribute format 
 Returns  : value of format
 Args     : newval of format (optional)

=cut

sub formatfile{

  my ($self, $value) = @_;
  $self->{'formatfile'} = $value if defined $value;
  
  return $self->{'formatfile'};
}

=head2 swarm

 Title    : swarm
 Usage    : $obj->swarm([$newval])
 Function : get/set method for attribute swarm 
 Returns  : value of swarm
 Args     : newval of swarm (optional)

=cut

sub swarm{

  my ($self, $value) = @_;
  $self->{'swarm'} = $value if defined $value;
  
  return $self->{'swarm'};
}

=head2 fromfile

 Title    : fromfile
 Usage    : $obj->fromfile([$newval])
 Function : get/set method for attribute fromfile 
 Returns  : value of fromfile
 Args     : newval of fromfile (optional)

=cut

sub fromfile{

  my ($self, $value) = @_;
  $self->{'fromfile'} = $value if defined $value;
  
  return $self->{'fromfile'};
}

=head2 fromquery

 Title    : fromquery
 Usage    : $obj->fromquery([$newval])
 Function : get/set method for attribute fromquery 
 Returns  : value of fromquery
 Args     : newval of fromquery (optional)

=cut

sub fromquery{

  my ($self, $value) = @_;
  $self->{'fromquery'} = $value if defined $value;
  
  return $self->{'fromquery'};
}

=head2 fromdb

 Title    : fromdb
 Usage    : $obj->fromdb([$newval])
 Function : get/set method for attribute fromdb 
 Returns  : value of fromdb
 Args     : newval of fromdb (optional)

=cut

sub fromdb{

  my ($self, $value) = @_;
  $self->{'fromdb'} = $value if defined $value;
  
  return $self->{'fromdb'};
}

=head2 todb

 Title    : todb
 Usage    : $obj->todb([$newval])
 Function : get/set method for attribute todb 
 Returns  : value of todb
 Args     : newval of todb (optional)

=cut

sub todb{

  my ($self, $value) = @_;
  $self->{'todb'} = $value if defined $value;
  
  return $self->{'todb'};
}

=head2 togff

 Title    : togff
 Usage    : $obj->togff([$newval])
 Function : get/set method for attribute togff 
 Returns  : value of togff
 Args     : newval of togff (optional)

=cut

sub togff{

  my ($self, $value) = @_;
  $self->{'togff'} = $value if defined $value;
  
  return $self->{'togff'};
}

=head2 tofasta

 Title    : tofasta
 Usage    : $obj->tofasta([$newval])
 Function : get/set method for attribute tofasta 
 Returns  : value of tofasta
 Args     : newval of tofasta (optional)

=cut

sub tofasta{

  my ($self, $value) = @_;
  $self->{'tofasta'} = $value if defined $value;
  
  return $self->{'tofasta'};
}

=head2 dir

 Title    : dir
 Usage    : $obj->dir([$newval])
 Function : get/set method for attribute dir 
 Returns  : value of dir
 Args     : newval of dir (optional)

=cut

sub dir{

  my ($self, $value) = @_;
  $self->{'dir'} = $value if defined $value;
  
  return $self->{'dir'};
}

=head2 new

  Arg [..]: $dbname, $host, $user, $pass, [$port]
  Example: $obj->get_from_db('dbname','localhost','user','3301','passwd');
  Description: This method take in input the parameters to connect to the source database. If some param are common with the new database it will take them form the config file. The query is also specified in the config file. 
  Returntype: trap hash
  Exceptions:
  Caller:
  Status: 

=cut

sub get_from_db{
	my ($self, $dbname, $host, $user, $port, $passwd) = @_;
	my $count = 0;
	$host = $self->fetch->host unless defined $host;
	$user = $self->fetch->user unless defined $user;
	$dbname = $self->fetch->dbname unless defined $dbname;
	$port = $self->fetch->port unless defined $port;
	$passwd = $self->fetch->pass unless defined $passwd;
	my @traps;
	my $db = Bio::Unitrap::self->fetch->new(-dbname=>$dbname);
	my $query = $conf{'retrieve'}{'query'};
	if ($self->swarm){
		my $last_trap = $db->select_from_table($conf{'retrieve'}{'last_id'});
		my $last_id = $last_trap->{'trap_id'};
		if ($debug){$last_id = 1}
		for (my $i = 1; $i <= $last_id; $i++){
			$query .= " AND t.trap_id = $i";
			$debug && print STDERR "$query\n";
			my $res = $db->select_many_from_table($query);
			$query = $conf{'retrieve'}{'query'};
			foreach my $traphref (@{$res}){
				$count ++;
				my $seqdir = $traphref->{'race'};
				$seqdir =~ s/'//;
				my %trap;
				my $trap_name = $traphref->{'trap_name'};
				$debug && print STDERR "loading $trap_name\n";
				#trap table
				$trap{'name'} = $trap_name;
				$trap{'sequence'} =  $traphref->{'sequence'};
				my $gbID = $traphref->{'gb_id'};
				unless (defined $gbID){$gbID = 'NULL'}
				$trap{'gb_id'} = $gbID;
				my $gbLocus = $traphref->{'gb_locus'};
				unless (defined $gbLocus){$gbLocus = "na"}
				$trap{'gb_locus'} = $gbLocus;
				$trap{'mol_type'} = $traphref->{'mol_type'};		
				$trap{'sequencing_direction'} = $seqdir;
				$trap{'seq_length'} = $traphref->{'seq_length'};
				$trap{'seq_length_not_N'} = $traphref->{'seq_length_not_N'};
				$trap{'max_frag_length_N_splitted'} = $traphref->{'max_frag_length_N_splitted'};
				my $xmaskedseq = $traphref->{'x_masked_seq'};
				unless (defined $xmaskedseq){$xmaskedseq = "na"}
				$trap{'x_masked_seq'} = $xmaskedseq;
				$trap{'nrepeat'} = $traphref->{'nrepeat'};
				$trap{'xrepeat'} = $traphref->{'xrepeat'};
				$trap{'x_percent_masked'} = $traphref->{'x_percent_masked'};
				$trap{'n_percent_masked'} = $traphref->{'n_percent_masked'};
				$trap{'gss_id'} = 0;
				#esclone table
				my $clone_id = $traphref->{'clone_id'};
				$debug && print STDERR "IN THIS CASE PROJECT ID =  ";
				$debug && print STDERR $traphref->{'project_name'}."\n";
				unless (defined $clone_id){
					if ($trap_name =~ /Ayu/ || $trap_name =~ /CMHD/ || $traphref->{'project_name'} eq "sanger"){$clone_id = $trap_name}
					elsif($traphref->{'project_name'} eq "ggtc"){
					$clone_id = substr($trap_name,2)}
				}
				$trap{'clone_id'}= $clone_id;
				$trap{'vector_name'} = $traphref->{'vector_name'};
				$trap{'cell_line'} = $traphref->{'cell_line'};		
				$trap{'vector_type'} = $traphref->{'vector_type'};
				$trap{'strain'} = "na";
				#project table
				$trap{'project_name'} = $traphref->{'project_name'};
				$trap{'project_url'} = $traphref->{'url'};
				$trap{'project_wordkey'} = $traphref->{'wordkey'};
				$trap{'private'} = $traphref->{'private'};
				#analysis table
				$trap{'analysis'}{'name'} = "from $dbname database";
				$trap{'analysis'}{'db'} = $dbname;
				$trap{'analysis'}{'db_version'} = "Jul09";
				$trap{'analysis'}{'db_file'} = "na";
				$trap{'analysis'}{'program'} = "na";
				$trap{'analysis'}{'program_version'} = "na";
				$trap{'analysis'}{'description'} = "Retrieve traps from $dbname database by ". ref $self;		
				#trapadditional table
				$trap{'additional'}{'logic_name'} = "imported";
				$trap{'additional'}{'note'} = $traphref->{'note'};
				$trap{'additional'}{'processing'} = $traphref->{'processing'};
				$trap{'additional'}{'user'} = $conf{'dbaccess'}{'user'};
				$trap{'additional'}{'verified'} = $traphref->{'verified'};
				$trap{'additional'}{'quality'} = $traphref->{'quality'};
				$trap{'additional'}{'tigem_exref'} = $traphref->{'old_tigem_trap_name'};
				$trap{'additional'}{'label'} = "imported";
				$trap{'additional'}{'comment'} = "imported from previous version";
				#trapcheck
				$trap{'mapped'} = $traphref->{'mapped'};
				$trap{'checked'} = $traphref->{'checked'};
				$trap{'splk'} = $traphref->{'splk'};
				
				$trap{'fromdb'}{'db'} = $db;
				$trap{'fromdb'}{'oldtrap_id'} = $traphref->{'trap_id'};
				push (@traps, \%trap);				
			}
		}
	}
	return \@traps;
}


sub get_from_gff {
	my ($self, $file) = @_;	
	$debug && print STDERR "FILE $file\n";
	my $version = substr($file,-1);
	$debug && print STDERR "VERSION $version\n";
	# specify input via -fh or -file
	my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => $version) || die $!;
	$debug && print STDERR ref $gffio;
	$debug && print STDERR "\n";
	my $traps;

	# loop over the input stream
	while(my $feature = $gffio->next_feature()) {
		my @tags = $feature->get_all_tags;
		my ($trap_name) = $feature->get_tag_values('Sequence_tag_ID');
		$debug && print STDERR "TRAP  $trap_name\n";
		my $trap = $self->fetch->get_trap_by_name($trap_name);
		$debug && print STDERR Dumper $trap;
		next if $trap;
		$traps = $self->get_from_query($trap_name);
	}
	return $traps;
}

sub get_from_gff_feature {
	my ($self, $feature) = @_;	
	
	my $traps;
	my @tags = $feature->get_all_tags;
	my ($trap_name) = $feature->get_tag_values('Sequence_tag_ID');
	$debug && print STDERR "TRAP  $trap_name\n";
	$traps = $self->get_from_query($trap_name);
	return $traps;
}

sub _get_from_csv{
	my ($self, $file) = @_;
	$debug && print STDERR ref $self;
	$debug && print STDERR "->_get_from_csv\n";
	open (IN, $file) || die $!;
	my @traps;
	while (my $row = <IN>){
		chomp $row;
		my $trap;
		$row =~ s/"//g;
		#`trap_id` , `trap_name` , `sequence` , `gb_id` , `gb_locus` , `gss_id` , `mol_type` , `seq_length` , `seq_length_not_N` , `max_frag_length_N_splitted` , `x_masked_seq` , `nrepeat` , `xrepeat` , `paired_tag_id` , `x_percent_masked` , `n_percent_masked` , `sequencing_direction` , t.esclone_id
		my @f = split /;/, $row;
		#trap table
		$trap->{'name'} = $f[1];
		$trap->{'sequence'} =  $f[2];
		my $gbID = $f[3];
		unless (defined $gbID){$gbID = 'NULL'}
		$trap->{'gb_id'} = $gbID;
		my $gbLocus = $f[4];
		unless (defined $gbLocus){$gbLocus = "na"}
		$trap->{'gb_locus'} = $gbLocus;
		$trap->{'mol_type'} = $f[6];		
		$trap->{'sequencing_direction'} = $f[16];
		$trap->{'seq_length'} = $f[7];
		$trap->{'seq_length_not_N'} = $f[8];
		$trap->{'max_frag_length_N_splitted'} = $f[9];
		my $xmaskedseq = $f[10];
		unless (defined $xmaskedseq){$xmaskedseq = "na"}
		$trap->{'x_masked_seq'} = $xmaskedseq;
		$trap->{'nrepeat'} = $f[11];
		$trap->{'xrepeat'} = $f[12];
		$trap->{'paired_id'} = $f[13]; 
		$trap->{'x_percent_masked'} = $f[14];
		$trap->{'n_percent_masked'} = $f[15];
		$trap->{'gss_id'} = $f[5];
		$trap->{'esclone_id'} = $f[16];
		push (@traps, $trap);
	}
	return \@traps;
}

sub _get_from_fasta{
	my ($self, $file) = @_;
	my $fasta_hash = $self->utilityfile->fromFastaFileToHash($file,$debug);
	my @traps;
	foreach my $fastaheader (keys %{$fasta_hash}){
		my $sequence = $fasta_hash->{$fastaheader}{'sequence'};
		my $seq = Bio::Seq->new(-seq=>$sequence);
		my ($skipgi, $gi,$skipgb,$gb,$gssid,$name,$clone_id,$moltype,$strain,$cell_line,$vector_name) = split /\|/, $fastaheader;
		my $trap;
		#trap table
		$trap->{'name'} = $name;
		$debug && print STDOUT "$name\n";
		$trap->{'sequence'} =  $sequence;
		my $gbID = undef;
		unless (defined $gbID){$gbID = 'NULL'}
		$trap->{'gb_id'} = $gbID;
		my $gbLocus = undef;
		unless (defined $gbLocus){$gbLocus = "na"}
		$trap->{'gb_locus'} = $gbLocus;
		$moltype = "genomic DNA";
		$trap->{'mol_type'} = $moltype;
		$trap = $self->_edit_seq_info($seq,$trap);
		$trap->{'gss_id'} = 0;
		$trap->{'project_name'} = "EUCOMM";
		$trap->{'sequencing_direction'} = substr($trap->{'name'},-4,1);
		$trap->{'splk'} = 1;
		$trap->{'strain'} = $strain;
		$trap->{'clone_id'} = $clone_id;
		$trap->{'cell_line'} = $cell_line;
		$trap->{'vector_name'} = $vector_name;
		$trap->{'vector_type'} = undef;
		$trap->{'analysis'}{'name'} = "from ".$self->formatfile." file";
		$trap->{'analysis'}{'db'} = $self->formatfile;
		$trap->{'analysis'}{'db_version'} = "na";
		$trap->{'analysis'}{'db_file'} = $file;
		$trap->{'analysis'}{'program'} = "na";
		$trap->{'analysis'}{'program_version'} = "na";
		$trap->{'analysis'}{'description'} = "Retrieve traps from ".$self->formatfile." file by ". ref $self;
		push (@traps, $trap);
	}
	return \@traps;
}
=head2 get_from_query

  Arg [..]: Take a set of named argumnts from a config file
  Example: 
  Description:
  Returntype:
  Exceptions:
  Caller:
  Status: 

=cut

sub get_from_query {
	my ($self,$value) = @_;
	$debug && print STDERR ref $self;
	$debug && print STDERR "->get_from_query - _genbank_exec($value)\n";
	my $seqio = $self->_genbank_exec($value);
	$debug && print STDERR "->get_from_query - _parse_genbank_file\n";
	my $traps = $self->_parse_genbank_file($value,$seqio);
	return $traps;

}


=head2 get_from_file

  Arg [..]: Take a set of named argumnts from a config file
  Example: 
  Description:
  Returntype:
  Exceptions:
  Caller:
  Status: 

=cut

sub get_from_file{
	my ($self, $src_file) = @_;
	my $format = $self->formatfile;
	my $traps;
    my $dir = $self->dir;
	if ($src_file){ 
		$debug && print STDERR "\n\n <<<<<<<<< $src_file <$format> >>>>>>>>>>>>\n\n";
		if ($format eq 'GenBank'){
			my $seqio = Bio::SeqIO->new(-format => 'genbank', -file => $src_file);
			$traps = $self->_parse_genbank_file($seqio,$src_file);
		} 
		elsif ($format eq 'csv'){$traps = $self->_get_from_csv($src_file)}
		elsif ($format eq 'gff'){$traps = $self->get_from_gff($src_file)}
		elsif($format eq 'dbGSS'){warning("Sorry method for $format not implemented yet"); return undef;}
		elsif($format eq 'fasta'){$traps = $self->_get_from_fasta($src_file,$debug)}
		else{throw("ERROR: no format was specified. Please check your configuration\n")} 
		return $traps;
	}
	
	my @alltraps;
	
	$debug && print STDERR "WORK on DIR $dir\n"; 
	opendir(DIR,$dir) || die "Error while opening dir $dir: $!";
		
	while (my $file = readdir(DIR)){
	  next if $file =~ /^\./; 
	  next unless $file =~ /gff/; 
	  my $dir = $self->dir; 
	  my $src_file = File::Spec->catfile($dir,$file);
	  
	  if ($format eq 'GenBank'){
		  my $seqio = Bio::SeqIO->new(-format => 'genbank', -file => $src_file);
		  $traps = $self->_parse_genbank_file($seqio,$src_file);
	  } 
	  elsif ($format eq 'csv'){$traps = $self->_get_from_csv($src_file)}
	  elsif ($format eq 'gff'){$traps = $self->get_from_gff($src_file)}
	  elsif($format eq 'dbGSS'){warning("Sorry method for $format not implemented yet"); return undef;}
	  elsif($format eq 'fasta'){$traps = $self->_get_from_fasta($file,$dir)}
	  else{throw("ERROR: no format was specified. Please check your configuration\n")} 
	  next unless $traps;
	  foreach my $t (@{$traps}){push(@alltraps,$t)}	  
	}	
	return \@alltraps;
}

=head2 new

  Arg [..]: Take a set of named argumnts from a config file
  Example: 
  Description:
  Returntype:
  Exceptions:
  Caller:
  Status: 

=cut

sub _parse_genbank_file{
  
  my ($self,$value,$seqio,$file) = @_; 
  unless (defined $seqio){return 0}
  my @traps;
  while (my $seq = $seqio->next_seq) {
    my %trap;
    %trap = %{$self->_get_project($seq,\%trap)};

    %trap = %{$self->_get_name($seq,\%trap)}; # get name 
    next if $trap{'name'} ne $value;
    
    if ($self->update){
		my $sql_check = qq{select * from trap where trap_name = \"$trap{'name'}\"};
		my $res = $self->fetch->select_from_table($sql_check);
		my $trap_dbID = $res->{'trap_id'};
		if ($trap_dbID){
			push (@traps,$res);
			next;
		}
    }
    %trap = %{$self->_edit_seq_info($seq, \%trap)};
    %trap = %{$self->_get_features($seq, \%trap)};
	%trap = %{$self->_get_common_name(\%trap)}; 

    
    %trap = %{$self->_get_seq_direction($seq, \%trap)}; #get direction and splinkerette tag
    $trap{'gb_id'} = $seq->primary_id;
    $trap{'gb_locus'} = $seq->accession_number;
    $trap{'sequence'} =  $seq->seq;
    $trap{'gss_id'} = 0;
    $trap{'user'} = $conf{'dbaccess'}{'user'};
    $trap{'analysis'}{'name'} = "from ".$self->formatfile." file";
    $trap{'analysis'}{'db'} = $self->formatfile;
    $trap{'analysis'}{'db_version'} = "na";
    $trap{'analysis'}{'db_file'} = $file;
    $trap{'analysis'}{'program'} = "na";
    $trap{'analysis'}{'program_version'} = "na";
    $trap{'analysis'}{'description'} = "Retrieve traps from ".$self->formatfile." file by ". ref $self;
    
    push (@traps,\%trap);
  }
  return \@traps;
}


=head2 new

  Arg [..]: Take a set of named argumnts from a config file
  Example: 
  Description:
  Returntype:
  Exceptions:
  Caller:
  Status: 

=cut

sub _get_seq_direction{
  my ($self, $seq, $trap) = @_;
  my %trap = %{$trap};
  my $direction = "na";
  my $splk = 0;
  # $seq is Bio::Seq object; $seq->annotation method returns a Bio::AnnotationCollectionI
  # get_Annotation returns a list of Bio::AnnotationI
  
  foreach my $value ($seq->annotation->get_Annotations("comment")){
    my $comment = $value->as_text;
    my $tag_name = $value->tagname;
    if ($comment =~ /5'/){$direction = "5";} 
    elsif ($comment =~ /3'/){$direction = "3";} 
    if ($comment =~ /Splinkerette/){$splk = 1;}
  }
  if ($trap->{'project_name'} eq 'stanford'){
  	#The sequence tag is generated by 5' or 3' Splinkerette-adaptor PCR (the trap name end with "S").
	#The sequence tag is generated by 3' or 5' RACE (the trap name end with 3 or 5).
	#In some case the trap name end with "-I", in some of these case the sequence tag is generated by inverse PCR and in other by 5' race.
  	if ($trap->{'name'} =~ /-(\d)S$/){$direction = $1}
  	if ($trap->{'name'} =~ /-(\d)$/){$direction = $1}
  	if ($trap->{'name'} =~ /I$/ && $direction eq 'na' ){$direction = "I"}
  }  
  elsif ($trap->{'project_name'} eq 'tigm'){
	  my $tag = substr($trap->{'name'}, -2);
	  if ($tag =~/R/){$direction = "R"}
	  elsif ($tag =~ /F/){$direction = "F"}
  }
  elsif($trap->{'project_name'} eq 'eucomm'){
	  $direction = substr($trap->{'name'},-4,1);
	  $splk = 1;
  }
  elsif($trap->{'project_name'} eq 'escells'){
	  my $tag = substr($trap->{'name'}, -2);
	  if ($tag =~/L/){$direction = "F"}
	  elsif ($tag =~ /R/){$direction = "R"}
  }
  
  $trap{'sequencing_direction'} = $direction;
  $trap{'splk'} = $splk;
  return \%trap;
}

=head2 new

  Arg [..]: Take a set of named argumnts from a config file
  Example: 
  Description:
  Returntype:
  Exceptions:
  Caller:
  Status: 

=cut

sub _get_project{
  my ($self, $seq, $trap) = @_;
  my %trap = %{$trap};
  my $www;
  my ($pname, $url,$wordkey);
  my ($user,$address);
  my ($annotation) = $seq->annotation->get_Annotations("comment");
  my ($reference)= $seq->annotation->get_Annotations("reference");
  my $comment = $annotation->as_text if $annotation;   
  my @anns = split / /, $comment;
  #$debug && print STDERR Dumper @anns;
  for (my $i = 0; $i < scalar @anns; $i++){
    if ($anns[$i] eq 'Contact:'){
    	if ($anns[$i + 1] eq 'Richard'){ #should be tigm
    		$url = $anns[$i + 18];
    		($www,$pname) = split /\./, $url;
    		$wordkey = $pname;}
    	elsif ($anns[$i + 1] eq 'Stanford'){ #should be cmhd
    		$url = $reference->title; 
    		($www,$pname) = split /\./, $url;
    		$wordkey = $pname;
    	}
    	elsif ($anns[$i + 1] eq 'European'){ #should be eucomm
    		$url = $anns[$i + 13];
    		($www,$pname) = split /\./, $url;
    		$wordkey = $pname;
    	}
    	elsif ($anns[$i + 1] eq 'Soriano'){ #should be fhcrh
    		$url = $anns[$i + 65];
    		($www,$pname) = split /\./, $url;
    		$wordkey = $pname;
    	}
    	elsif ($anns[$i + 1] eq 'Zambrowicz'){ #should be lexicon
    		$pname = $anns[$i + 4];
    		$wordkey = "129Sv/Ev";
    	}
    	elsif ($anns[$i + 1] eq 'Ruley'){ #should be vanderbilt
    		$wordkey = lc($anns[$i + 1]);
    		$pname = "vanderbilt";
    	}
    	elsif ($anns[$i + 1] eq 'Exchangeable'){ #should be egtc
    		$pname = "egtc";
    	
    	}
    	elsif ($anns[$i + 1] eq 'Hicks'){ #should be escells
    		$pname = "ESCELLS";
  			$wordkey = "PST";
  			$url = "www.EScells.ca";	
    	}
    	else{ #should be tigem, ggtc, baygenomics 
    	$pname = $anns[$i + 1];}
    }
    elsif($anns[$i] =~ /mail:/){
      my ($user_name,$address) = split /\@/, $anns[$i+1];
      my ($user,$name) = split /\./, $user_name;
      #unless (defined $pname) = split /\./, $address;
      unless (defined $url){$url = "www.".$address;} 
	  unless (defined $pname){$pname = $user;} #works for egtc
	  unless (defined $wordkey){$wordkey = $pname;}
      $i+=2;
    }
    else{
    	next;
    }
  }
  unless (defined $pname){
  	if ($reference->authors =~ /Ishida/){
  		$pname = "NAIStrap";
  		$wordkey = "naist";
  		$url = "http://www2.brc.riken.jp/lab/mouse_es/index.html.en";
  	}
  	elsif($reference->title =~ /EScells/){
  		$pname = "ESCELLS";
  		$wordkey = "PST";
  		$url = "www.EScells.ca";
  	}
  }
  $trap{'comment'} = $comment; 
  
  $trap{'project_name'} = lc($pname);
  $trap{'project_url'} = $url;
  $trap{'project_wordkey'} = $wordkey;
  $trap{'private'} = 0;
  #$debug && print STDERR  Dumper %trap;
  return \%trap;
}

=head2 new

  Arg [..]: Take a set of named argumnts from a config file
  Example: 
  Description:
  Returntype:
  Exceptions:
  Caller:
  Status: 

=cut

sub _get_name{
  my ($self, $seq, ,$trap) = @_;
  my %trap = %{$trap};
  my @desc = split (" ", $seq->desc());
  
  #$debug && print STDERR Dumper @desc;
  my $trap_name = $desc[0];
  my $p_name = $trap{'project_name'};
  if ($trap_name eq 'Mus'){
    if ($p_name eq 'egtc') {
      $trap_name = $desc[8];
      if ($trap_name eq 'genomic'){$trap_name = $desc[7];}
      if($trap_name =~ ":"){
	    my @words = split(":", $trap_name);
	    $trap_name = $words[1];
      }
    }
    elsif($p_name eq 'tigm'){$trap_name = $desc[3]; }
    elsif($p_name eq 'naistrap'){$trap_name = $desc[7];}
  }
  $trap_name =~ s/,// if $trap_name =~ /,/;
  $trap{'name'} = $trap_name;
  #$debug && print STDERR Dumper %trap;
  return \%trap;
}

=head2 new

  Arg [..]: Take a set of named argumnts from a config file
  Example: 
  Description:
  Returntype:
  Exceptions:
  Caller:
  Status: 

=cut

sub _get_common_name{
  my ($self, $trap) = @_;
  my %trap = %{$trap};
  my $trap_name = $trap{'name'};
  my $pname = $trap{'project_name'};
  my $common;
    if ($trap_name =~ /^(3|5)S/){$common = substr($trap_name, 2);}
    elsif ($trap_name =~ /(3|5)SPK$/){($common) = split /\./,$trap_name;}
    elsif ($trap_name =~ /^IST/){$common = substr($trap_name,0,-2)}
    elsif($trap_name =~ /^PST/ && $pname eq 'escells'){$common = substr($trap_name,0,-3)}
    elsif($pname eq 'tigm'){$common = substr($trap_name,-2)}
    else{$common = $trap_name}
  $trap{'common'} = $common;
  return \%trap;
}

sub retrieve_common_name{
  my ($self, $trap) = @_;
  my $trap_name = $trap->{'trap_name'};
  my $pname = $trap->{'project_name'};
  my $common;
    if ($trap_name =~ /^(3|5)S/){$common = substr($trap_name, 2);}
    elsif ($trap_name =~ /(3|5)SPK$/){($common) = split /\./,$trap_name;}
    elsif ($trap_name =~ /^IST/){$common = substr($trap_name,0,-2)}
    elsif($trap_name =~ /^PST/ && $pname eq 'escells'){$common = substr($trap_name,0,-3)}
    elsif($pname eq 'tigm'){$common = substr($trap_name,-2)}
    else{$common = $trap_name}
  return $common;
}

=head2 new

  Arg [..]: Take a set of named argumnts from a config file
  Example: 
  Description:
  Returntype:
  Exceptions:
  Caller:
  Status: 

=cut

sub _edit_seq_info{
 my ($self, $seq, $trap) = @_;
  my %trap = %{$trap};
  my $seq_length = $seq->length;
  my $nrepeat = $seq->seq =~ tr/N/N/;
  my $seq_length_not_N = $seq_length - $nrepeat;
  my $max_frag_length_N_splitted = 0;
  if ($seq_length_not_N < $seq_length) {
    my @fragments = split (/N+/, $seq->seq);
    foreach my $fragment (@fragments) {
      my $frag_length = length($fragment);
      # check for the maximum length of non repated fragments
      if ($frag_length > $max_frag_length_N_splitted){$max_frag_length_N_splitted = $frag_length;}
    }
  } 
  else {$max_frag_length_N_splitted = $seq_length;}
  
  $trap{'max_frag_length_N_splitted'} = $max_frag_length_N_splitted;
  $trap{'seq_length_not_N'} = $seq_length_not_N;
  $trap{'seq_length'} = $seq_length;
  $trap{'nrepeat'} = $nrepeat;
  $trap{'n_percent_masked'} = ($trap{'nrepeat'}/$trap{'seq_length'})*100;
  
  $trap{'x_maskeq_seq'} ="na";
  $trap{'xrepeat'} = 0;
  $trap{'x_percent_masked'} = 0;
  return \%trap;
}

=head2 new

  Arg [..]: Take a set of named argumnts from a config file
  Example: 
  Description:
  Returntype:
  Exceptions:
  Caller:
  Status: 

=cut

sub _get_features{
  my ($self, $seq, $trap) = @_;
  my %trap = %{$trap};
  my @features = $seq->get_all_SeqFeatures;
 
  my ($mol_type, $note, $vector_name,$clone_id, $cell, $strain, $gss_name);
  foreach my $feat (@features) {
    if ($feat->primary_tag eq 'source') {
    	### Extracting the mol_type
		($mol_type) = $feat->get_tag_values('mol_type');
		
      	if ($feat->has_tag('note')) {

			### Extracting the vector_name
			($note) = $feat->get_tag_values('note');
			#'Cell line ID: 21-W550; Gene trap; Vector: pU-21W'
			my $other;
			my @notes = split /;/, $note;
			foreach my $field (@notes){
				if ($field =~ /:/){
					my ($tag, $desc) = split /:/, $field;
					#$debug && print STDOUT "Got $field TAG $tag and DESC $desc\n";
					if ($tag =~ 'Vector'){
						$desc =~ s/ //g;
						$desc =~ s/\n//g;
						$vector_name = $desc;
					}
					elsif($tag eq 'GSS name'){
						$desc =~ s/ //g;
						$desc =~ s/\n//g;
						$gss_name = $desc;
					}
				}
			}
		
		}
		### Extracting the clone_name

		($clone_id) = $feat->get_tag_values('clone') if $feat->has_tag('clone');
		unless ($clone_id) {
				if ($seq->desc =~ /clone\: (.+)\,/) {
				$clone_id = $1;
			}
			else {
				print STDOUT "Could not find clone_id for trap $trap->{'name'}\n";
			}	
		}
		### Extracting the cell line
		($cell) = $feat->get_tag_values('cell_line') if $feat->has_tag('cell_line');
		unless($cell){
			($cell)=$feat->get_tag_values('clone_lib') if $feat->has_tag('clone_lib');
			unless ($cell){
				($cell) = $feat->get_tag_values('cell_type')  if $feat->has_tag('cell_type');
			}
		}
		
		### Extracting the strain
		($strain) = $feat->get_tag_values('strain')  if $feat->has_tag('strain');
      
    }
  }

  $trap{'mol_type'} = $mol_type;
  $trap{'note'} = $note;
  $trap{'clone_id'}= $clone_id;
  $trap{'vector_name'} = $vector_name;
  $trap{'cell_line'} = $cell;
  $trap{'strain'} = $strain;
  $trap{'name'} = $gss_name if defined $gss_name;
  
  return \%trap;
}


=head2 NAME genbank_exec
  
    arguments[2]: keywords (string) and date (date)
    caller: _get_from_genbank
    return: Bio::Seq
    exception: none
    
=cut

sub _genbank_exec () {
	my ($self,$value) = @_;
	if (defined $value){$base_keyword = $value}
	my $dbh = Bio::DB::GenBank->new();
	my $query = Bio::DB::Query::GenBank->new (      -query   =>  $base_keyword,
							-db      => 'nucgss', 
							-mindate => $self->mindate);
	
	my $seqio;
	eval{$seqio = $dbh->get_Stream_by_query($query)};
	
	$self->{'genbank_exec'} = $seqio;
	return $self->{'genbank_exec'};
}

1;