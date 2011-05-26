#!/usr/bin/perl -w

=head1 Bio::Unitrap::LoadTrap

=head2 Authors

=head3 Created by

Vincenza Maselli v.maselli@ucl.ac.uk

=head2 Description

This module read information from file/web/local database and load the table
`project`, `product`, `genbank_files`, `trap`

=head2 Usage

=cut

package Bio::Unitrap::Utils::TextUtility;


use strict; 
use DBI; 
use Carp;
use Data::Dumper; 
use vars qw(@ISA);
use Const::Fast;
use IO::File;

require "$ENV{'Unitrap'}/unitrap_conf.pl";
my %conf =  %::conf;
my $ucsc_db_name = $conf{'default'}{'db_name'};
use Bio::Unitrap::Utils::Argument qw(rearrange); 
use Bio::Unitrap::Utils::Exception qw(throw warning deprecate);
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use Bio::SeqIO;
use File::Spec;
@ISA = qw(Bio::Root::Root);

use Bio::Unitrap::Fetch;
use Bio::Unitrap::Utils::DBSQL::DB;
my $fetch = Bio::Unitrap::Fetch->new;
my $db = Bio::Unitrap::Utils::DBSQL::DB->new;
my $tmpdir = $conf{'default'}{'tmp_dir'};
my $time = localtime;
$time =~ s/ /-/g;
$time =~ s/:/-/g;


sub get_unitrap_gff{
	my ($self) = @_;
	my $file = File::Spec->catfile($tmpdir, $time."unitrap.gff");
	my $gwriter = new Bio::Tools::GFF(  -gff_version => 3,
										-file        => ">>$file");
										
	my @gfeat;									
										
	my $unitrap_query = qq{SELECT * FROM unitrap};
		
	foreach my $unitrap (@{$db->select_many_from_table($unitrap_query)}){
		next if $unitrap->{'chr'} =~ /NT/;
		my $start = $unitrap->{'start'};
		my $end = $unitrap->{'end'};
		if ($start == $end){
			$start = $end - 50;
			$end = $start + 50;
		}
		my $feat = Bio::SeqFeature::Generic->new (   -seq_id   => "chr".$unitrap->{'chr'},
												 -source_tag   => "UniTrap",
												 -primary      => "unitrap", 
												 -start        => $unitrap->{'start'}, 
												 -end          => $unitrap->{'end'},
												 -strand	   => $unitrap->{'strand'},
												 -tag          => {
																	ID => $unitrap->{'accession'}
																  } 
												);
		push(@gfeat, $feat);
		
		my $unitrap_id = $unitrap->{'unitrap_id'};
		my $trap_unitrap_query = qq{SELECT trap_id, trapblock_id FROM trap_unitrap WHERE unitrap_id = $unitrap_id};
		foreach my $trap_unitrap (@{$db->select_many_from_table($trap_unitrap_query)}){
			my $trap_id = $trap_unitrap->{'trap_id'};
			my $trapmap_id = $trap_unitrap->{'trapblock_id'};
			my $trap_query = qq{SELECT t.*, tm.*  FROM trap t, trapmap tm, trapblock tb WHERE tm.trapmap_id = tb.trapmap_id AND tm.trap_id = t.trap_id AND t.trap_id = $trap_id AND tm.trapmap_id = $trapmap_id};
			my $trap = $db->select_from_table($trap_query);
			
			my $feat = Bio::SeqFeature::Generic->new (   -seq_id   => "chr".$unitrap->{'chr'},
												 -source_tag   => "UniTrap",
												 -primary      => "trap hit", 
												 -start        => $trap->{'start'}, 
												 -end          => $trap->{'end'},
												 -score        => $trap->{'score'},
												 -strand     => $unitrap->{'strand'},
												 -tag          => {
																	ID => $unitrap->{'accession'},
																	NAME => $trap->{'trap_name'}
																  } 
												);
			push(@gfeat, $feat);
			
			# my $trapblock_query = qq{SELECT * FROM trapblock WHERE trapmap_id = $trapmap_id};
# 			foreach my $trapblock (@{$db->select_many_from_table($unitrap_query)}){
# 				
# 				my $feat = Bio::SeqFeature::Generic->new (   -seq_id   => "chr".$unitrap->{'chr'},
# 												 -source_tag   => "UniTrap",
# 												 -primary      => "trap block", 
# 												 -start        => $trapblock->{'start'}, 
# 												 -end          => $trapblock->{'end'},
# 												 -score        => $trapblock->{'perc_ident'},
# 												 -strand     => $unitrap->{'strand'},
# 												 -tag          => {
# 																	ID => $trap->{'trap_name'},
# 																	NAME => $trapblock 
# 																  } 
# 												);
# 			push(@gfeat, $feat);
# 			
# 			
# 			}
			
		}	
	}	
	$gwriter->write_feature(@gfeat);
	return $file;
}

sub get_all_trap_fasta{
	my ($self) = @_;
	
	my $file = $time."trapseq.fa";
	my $sql = qq{select t.*, e.* from trap t, esclone e where t.esclone_id=e.esclone_id};
	foreach my $trap (@{$db->select_many_from_table($sql)}){
		my $gb = $trap->{'gb_id'};
	  	unless ($gb){$gb = "na"}
	  	
	  	my $gbacc = $trap->{'gb_locus'};
	  	
	  	my $gssid = $trap->{'gss_id'};
	  	unless ($gssid){$gssid = "na"}
	  	
	  	my $clone_id = $trap->{'clone_id'};
	  	unless ($clone_id) {$clone_id = "na"}
	  	
	  	my $strain = $trap->{'strain'};
	  	unless ($strain){$strain = "na"}
	  	$strain =~ s/ //g;
	  	
	  	my $cell_line = $trap->{'cell_line'};
	  	unless ($cell_line){$cell_line = "na"}
	  	$cell_line =~ s/ //g;
	  	
	  	my $vector_name = $trap->{'vector_name'};
	  	unless ($vector_name){$vector_name = "na"}
	  	$vector_name =~ s/ //g;
	  	
	  	my $mol_type = $trap->{'mol_type'};
	  	$mol_type =~ s/ /_/g;
	  	
	  	my $trapname = "ti|".$trap->{'trap_id'}."|gi|".$trap->{'gb_id'}."|gb|".$trap->{'gb_locus'}."|".$trap->{'gss_id'}."|".$trap->{'trap_name'}."|".$trap->{'clone_id'}."|".$trap->{'mol_type'}."|".$trap->{'strain'}."|".$trap->{'cell_line'}."|".$trap->{'vector_name'};
	  	my $sequence = $trap->{'sequence'};
	  	Bio::Unitrap::Utils::File->writeFasta($file,$tmpdir,1,$trapname,$sequence);
	}
}

sub get_all_trap_by_project_fasta{
	my ($self, $project_id) = @_;
	
	my $file = $project_id."trapseq.fa";
	
	my $sql = qq{select t.*, e.* from trap t, esclone e where t.esclone_id=e.esclone_id and e.project_id = $project_id};
	#print "$sql\n";
	#print $db->dbname,"\n";
	foreach my $trap (@{$db->select_many_from_table($sql)}){
		
	  	my $gb = $trap->{'gb_id'};
	  	unless ($gb){$gb = "na"}
	  	
	  	my $gbacc = $trap->{'gb_locus'};
	  	
	  	my $gssid = $trap->{'gss_id'};
	  	unless ($gssid){$gssid = "na"}
	  	
	  	my $clone_id = $trap->{'clone_id'};
	  	unless ($clone_id) {$clone_id = "na"}
	  	
	  	my $strain = $trap->{'strain'};
	  	unless ($strain){$strain = "na"}
	  	$strain =~ s/ //g;
	  	
	  	my $cell_line = $trap->{'cell_line'};
	  	unless ($cell_line){$cell_line = "na"}
	  	$cell_line =~ s/ //g;
	  	
	  	my $vector_name = $trap->{'vector_name'};
	  	unless ($vector_name){$vector_name = "na"}
	  	$vector_name =~ s/ //g;
	  	
	  	my $mol_type = $trap->{'mol_type'};
	  	$mol_type =~ s/ /_/g;
	  	
	  	
	  	
	  	my $trapname = "ti|".$trap->{'trap_id'}."|gi|".$gb."|gb|".$gbacc."|".$gssid."|".$trap->{'trap_name'}."|".$clone_id."|".$mol_type."|".$strain."|".$cell_line."|".$vector_name;
	  	
	  	my $sequence = $trap->{'sequence'};
	  	Bio::Unitrap::Utils::File->writeFasta($file,$tmpdir,1,$trapname,$sequence);
	  	
	}
}

sub get_trap_by_id_fasta{
	my ($self, $trap_id) = @_;
	
	my $file_name = $trap_id."trapseq.fa";
	
	my $sql = qq{select t.*, e.* from trap t, esclone e where t.esclone_id=e.esclone_id and t.trap_id = $trap_id};
	#print "$sql\n";
	#print $db->dbname,"\n";
	my $trap = $db->select_from_table($sql);
		
	my $gb = $trap->{'gb_id'};
	unless ($gb){$gb = "na"}
	
	my $gbacc = $trap->{'gb_locus'};
	
	my $gssid = $trap->{'gss_id'};
	unless ($gssid){$gssid = "na"}
	
	my $clone_id = $trap->{'clone_id'};
	unless ($clone_id) {$clone_id = "na"}
	
	my $strain = $trap->{'strain'};
	unless ($strain){$strain = "na"}
	$strain =~ s/ //g;
	
	my $cell_line = $trap->{'cell_line'};
	unless ($cell_line){$cell_line = "na"}
	$cell_line =~ s/ //g;
	
	my $vector_name = $trap->{'vector_name'};
	unless ($vector_name){$vector_name = "na"}
	$vector_name =~ s/ //g;
	
	my $mol_type = $trap->{'mol_type'};
	$mol_type =~ s/ /_/g;
	
	my $trapname = "ti|".$trap->{'trap_id'}."|gi|".$gb."|gb|".$gbacc."|".$gssid."|".$trap->{'trap_name'}."|".$clone_id."|".$mol_type."|".$strain."|".$cell_line."|".$vector_name;
	
	my $sequence = $trap->{'sequence'};
	my $file = Bio::Unitrap::Utils::File->writeFasta($file_name,$tmpdir,1,$trapname,$sequence);
	return $file;  	
}

sub get_trap_by_id_range_fasta{
	my ($self, $id, $idcount) = @_;
	my $last = $id + $idcount;
	my $file_name = $id."_".$last."trapseq.fa";
	my $file;
	for (my $i = $id; $i <= $last; $i ++){
		my $trapmaps = $fetch->get_trapmap_by_trap_id($i);
		my $current = `date`;
		print "-- $current ";
		if (scalar @{$trapmaps}){
			$last ++;
			next;
		}
		else {
			my $current = `date`;
			print "\t writing fasta for $i\n";
		}
		my $sql = qq{select t.*, e.* from trap t, esclone e where t.esclone_id=e.esclone_id and t.trap_id = $i};
		#print "$sql\n";
		#print $db->dbname,"\n";
		my $trap = $db->select_from_table($sql);
			
		my $gb = $trap->{'gb_id'};
		unless ($gb){$gb = "na"}
		
		my $gbacc = $trap->{'gb_locus'};
		
		my $gssid = $trap->{'gss_id'};
		unless ($gssid){$gssid = "na"}
		
		my $clone_id = $trap->{'clone_id'};
		unless ($clone_id) {$clone_id = "na"}
		
		my $strain = $trap->{'strain'};
		unless ($strain){$strain = "na"}
		$strain =~ s/ //g;
		
		my $cell_line = $trap->{'cell_line'};
		unless ($cell_line){$cell_line = "na"}
		$cell_line =~ s/ //g;
		
		my $vector_name = $trap->{'vector_name'};
		unless ($vector_name){$vector_name = "na"}
		$vector_name =~ s/ //g;
		
		my $mol_type = $trap->{'mol_type'};
		$mol_type =~ s/ /_/g;
		
		my $trapname = "ti|".$trap->{'trap_id'}."|gi|".$gb."|gb|".$gbacc."|".$gssid."|".$trap->{'trap_name'}."|".$clone_id."|".$mol_type."|".$strain."|".$cell_line."|".$vector_name;
		
		my $sequence = $trap->{'sequence'};
		$file = Bio::Unitrap::Utils::File->writeFasta($file_name,$tmpdir,1,$trapname,$sequence);
	}
	return $file;  	
}


sub get_unitrap_bed{
	my ($self) = @_;
	my $unitraps = $fetch->get_all_unitrap;
	my $dir = "./";
	my $unitrap_file = File::Spec->catfile($dir,"unitrap.bed");
	my $trap_file = File::Spec->catfile($dir,"trap.bed");
	my $unitrapsorted_file = File::Spec->catfile($dir,"unitrapsorted.bed");
	my $trapsorted_file = File::Spec->catfile($dir,"trapsorted.bed");
	my $unitrap_bbfile = File::Spec->catfile($dir,"unitrap.bb");
	my $trap_bbfile = File::Spec->catfile($dir,"trap.bb");
	my $chrsizes = File::Spec->catfile("$ENV{'Unitrap'}/data/",$ucsc_db_name.".chrom.sizes");
	open (UNI, ">$unitrap_file");
	open (TRP, ">$trap_file");
	
	system "$ENV{'Unitrap'}/scripts/utils/fetchChromSizes $ucsc_db_name > $chrsizes";
	open (CHR, $chrsizes);
	my %chrh;
	while (my $r = <CHR>){
		chomp $r;
		my ($chr, $size) = split /\t/, $r;
		$chrh{$chr} = $size;
	}
	foreach my $unitrap (@{$unitraps}){
	
		my $chr = $unitrap->{'chr'};# - The name of the chromosome (e.g. chr3, chrY)
		next if $chr =~ /NT/;
		next if $chr =~/MT/;
		my $start = $unitrap->{'start'};# - The starting position of the feature in the chromosome. The first base in a chromosome is numbered 0.
		my $end = $unitrap->{'end'}; # - The ending position of the feature in the chromosome. 
		#if ($start == $end){$end = $start + 1}
		my $name = $unitrap->{'accession'};# - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
		my $score = 0; # - A score between 0 and 1000.
		my $gene = $fetch->get_region_by_id($unitrap->{'region_id'});
		 
		my $strand = $gene->{'region_strand'};# - Defines the strand - either '+' or '-'.
 		if ($strand eq '-1'){$strand = '-'}
 		else{$strand = '+';}
		print UNI "$chr\t$start\t$end\t$name\t$score\t$strand\n";

  		my $trap_unitraps = $fetch->get_traps_by_unitrap_id($unitrap->{'unitrap_id'});
  		my @starts;

  		my @ends;
  		my @sizes;
  		my %seen;
  		foreach my $trap_unitrap (@{$trap_unitraps}){
  			my $tstrand = $trap_unitrap->{'strand'};# - Defines the strand - either '+' or '-'.
			if ($tstrand eq '-1'){$tstrand = '-'}
			else{$tstrand = '+';}
  			my $str = "$chr\t".$trap_unitrap->{'start'}."\t".$trap_unitrap->{'end'}."\t".$name."-".$trap_unitrap->{'trap_name'}."($tstrand)\t$score\t$tstrand";
  			
  			next if $seen{$str};
  			print STDERR "$str\n";
  			$seen{$str} ++;
  			if ($trap_unitrap->{'start'} < 0){
  			print STDERR " Case 1 $str	\n";
  			}elsif($trap_unitrap->{'end'} > $chrh{"$chr"}){
  				print STDERR "case 2 $str\n";
  			}elsif($trap_unitrap->{'end'} - $trap_unitrap->{'start'} > 1000000 ){
  				print STDERR "case 3 $str\n";
  			}
  			elsif($trap_unitrap->{'end'} == $trap_unitrap->{'start'} ){
  				print STDERR "case 4 $str\n";
  			}
  			else{
  			print TRP "$str\n";
  			}
  		}
  
  	}
  	print "sorting files\n";
  	system "sort -k1,1 -k2,2n $trap_file > $trapsorted_file";
  	system "sort -k1,1 -k2,2n $unitrap_file > $unitrapsorted_file";
	print "converting files \n";
  	system "$ENV{'Unitrap'}/scripts/utils/bedToBigBed $trapsorted_file $chrsizes $trap_bbfile";
  	system "$ENV{'Unitrap'}/scripts/utils/bedToBigBed $unitrapsorted_file $chrsizes $unitrap_bbfile";
  	
  	print "track type=bigBed name=\"TRAP\" description=\"A bigBed trap file\" 	bigDataUrl=http://www2.cancer.ucl.ac.uk/stupkalab/idcc-unitrap/$trap_bbfile\n";
  	print "track type=bigBed name=\"UNITRAP\" description=\"A bigBed unitrap file\" 	bigDataUrl=http://www2.cancer.ucl.ac.uk/stupkalab/idcc-unitrap/$unitrap_bbfile\n";
}

sub splitme {
	my ($self, $in_file) = @_;
	const my $NEW_RECORD_PATTERN => qr/^\QBLASTN 2.2.21 [Jun-14-2009]\E$/;

	my $ifh = IO::File->new( $in_file, O_RDONLY )
		or die "open $in_file: $!";
	
	my $out_file_num = 0;
	
	my $line = $ifh->getline;
	
	while ( $line ) {
		die "Unexpected input: $line" unless $line =~ $NEW_RECORD_PATTERN;
		$out_file_num++;
		my $out_file_name = $in_file . '.' . $out_file_num;
		warn "Writing $out_file_name\n";
		my $ofh = IO::File->new( $out_file_name, O_CREAT|O_EXCL|O_RDWR, 0644 )
			or die "open $out_file_name: $!";
		do {
			$ofh->print( $line );
			$line = $ifh->getline;
		} while $line and $line !~ $NEW_RECORD_PATTERN;
	}
	return 1;
}


1;
