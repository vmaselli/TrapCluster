=head1 Bio::File

=head2 Authors

=head3 Modified by

             Vincenza Maselli
             v.maselli@ucl.ac.uk

=head3 Created by

              Guglielmo Roma
              guglielmoroma@libero.it

=head2 Description
        
             This module have method to handle with files
             
=head2 Example

	  my $file = Bio::File->writeDataToFile("my file","my dir", "my data", 1);
	  my $read_content = Bio::Unitrap::File->readFileData($file,"my dir",1);
	  my $file = Bio::File->appendDataToFile($file,"my dir", "my new data",, 1);

	  my %hash = Bio::File->fromFastaFileToHash($file,$dir,1);
	  my %hash = Bio::File->fromFastaFormatToHash($read_content,1); 
            
=cut


package Bio::Unitrap::Utils::File;

use strict;
use vars qw(@ISA);
use File::Spec;
use Bio::SeqIO;
@ISA = qw(Bio::Root::Root);


sub new{
	my $caller = shift;

 	my $class = ref($caller) || $caller;
  	my $self = $class->SUPER::new(@_);

	return $self;
}


sub readFileData () {
    my ($self, $file, $dir,$debug) = @_;	
	my $content;
	my $file_path = File::Spec->catfile($dir,$file);
	open (FILE, ($file_path)) || die "$! $file_path" ;
	while (defined (my $line=<FILE>)) {$content .= $line;}
	close (FILE);
	return $content;
}

sub writeDataToFile () {
    my ($self, $file,$dir, $data, $debug) = @_;
	my $file_path = File::Spec->catfile($dir,$file);
	#print ref $self;
	#print "->writeDataToFile($file)\n";
	#print "open FILE $file_path\n";
	open (FILE, ">$file_path")|| return undef;
	print FILE "$data";
	close (FILE);
    return $file;
}

sub appendDataToFile () {
    my ($self, $file, $dir,$data, $debug) = @_;
	my $file_path = File::Spec->catfile($dir,$file);
	open (FILE, ">>$file_path")|| return undef;
	print FILE "$data";
	close (FILE);
    return $file;	
}

sub fromFastaFileToHash () {
	my ($self, $file,$dir, $debug) = @_;
	my $file_path = File::Spec->catfile($dir,$file);
	my $cnt = $self->readFileData ($file, $dir, $debug);	
	my %hash = $self->fromFastaFormatToHash ($cnt, $debug);
	return \%hash;
}


sub writeFasta {
	my ($self, $file,$dir, $debug,$id, $seq) = @_;
	my $file_path = File::Spec->catfile($dir,$file);
	open (FILE, ">>$file_path")|| return undef;
	print FILE ">$id\n";
	my $strlen = length $seq;
	my $len = 80;
	my $bin = int($strlen/$len);
	my $lastlen = $strlen - ($len*$bin);
	#print "SEQ LEN $strlen LEN $len BIN $bin LASTLEN $lastlen\n";
	my $lastoffset = 0;
	for (my $i = 0; $i < $bin; $i ++){
	
		my $offset = $len * $i;
		my $expansion = $len;
		
		#print "I = $i OFFSET $offset LEN $expansion\n";
		my $substr = substr($seq,$offset,$expansion);
		print FILE "$substr\n";
		
		$lastoffset = $expansion * ($i+1);
	
	}
	if ($lastlen){
		#print "OFFSET $lastoffset LEN $lastlen\n";
		my $substr = substr($seq,$lastoffset,$lastlen);
		print FILE "$substr\n";
	}
	
	close (FILE);
    return $file_path;
}

# Convert a multifasta text in a hash of sequences. the Key is the sequence name and the value is the sequence
sub fromFastaFormatToHash () {
    my ($self, $text, $debug) = @_;	
	my @array = split (/>/, $text);
	my ($sequence, %hash);
	foreach my $seq (@array) {
		if ($seq) {
			my @lines = split ("\n", $seq);
			my $description = shift @lines;
			
			my @words = split (" ", $description);
			my $name = $words[0];
			my $desc = $words[1];
			
			foreach my $line (@lines) {$sequence = $sequence.$line;}
			
			$name =~ s/\,//g;
			
			$hash {$name}{'sequence'} = $sequence;
			$hash {$name}{'desc'} = $desc;
			
			undef $name;
			undef $sequence;
			undef $desc;
		}
	}
	return %hash; 
}

#remove comment and strange things from sql file
sub purge_sql{
	my ($self, $dir) = @_;
	
	opendir(DIR, $dir);
	while (my $file = readdir(DIR)){
		next if $file =~ /^\./;
		my $filein = File::Spec->catfile($dir,$file);
		my $tmp = $filein.".tmp";
		open(IN, $filein);
		`mv $filein $tmp`;
		open(OUT,">".$filein);
		while (<IN>){
			my $row;
			if ($_ =~ /^-/){next;}
			elsif ($_ =~/^\//){next;}
			elsif ($_ =~ /^SET/){next;}
			elsif($_ =~ /ENGINE/){$row=") ENGINE=MyISAM;\n"}
			else{$row = $_}
			print OUT $row;
		}
	}

}

1;