#! /usr/bin/perl -w 

=head2 Authors

=head3 Created by
            
             Vincenza Maselli
             v.maselli@cancer.ucl.ac.uk

=head2 Description
            

=head2 Usage

            

           
=head2 CODE BEGIN

=cut

BEGIN{

  print "Reading settings from $ENV{'Unitrap'}/unitrap_conf.pl\n";
  require "$ENV{'Unitrap'}/unitrap_conf.pl";
  
}
use DBI;
use Getopt::Long;
use Data::Dumper;
use File::Spec;

use Bio::Unitrap::LoadTrap;
use Bio::Unitrap::Fetch;

my $load = Bio::Unitrap::LoadTrap->new;
my $fetch = Bio::Unitrap::Fetch->new;

my %conf =  %::conf;

my ($file,$help);
my $opt = &GetOptions( 
	    "file|f=s" =>    	\$file,
	    "help|h" => \$help
	    );


my $USAGE = "./unitrap_history.pl [-file] [-h help] ";

if ($help){print STDERR "$USAGE
         		  -f file : 
                  -h : this help\n";
			exit;
}

unless ($file){
	print "please enter a file name\n";
	$file = <STDIN>;
}

print "QUI\n";
my $out ="$ENV{'HOME'}/Data/to_test.tab";
open (OUT, ">>$out") || die $!;
open (IN, $file) || die $!;
my $tag = 1;
print "QUO\n";
while (my $row = <IN>){
	chomp $row;
	if ($tag){$tag =0;next;}
	print "QUA\n";
	my ($trap_name,$oldID) = split /;/, $row;
	$trap_name =~ s/"//g;
	$oldID =~ s/"//g;
	
	my $newID = $fetch->get_unitrap_acc_by_trap($trap_name);
	unless ($newID){
		print OUT "$trap_name had an old unitrap id $oldID but not a new one\n";
		next;
	}
	print "$row = <$trap_name>\t<$oldID>\t<$newID>\n";
	$load->load_unitrap_history($newID,$oldID);
}
close(IN);
close(OUT);

END{
print "END\n";
}	