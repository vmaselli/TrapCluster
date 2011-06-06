BEGIN {
    # This section contains the library configuration
    package main;

    %conf = (
		# databases access parameters
		# Details of the Ensembl Database as well as the Unitrap database to be created
		'dbaccess' =>{
						'enshost' => "ensembldb.ensembl.org",
						'ensuser' => "anonymous",
						'enspass' => "",
						'ensport'=>5306,
						'ensdbname' => "mus_musculus_core_61_37n",
						'ensestdbname' => "mus_musculus_otherfeatures_61_37n", 
						'ens_hs_dbname' => "homo_sapiens_core_61_37f", 
						'martdbname' => "ensembl_mart_61",
						'host' => "localhost",
						'user' => "mysql-dev",
						#'dbname' => "UniTrapPipeline",
						'dbname' => 'trapcluster',
						'pass' => "dEvEl0pEr",
						'port' => 3306
					},
        # default variables 
	    'default' =>{
						'hit_db' => "NCBIM37",#Version of the mouse genome to use
						'ensembl_version' =>"61",
						'ensembl_base_url' => "http://www.ensembl.org/Mus_musculus/",
						'db_name' => "mm9", #UCSC mouse db version
						'tmp_dir' => "/tmp/",
						'mysql_path' => '/usr/bin/',
						'version' => 'blat',
						'db_file' => '/mnt/store1/comparative_genomics/suzana/indexes/mm9.2bit',
						'program' => 'Blat',
						'program_version' => 'BLAT gfClient v. 34',
						'tmpdir' =>'/tmp/'
					},
        # user defined variables
	    'global' =>{
						'debug' => 1,
						'debugSQL' => 0,
						'swarm' => 0, # set 1 to run in parallel mode
					},
		'retrieve' =>{
						'update' => 0, #Set to 1 if you are performing an update on an existing db
						'fromfile' => 1, #Set to 1 if starting data is from files (Genbank, FASTA, GSS) 
						'fromquery' => 0, #Set to 1 if starting from an NCBI query
						'mindate' => "2000/01/01", #Date from which NCBI query should start
						'fromdb' => 0, #Set to 1 to start from a Unitrap DB already populated with traps
						'fromcsv' =>0, #Set to 1 to start from CSV files 
						#'dir' => "$ENV{'HOME'}/data/fasta", # file directory
						#'dir' => "$ENV{'HOME'}/Data/Projects/Trap/gb", # test genbank files directory
						#'dir' => "$ENV{'HOME'}/Data/Projects/Trap/dbGSS", # gss files directory
						#'dir' => "$ENV{'HOME'}/Data/Projects/Trap/gff/", # gff files directory	
						'dir'  => "$ENV{'HOME'}/Data/Projects/Trap/mgi/",
						#'format' => 'GenBank',
						#'format' => 'dbGSS',
						#'format' =>'fasta',
						#'format' => 'csv',
						'format' => 'gff',
						'base_keyword' => "PST26490",
						'dbname' =>'otherdbname'
						
					},
		'loading' =>{
						'todb' => 1, #Set to 1 to write results to DB
						'togff' => 0, #Set to 1 to write results to GFF
						'tofasta' => 0,#Write to FASTA
						'complete' => 1 #Perform the complete pipeline or only partial
					},
		#Configuration for BLAST parameters
		'application'=>{
						'genomic_blast' =>
										  {
													  'e' => 10,
													  'G' => 5,
													  'E' => 2,
													  'q' => -3,
													  'r' => 1,
													  'v' => 500,
													  'b' => 250,
													  'g' => 'T',
													  'W' => 11,
													  'T' => 'F',
													  'n' => 'F',
													  'K' => 0,
													  'a' => 4
													 },
						'mrna_blast' => 		{
													  'e' => 10,
													  'G' => 5,
													  'E' => 2,
													  'q' => -3,
													  'r' => 1,
													  'v' => 500,
													  'b' => 250,
													  'g' => 'T',
													  'W' => 11,
													  'T' => 'F',
													  'n' => 'F',
													  'K' => 0,													  
													  'a' => 4
													 },
						'genomic_blat' => 0,
						'mrna_blat' => 0,
						'wrap_dir' =>"$ENV{'TrapCluster'}/wrap" #Directory where shell scripts for cluster jobs are written
		
		},
		#Configuration for which annotations should be performed
		'annotation' =>{
						'do'=>1, #Perform annotation (if 0, only mapping)
						'do_ensgene' => 1, #Check Ensembl gene annotation
						'do_ensestgene' => 1, #Check Ensembl EST gene annotation
						'do_genescan' => 1, #Check GENSCAN annotation
						'do_unigene' => 1, #Check Unigene
						'do_mouse_cdna' => 1, #Check Mouse cDNAs
						'do_ensest' => 1, #Check Ensembl ESTs
						'do_tclg' => 0, #Build Trap Cluster Genes ***NOT IMPLEMENTED YET***
						'do_ensrepeats' => 1, #Check repeats
						'unigene_coverage' => 0,
						'unigene_perc_id' => 96,
						'cdna_coverage' => 66,
						'cdna_perc_id' => 96,
						'est_coverage' => 66,
						'est_perc_id' => 96
		}

    );
}
1;
