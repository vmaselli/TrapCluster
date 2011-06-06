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
						'ensdbname' => "mus_musculus_core_57_37j",
						'ensestdbname' => "mus_musculus_otherfeatures_57_37j", 
						'ens_hs_dbname' => "homo_sapiens_core_57_37b", 
						'martdbname' => "ensembl_mart_57",
						'host' => "localhost",
						'user' => "root",
						#'dbname' => "UniTrapPipeline",
						'dbname' => 'unitrap',
						'pass' => "",
						'port' => 3306
					},
        # default variables 
	    'default' =>{
						'hit_db' => "NCBIM37",#Version of the mouse genome to use
						'ensembl_version' =>"59",
						'ensembl_base_url' => "http://www.ensembl.org/Mus_musculus/",
						'ucsc_db_name' => "mm9", #UCSC mouse db version
						'tmp_dir' => "/tmp/",
						'mysql_path' => '/usr/local/mysql/bin/',
						'version' => 'ab-blast 2.2.6'
					},
        # user defined variables
	    'global' =>{
						'debug' => 1,
						'debugSQL' => 1,
						'swarm' => 1, # set 1 to run in parallel mode
					},
		'retrieve' =>{
						'update' => 0, #Set to 1 if you are performing an update on an existing db
						'fromfile' => 0, #Set to 1 if starting data is from files (Genbank, FASTA, GSS) 
						'fromquery' => 1, #Set to 1 if starting from an NCBI query
						'mindate' => "2000/01/01", #Date from which NCBI query should start
						'fromdb' => 0, #Set to 1 to start from a Unitrap DB already populated with traps
						'fromcsv' =>1, #Set to 1 to start from CSV files 
						'dir' => "$ENV{'HOME'}/Data/Projects/Trap/GenBank", # Genbank file directory
						#'dir' => "$ENV{'HOME'}/Data/Projects/Trap/Test/GenBank", # test genbank files directory
						#'dir' => "$ENV{'HOME'}/Data/Projects/Trap/dbGSS", # gss files directory
						'format' => 'GenBank',
						#'format' => 'dbGSS',
						#'format' =>'fasta',
						#'format' => 'csv',
						'base_keyword' => "gene trap and mouse AND GSS and ",
						#This should move to the API
						'query' => "SELECT t.*, p.* FROM `trap` t, `project` p WHERE t.project_id = p.project_id ",
						'last_id' => "SELECT trap_id FROM trap order by trap_id DESC limit 1"
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
						'wrap_dir' =>"$ENV{'Unitrap'}/wrap" #Directory where shell scripts for cluster jobs are written
		
		},
		#Configuration for which annotations should be performed
		'annotation' =>{
						'do'=>1, #Perform annotation (if 0, only mapping)
						'do_ensgene' => 1, #Check Ensembl gene annotation
						'do_ensestgene' => 0, #Check Ensembl EST gene annotation
						'do_genescan' => 0, #Check GENSCAN annotation
						'do_unigene' => 0, #Check Unigene
						'do_mouse_cdna' => 0, #Check Mouse cDNAs
						'do_ensest' => 0, #Check Ensembl ESTs
						'do_tclg' => 0, #Build Trap Cluster Genes ***NOT IMPLEMENTED YET***
						'do_ensrepeats' => 0 #Check repeats
		}

    );
}
1;
