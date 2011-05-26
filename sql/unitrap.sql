-- MySQL dump 10.11
--
-- Host: localhost    Database: unitrap
-- ------------------------------------------------------
-- Server version	5.1.51

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `analysis`
--

DROP TABLE IF EXISTS `analysis`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `analysis` (
  `analysis_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `created` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `logic_name_id` int(10) unsigned NOT NULL,
  `db` varchar(255) DEFAULT NULL,
  `db_version` varchar(255) DEFAULT NULL,
  `db_file` varchar(255) DEFAULT NULL,
  `program` varchar(255) DEFAULT NULL,
  `program_version` varchar(255) DEFAULT NULL,
  `parameters` varchar(255) DEFAULT NULL,
  `description` text,
  PRIMARY KEY (`analysis_id`),
  KEY `logic_name_id_idxfk` (`logic_name_id`)
) ENGINE=MyISAM AUTO_INCREMENT=2523813 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `esclone`
--

DROP TABLE IF EXISTS `esclone`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `esclone` (
  `esclone_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `design_type` varchar(255) DEFAULT NULL,
  `clone` varchar(255) DEFAULT NULL,
  `project_id` int(10) unsigned NOT NULL DEFAULT '0',
  `clone_genbank_id` int(10) unsigned DEFAULT NULL,
  `cell_line` varchar(40) DEFAULT NULL,
  `vector_name` varchar(40) NOT NULL,
  `vector_type` varchar(40) NOT NULL,
  `strain` varchar(40) NOT NULL,
  PRIMARY KEY (`esclone_id`),
  UNIQUE KEY `design_type` (`design_type`,`clone`,`project_id`,`clone_genbank_id`,`cell_line`,`strain`),
  KEY `design_type_idx` (`design_type`),
  KEY `clone_idx` (`clone`),
  KEY `clone_genbank_id_idx` (`clone_genbank_id`),
  KEY `project_id` (`project_id`),
  FULLTEXT KEY `vector_name` (`vector_name`,`vector_type`)
) ENGINE=MyISAM AUTO_INCREMENT=717726 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `genbank_features`
--

DROP TABLE IF EXISTS `genbank_features`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `genbank_features` (
  `genbank_features_id` int(22) unsigned NOT NULL AUTO_INCREMENT,
  `tagname` varchar(255) DEFAULT NULL,
  `esclone_id` int(11) DEFAULT NULL,
  `Comment` text NOT NULL,
  PRIMARY KEY (`genbank_features_id`),
  UNIQUE KEY `tagname` (`tagname`,`esclone_id`),
  KEY `genbank_features_id_idx` (`genbank_features_id`),
  KEY `tagname_idx` (`tagname`),
  KEY `esclone_id_idx` (`esclone_id`),
  FULLTEXT KEY `Comment` (`Comment`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `insertion`
--

DROP TABLE IF EXISTS `insertion`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `insertion` (
  `insertion_id` int(11) NOT NULL AUTO_INCREMENT,
  `trap_id` int(11) NOT NULL,
  `trapmap_id` int(11) NOT NULL,
  `trapblock_id` int(11) NOT NULL,
  `gene_id` int(11) NOT NULL,
  `putative_insertion_start` int(10) DEFAULT NULL,
  `putative_insertion_end` int(10) DEFAULT NULL,
  `insertion_case` int(1) DEFAULT NULL COMMENT 'logic_name_id',
  `insertion_ambiguous` tinyint(1) NOT NULL DEFAULT '0',
  `region_start` int(11) NOT NULL,
  `region_end` int(11) NOT NULL,
  PRIMARY KEY (`insertion_id`),
  UNIQUE KEY `type` (`trap_id`,`trapmap_id`,`trapblock_id`,`gene_id`),
  KEY `insertion_case` (`insertion_case`),
  KEY `insertion_ambiguous` (`insertion_ambiguous`),
  KEY `trap_id` (`trap_id`,`trapmap_id`,`trapblock_id`),
  KEY `transcript_id` (`gene_id`),
  KEY `region_start` (`region_start`,`region_end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `logic_name`
--

DROP TABLE IF EXISTS `logic_name`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `logic_name` (
  `logic_name_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  PRIMARY KEY (`logic_name_id`),
  UNIQUE KEY `logic_name_id` (`logic_name_id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=MyISAM AUTO_INCREMENT=28 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `mutated_gene`
--

DROP TABLE IF EXISTS `mutated_gene`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `mutated_gene` (
  `mutated_gene_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_id` int(10) unsigned NOT NULL,
  `gene_type` varchar(40) DEFAULT NULL,
  `unitrap_num` int(4) NOT NULL DEFAULT '0',
  `transcript_num` int(4) NOT NULL DEFAULT '0',
  `finished` tinyint(1) NOT NULL DEFAULT '0',
  `word` text,
  `disease` mediumtext,
  `human_ortholog` tinyint(1) NOT NULL DEFAULT '0',
  `omim` tinyint(1) NOT NULL DEFAULT '0',
  `mutagenicity` varchar(20) NOT NULL,
  `public` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`mutated_gene_id`),
  UNIQUE KEY `region_id` (`region_id`,`gene_type`,`unitrap_num`,`transcript_num`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `mutated_gene_exon`
--

DROP TABLE IF EXISTS `mutated_gene_exon`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `mutated_gene_exon` (
  `mutated_gene_exon_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_id` int(10) unsigned DEFAULT NULL,
  `exon_seq` mediumtext,
  `exon_rank` int(4) NOT NULL DEFAULT '0',
  `mutated_gene_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`mutated_gene_exon_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `mutated_gene_exon_primer`
--

DROP TABLE IF EXISTS `mutated_gene_exon_primer`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `mutated_gene_exon_primer` (
  `mutated_gene_exon_primer_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `mutated_gene_exon_id` int(10) unsigned NOT NULL DEFAULT '0',
  `name` varchar(255) NOT NULL,
  `sequence` mediumtext NOT NULL,
  `start` int(10) NOT NULL DEFAULT '0',
  `end` int(10) NOT NULL DEFAULT '0',
  `length` int(10) NOT NULL DEFAULT '0',
  `genomic_start` int(10) NOT NULL DEFAULT '0',
  `genomic_end` int(10) NOT NULL DEFAULT '0',
  `genomic_strand` int(10) NOT NULL DEFAULT '0',
  `gc_perc` float NOT NULL DEFAULT '0',
  `tm` float NOT NULL DEFAULT '0',
  `type` varchar(255) NOT NULL,
  `best_primer` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`mutated_gene_exon_primer_id`),
  KEY `mutated_gene_exon_id_idxfk` (`mutated_gene_exon_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `mutated_gene_fragment`
--

DROP TABLE IF EXISTS `mutated_gene_fragment`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `mutated_gene_fragment` (
  `mutated_gene_fragment_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `mutated_gene_id` int(10) unsigned NOT NULL,
  `start` int(10) NOT NULL DEFAULT '0',
  `end` int(10) NOT NULL DEFAULT '0',
  `chr` varchar(40) NOT NULL,
  `enzyme` varchar(20) NOT NULL,
  `restriction_site` varchar(20) NOT NULL,
  `sequence` mediumtext NOT NULL,
  PRIMARY KEY (`mutated_gene_fragment_id`),
  UNIQUE KEY `mutated_gene_id` (`mutated_gene_id`,`start`,`end`,`chr`,`enzyme`,`restriction_site`),
  KEY `mutated_gene_id_idxfk` (`mutated_gene_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `mutated_transcript`
--

DROP TABLE IF EXISTS `mutated_transcript`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `mutated_transcript` (
  `mutated_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_id` int(10) unsigned DEFAULT NULL,
  `transcript_seq` mediumtext,
  `peptide_seq` mediumtext,
  `exon_num` int(4) NOT NULL DEFAULT '0',
  `uni_transcript_seq` mediumtext,
  `uni_peptide_seq` mediumtext,
  `mutated_gene_id` int(11) DEFAULT NULL,
  PRIMARY KEY (`mutated_transcript_id`),
  UNIQUE KEY `region_id` (`region_id`,`exon_num`,`mutated_gene_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `mutated_transcript_exon`
--

DROP TABLE IF EXISTS `mutated_transcript_exon`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `mutated_transcript_exon` (
  `mutated_transcript_exon_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_id` varchar(40) DEFAULT NULL,
  `mutated_transcript_id` int(10) unsigned DEFAULT NULL,
  `exon_seq` mediumtext,
  `peptide_seq` mediumtext,
  `exon_rank` int(4) NOT NULL DEFAULT '0',
  PRIMARY KEY (`mutated_transcript_exon_id`),
  UNIQUE KEY `region_id` (`region_id`,`mutated_transcript_id`,`exon_rank`),
  KEY `mutated_transcript_id_idxfk` (`mutated_transcript_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `omim`
--

DROP TABLE IF EXISTS `omim`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `omim` (
  `omim_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ortholog_id` int(10) unsigned NOT NULL,
  `accession` varchar(40) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`omim_id`),
  UNIQUE KEY `ortholog_id` (`ortholog_id`,`accession`,`description`),
  KEY `ortholog_id_idxfk` (`ortholog_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `ortholog`
--

DROP TABLE IF EXISTS `ortholog`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `ortholog` (
  `ortholog_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `hs_ensembl_id` varchar(40) DEFAULT NULL,
  `region_id` int(10) unsigned NOT NULL,
  `ort_name` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`ortholog_id`),
  UNIQUE KEY `hs_ensembl_id` (`hs_ensembl_id`,`region_id`,`ort_name`),
  KEY `region_id_idxfk` (`region_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `pcr`
--

DROP TABLE IF EXISTS `pcr`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `pcr` (
  `pcr_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `unitrap_id` int(10) unsigned NOT NULL DEFAULT '0',
  `hit_id` varchar(40) NOT NULL DEFAULT '0',
  `forward` varchar(40) NOT NULL DEFAULT '0',
  `forward_start` int(10) NOT NULL DEFAULT '0',
  `forward_strand` enum('+','-') DEFAULT NULL,
  `reverse` varchar(40) NOT NULL DEFAULT '0',
  `gc_perc` float DEFAULT NULL,
  `tm` float DEFAULT NULL,
  `reverse_name` varchar(10) DEFAULT NULL,
  `forward_end` int(10) NOT NULL DEFAULT '0',
  PRIMARY KEY (`pcr_id`),
  UNIQUE KEY `unitrap_id` (`unitrap_id`,`hit_id`,`forward`,`forward_start`,`forward_strand`,`reverse`,`gc_perc`,`tm`,`reverse_name`,`forward_end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `project`
--

DROP TABLE IF EXISTS `project`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `project` (
  `project_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `project_name` varchar(40) NOT NULL,
  `url` varchar(250) NOT NULL,
  `wordkey` varchar(255) NOT NULL,
  `private` tinyint(1) NOT NULL DEFAULT '0',
  `note` text NOT NULL,
  PRIMARY KEY (`project_id`)
) ENGINE=MyISAM AUTO_INCREMENT=16 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `protein_mutagenicity`
--

DROP TABLE IF EXISTS `protein_mutagenicity`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `protein_mutagenicity` (
  `protein_mutagenicity_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `unitrap_id` int(10) unsigned NOT NULL DEFAULT '0',
  `protein_seq` mediumtext,
  `mutagenicity` varchar(50) DEFAULT NULL,
  `not_mutated_exons_num` int(10) NOT NULL DEFAULT '0',
  `mutated_exons_num` int(10) NOT NULL DEFAULT '0',
  `not_mutated_domains_num` int(10) NOT NULL DEFAULT '0',
  `mutated_domains_num` int(10) NOT NULL DEFAULT '0',
  `not_mutated_domains` mediumtext,
  `mutated_domains` mediumtext,
  `rank` tinyint(1) NOT NULL DEFAULT '0',
  `longest_transcript` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`protein_mutagenicity_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `region`
--

DROP TABLE IF EXISTS `region`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `region` (
  `region_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_name` varchar(40) NOT NULL,
  `seq_id` varchar(255) NOT NULL,
  `region_strand` enum('-1','1') DEFAULT NULL,
  `region_start` int(10) unsigned DEFAULT NULL,
  `region_end` int(10) unsigned DEFAULT NULL,
  `refseq_id` varchar(40) NOT NULL,
  `external_name` varchar(255) DEFAULT NULL,
  `biotype` int(11) DEFAULT NULL COMMENT 'logic_name_id',
  `description` varchar(255) NOT NULL,
  `parent_id` int(11) DEFAULT NULL,
  `rank` int(11) DEFAULT NULL,
  PRIMARY KEY (`region_id`),
  UNIQUE KEY `region_name` (`seq_id`,`region_start`,`region_end`,`region_strand`,`region_name`),
  KEY `region_id_idx` (`region_id`),
  KEY `seq_id_idx` (`seq_id`),
  KEY `region_start_idx` (`region_start`),
  KEY `region_end_idx` (`region_end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `region_go`
--

DROP TABLE IF EXISTS `region_go`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `region_go` (
  `region_go_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_id` int(10) unsigned NOT NULL,
  `accession` varchar(20) NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  `source` varchar(20) DEFAULT NULL,
  `go_class` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`region_go_id`),
  UNIQUE KEY `region_id` (`region_id`,`accession`,`description`,`source`,`go_class`),
  KEY `region_id_idxfk` (`region_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `trap`
--

DROP TABLE IF EXISTS `trap`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `trap` (
  `trap_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `trap_name` varchar(40) NOT NULL,
  `sequence` mediumtext NOT NULL,
  `timestp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `gb_id` int(10) DEFAULT NULL,
  `gb_locus` varchar(30) DEFAULT NULL,
  `gss_id` int(11) NOT NULL DEFAULT '0',
  `mol_type` varchar(40) DEFAULT NULL,
  `seq_length` int(2) unsigned DEFAULT '0',
  `seq_length_not_N` int(2) unsigned DEFAULT '0',
  `max_frag_length_N_splitted` int(2) unsigned DEFAULT '0',
  `x_masked_seq` mediumtext,
  `nrepeat` int(2) DEFAULT '0',
  `xrepeat` int(2) DEFAULT '0',
  `paired_tag_id` int(10) unsigned NOT NULL,
  `x_percent_masked` int(2) DEFAULT '0',
  `n_percent_masked` int(2) DEFAULT '0',
  `sequencing_direction` enum('5','3','na','I','R','F') DEFAULT 'na',
  `esclone_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`trap_id`),
  UNIQUE KEY `trap_name` (`trap_name`,`gb_id`,`gb_locus`,`gss_id`,`mol_type`,`seq_length`,`seq_length_not_N`,`max_frag_length_N_splitted`,`nrepeat`,`xrepeat`,`paired_tag_id`,`x_percent_masked`,`n_percent_masked`,`sequencing_direction`,`esclone_id`),
  KEY `paired_tag_id_idx` (`paired_tag_id`),
  KEY `esclone_id_idx` (`esclone_id`)
) ENGINE=MyISAM AUTO_INCREMENT=902481 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `trap_error`
--

DROP TABLE IF EXISTS `trap_error`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `trap_error` (
  `trap_error_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `logic_name_id` int(10) unsigned NOT NULL,
  `table_id` int(10) unsigned NOT NULL DEFAULT '0',
  `note` text NOT NULL,
  `table_name` varchar(255) DEFAULT NULL,
  `timestamp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`trap_error_id`),
  UNIQUE KEY `table_id` (`table_id`,`table_name`,`logic_name_id`),
  KEY `table_id_idxfk` (`table_id`)
) ENGINE=MyISAM AUTO_INCREMENT=7991 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `trap_unitrap`
--

DROP TABLE IF EXISTS `trap_unitrap`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `trap_unitrap` (
  `trap_unitrap_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `trap_id` int(10) unsigned DEFAULT NULL,
  `trapblock_id` int(10) NOT NULL DEFAULT '0',
  `unitrap_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`trap_unitrap_id`),
  UNIQUE KEY `tum` (`trap_id`,`unitrap_id`,`trapblock_id`),
  KEY `trap_id_idxfk` (`trap_id`),
  KEY `unitrap_id_idxfk` (`unitrap_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `trapadditional`
--

DROP TABLE IF EXISTS `trapadditional`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `trapadditional` (
  `trapadditional_id` int(11) NOT NULL AUTO_INCREMENT,
  `processing` varchar(255) DEFAULT NULL,
  `quality` varchar(255) DEFAULT NULL,
  `comment` varchar(255) DEFAULT NULL,
  `verified` tinyint(1) DEFAULT NULL,
  `label` varchar(255) DEFAULT NULL,
  `tigem_ex_ref` varchar(40) NOT NULL COMMENT 'old tigem trap name',
  `trap_id` int(10) unsigned DEFAULT NULL,
  `user` varchar(255) DEFAULT NULL,
  `note` mediumtext,
  PRIMARY KEY (`trapadditional_id`),
  UNIQUE KEY `pair` (`trap_id`,`label`),
  KEY `trap_id` (`trap_id`),
  KEY `processing` (`processing`),
  KEY `quality` (`quality`),
  KEY `verified` (`verified`),
  FULLTEXT KEY `label` (`label`),
  FULLTEXT KEY `comment` (`comment`),
  FULLTEXT KEY `note` (`note`),
  FULLTEXT KEY `user` (`user`)
) ENGINE=MyISAM AUTO_INCREMENT=5 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `trapblock`
--

DROP TABLE IF EXISTS `trapblock`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `trapblock` (
  `trapblock_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `trapmap_id` int(10) unsigned NOT NULL DEFAULT '0',
  `start` int(10) NOT NULL DEFAULT '0',
  `end` int(10) NOT NULL DEFAULT '0',
  `strand` enum('1','-1') NOT NULL,
  `perc_ident` float DEFAULT '0',
  `accepted` tinyint(4) DEFAULT NULL,
  PRIMARY KEY (`trapblock_id`),
  UNIQUE KEY `trapmap_id` (`trapmap_id`,`start`,`end`,`strand`,`perc_ident`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `trapblock_annotation`
--

DROP TABLE IF EXISTS `trapblock_annotation`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `trapblock_annotation` (
  `trapblock_annotation_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `source_logic_name_id` int(11) DEFAULT NULL,
  `trapblock_id` int(10) unsigned NOT NULL DEFAULT '0',
  `rank` float NOT NULL DEFAULT '0',
  `display_name` varchar(40) DEFAULT NULL,
  `dS_s` int(10) DEFAULT NULL,
  `dE_e` int(10) DEFAULT NULL,
  `label_logic_name_id` int(11) DEFAULT NULL,
  `comment` text,
  `coverage` float DEFAULT NULL,
  `region_id` int(10) unsigned DEFAULT NULL,
  `flanking_exon_id` int(10) unsigned DEFAULT NULL,
  `type` tinyint(1) NOT NULL,
  PRIMARY KEY (`trapblock_annotation_id`),
  UNIQUE KEY `unique` (`trapblock_id`,`rank`,`display_name`,`dS_s`,`dE_e`,`label_logic_name_id`,`coverage`,`region_id`,`source_logic_name_id`,`comment`(350)),
  KEY `source_logic_name_id_idxfk` (`source_logic_name_id`),
  KEY `trapblock_id_idxfk` (`trapblock_id`),
  KEY `label_logic_name_id_idxfk` (`label_logic_name_id`),
  KEY `region_id_idxfk` (`region_id`),
  KEY `flanking_exon_id` (`flanking_exon_id`),
  KEY `intronic` (`type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `trapcheck`
--

DROP TABLE IF EXISTS `trapcheck`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `trapcheck` (
  `checkid` int(11) NOT NULL AUTO_INCREMENT,
  `trap_id` int(10) unsigned DEFAULT NULL,
  `mapped` tinyint(1) NOT NULL DEFAULT '0',
  `annotated` int(11) NOT NULL,
  `checked` tinyint(1) NOT NULL DEFAULT '0',
  `mapping_type` tinyint(1) DEFAULT NULL,
  `splk` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`checkid`),
  UNIQUE KEY `trap_id` (`trap_id`),
  UNIQUE KEY `trap_id_2` (`trap_id`,`mapped`,`checked`,`mapping_type`,`splk`),
  KEY `annotated` (`annotated`)
) ENGINE=MyISAM AUTO_INCREMENT=902481 DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `trapmap`
--

DROP TABLE IF EXISTS `trapmap`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `trapmap` (
  `trapmap_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `trap_id` int(10) unsigned NOT NULL DEFAULT '0',
  `hit_id` varchar(40) NOT NULL,
  `hit_db` varchar(40) NOT NULL,
  `start` int(10) NOT NULL DEFAULT '0',
  `end` int(10) NOT NULL DEFAULT '0',
  `strand` enum('1','-1') NOT NULL,
  `frac_aligned_query` float DEFAULT NULL,
  `num_hsps` int(10) NOT NULL DEFAULT '0',
  `frac_identical` float DEFAULT NULL,
  `score` float NOT NULL DEFAULT '0',
  `significance` double DEFAULT NULL,
  `chosen` tinyint(1) NOT NULL DEFAULT '0',
  `multiple` tinyint(1) DEFAULT NULL,
  `chosen_filter` tinyint(1) NOT NULL DEFAULT '0',
  `checked` tinyint(1) NOT NULL DEFAULT '0',
  `target_type` int(11) NOT NULL COMMENT 'logic_name ',
  `analysis_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`trapmap_id`),
  UNIQUE KEY `trap_id_2` (`trap_id`,`hit_id`,`hit_db`,`start`,`end`,`frac_aligned_query`,`num_hsps`,`frac_identical`,`score`,`significance`,`target_type`),
  KEY `target_type` (`target_type`),
  KEY `trap_id` (`trap_id`),
  KEY `analysis_id` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `trapmap_region`
--

DROP TABLE IF EXISTS `trapmap_region`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `trapmap_region` (
  `trapmap_region_id` int(11) NOT NULL AUTO_INCREMENT,
  `region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `type` varchar(40) NOT NULL,
  `annotation_ambiguous` tinyint(1) NOT NULL DEFAULT '1',
  `total_number_exons` int(11) NOT NULL,
  `number_trapblocks` int(10) unsigned DEFAULT NULL,
  `number_annotated_trapblocks` int(10) unsigned DEFAULT NULL,
  `overlap` int(10) NOT NULL,
  `trapmap_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`trapmap_region_id`),
  UNIQUE KEY `type` (`type`,`annotation_ambiguous`,`number_trapblocks`,`number_annotated_trapblocks`,`overlap`,`trapmap_id`,`region_id`),
  KEY `region_id_idx` (`region_id`),
  KEY `trapmap_id_idx` (`trapmap_id`),
  KEY `total_number_exons` (`total_number_exons`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `traprepeat`
--

DROP TABLE IF EXISTS `traprepeat`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `traprepeat` (
  `traprepeat_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `trap_id` int(10) unsigned DEFAULT NULL,
  `start` int(10) NOT NULL DEFAULT '0',
  `end` int(10) NOT NULL DEFAULT '0',
  `score` float DEFAULT NULL,
  `tag` varchar(25) DEFAULT NULL,
  PRIMARY KEY (`traprepeat_id`),
  UNIQUE KEY `trap_id` (`trap_id`,`start`,`end`,`score`,`tag`),
  KEY `trap_id_idxfk` (`trap_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `unitrap`
--

DROP TABLE IF EXISTS `unitrap`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `unitrap` (
  `unitrap_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `accession` varchar(40) NOT NULL,
  `hit_db` varchar(40) NOT NULL,
  `chr` varchar(40) NOT NULL DEFAULT '0',
  `start` int(10) NOT NULL DEFAULT '0',
  `end` int(10) NOT NULL DEFAULT '0',
  `ambiguous` tinyint(1) NOT NULL DEFAULT '0',
  `region_id` int(11) NOT NULL COMMENT 'gene_id',
  PRIMARY KEY (`unitrap_id`),
  UNIQUE KEY `accession_2` (`accession`),
  UNIQUE KEY `accession` (`accession`,`hit_db`,`chr`,`region_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `unitrap_history`
--

DROP TABLE IF EXISTS `unitrap_history`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `unitrap_history` (
  `history_id` int(11) NOT NULL AUTO_INCREMENT,
  `old_unitrap_id` varchar(255) NOT NULL,
  `new_unitrap_id` varchar(255) NOT NULL,
  PRIMARY KEY (`history_id`),
  KEY `old_unitrap_id` (`old_unitrap_id`,`new_unitrap_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='historic link between unitrap accession on the web site and ';
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `xref`
--

DROP TABLE IF EXISTS `xref`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `xref` (
  `xref_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_id` int(10) unsigned NOT NULL,
  `dbname` varchar(40) DEFAULT NULL,
  `accession` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`xref_id`),
  UNIQUE KEY `region_id` (`region_id`,`dbname`,`accession`),
  KEY `region_id_idxfk` (`region_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2010-12-06 14:08:03
