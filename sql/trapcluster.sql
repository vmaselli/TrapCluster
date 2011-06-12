-- phpMyAdmin SQL Dump
-- version 3.2.2.1deb1
-- http://www.phpmyadmin.net
--
-- Host: localhost
-- Generation Time: Jun 09, 2011 at 11:32 AM
-- Server version: 5.1.37
-- PHP Version: 5.2.10-2ubuntu6.10

SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: `trapcluster`
--

-- --------------------------------------------------------

--
-- Table structure for table `analysis`
--

CREATE TABLE IF NOT EXISTS `analysis` (
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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=956026 ;

-- --------------------------------------------------------

--
-- Table structure for table `esclone`
--

CREATE TABLE IF NOT EXISTS `esclone` (
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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=1010170 ;

-- --------------------------------------------------------

--
-- Table structure for table `genbank_features`
--

CREATE TABLE IF NOT EXISTS `genbank_features` (
  `genbank_features_id` int(22) unsigned NOT NULL AUTO_INCREMENT,
  `tagname` varchar(255) DEFAULT NULL,
  `esclone_id` int(11) DEFAULT NULL,
  `Comment` text NOT NULL,
  PRIMARY KEY (`genbank_features_id`),
  UNIQUE KEY `tagname` (`tagname`,`esclone_id`),
  KEY `genbank_features_id_idx` (`genbank_features_id`),
  KEY `tagname_idx` (`tagname`),
  KEY `esclone_id_idx` (`esclone_id`),
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `insertion`
--

CREATE TABLE IF NOT EXISTS `insertion` (
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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=189078 ;

-- --------------------------------------------------------

--
-- Table structure for table `logic_name`
--

CREATE TABLE IF NOT EXISTS `logic_name` (
  `logic_name_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  PRIMARY KEY (`logic_name_id`),
  UNIQUE KEY `logic_name_id` (`logic_name_id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=56 ;

-- --------------------------------------------------------

--
-- Table structure for table `maxicluster`
--

CREATE TABLE IF NOT EXISTS `maxicluster` (
  `maxicluster_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `sequence` mediumtext NOT NULL,
  `accession` varchar(40) NOT NULL DEFAULT '',
  `mol_type` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`maxicluster_id`),
  UNIQUE KEY `acc` (`accession`),
  KEY `mol` (`mol_type`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=2 ;

-- --------------------------------------------------------

--
-- Table structure for table `maxiclusterblock`
--

CREATE TABLE IF NOT EXISTS `maxiclusterblock` (
  `maxiclusterblock_id` int(11) NOT NULL AUTO_INCREMENT,
  `maxiclustermap_id` int(11) NOT NULL,
  `start` int(11) NOT NULL,
  `end` int(11) NOT NULL,
  `strand` enum('1','-1') NOT NULL,
  `sequence` text NOT NULL,
  `checked` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`maxiclusterblock_id`),
  KEY `maxiclustermap_id` (`maxiclustermap_id`,`start`,`end`),
  KEY `checked` (`checked`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=55 ;

-- --------------------------------------------------------

--
-- Table structure for table `maxiclustermap`
--

CREATE TABLE IF NOT EXISTS `maxiclustermap` (
  `maxiclustermap_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `hit_id` varchar(40) NOT NULL,
  `hit_db` varchar(40) NOT NULL,
  `start` int(10) NOT NULL DEFAULT '0',
  `end` int(10) NOT NULL DEFAULT '0',
  `strand` enum('1','-1') NOT NULL,
  `accession` varchar(255) NOT NULL,
  `analysis_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`maxiclustermap_id`),
  UNIQUE KEY `trap_id_2` (`hit_id`,`hit_db`,`start`,`end`),
  KEY `analysis_id` (`analysis_id`),
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=6 ;

-- --------------------------------------------------------

--
-- Table structure for table `mutated_gene`
--

CREATE TABLE IF NOT EXISTS `mutated_gene` (
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
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `mutated_gene_exon`
--

CREATE TABLE IF NOT EXISTS `mutated_gene_exon` (
  `mutated_gene_exon_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_id` int(10) unsigned DEFAULT NULL,
  `exon_seq` mediumtext,
  `exon_rank` int(4) NOT NULL DEFAULT '0',
  `mutated_gene_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`mutated_gene_exon_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `mutated_gene_exon_primer`
--

CREATE TABLE IF NOT EXISTS `mutated_gene_exon_primer` (
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
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `mutated_gene_fragment`
--

CREATE TABLE IF NOT EXISTS `mutated_gene_fragment` (
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
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `mutated_transcript`
--

CREATE TABLE IF NOT EXISTS `mutated_transcript` (
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
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `mutated_transcript_exon`
--

CREATE TABLE IF NOT EXISTS `mutated_transcript_exon` (
  `mutated_transcript_exon_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_id` varchar(40) DEFAULT NULL,
  `mutated_transcript_id` int(10) unsigned DEFAULT NULL,
  `exon_seq` mediumtext,
  `peptide_seq` mediumtext,
  `exon_rank` int(4) NOT NULL DEFAULT '0',
  PRIMARY KEY (`mutated_transcript_exon_id`),
  UNIQUE KEY `region_id` (`region_id`,`mutated_transcript_id`,`exon_rank`),
  KEY `mutated_transcript_id_idxfk` (`mutated_transcript_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `omim`
--

CREATE TABLE IF NOT EXISTS `omim` (
  `omim_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ortholog_id` int(10) unsigned NOT NULL,
  `accession` varchar(40) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`omim_id`),
  UNIQUE KEY `ortholog_id` (`ortholog_id`,`accession`,`description`),
  KEY `ortholog_id_idxfk` (`ortholog_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `ortholog`
--

CREATE TABLE IF NOT EXISTS `ortholog` (
  `ortholog_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `hs_ensembl_id` varchar(40) DEFAULT NULL,
  `region_id` int(10) unsigned NOT NULL,
  `ort_name` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`ortholog_id`),
  UNIQUE KEY `hs_ensembl_id` (`hs_ensembl_id`,`region_id`,`ort_name`),
  KEY `region_id_idxfk` (`region_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `pcr`
--

CREATE TABLE IF NOT EXISTS `pcr` (
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
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `project`
--

CREATE TABLE IF NOT EXISTS `project` (
  `project_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `project_name` varchar(40) NOT NULL,
  `url` varchar(250) NOT NULL,
  `wordkey` varchar(255) NOT NULL,
  `private` tinyint(1) NOT NULL DEFAULT '0',
  `note` text NOT NULL,
  PRIMARY KEY (`project_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=24 ;

-- --------------------------------------------------------

--
-- Table structure for table `protein_mutagenicity`
--

CREATE TABLE IF NOT EXISTS `protein_mutagenicity` (
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
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `region`
--

CREATE TABLE IF NOT EXISTS `region` (
  `region_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_name` varchar(40) NOT NULL,
  `seq_id` varchar(255) DEFAULT NULL,
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
  UNIQUE KEY `region_name` (`seq_id`,`region_start`,`region_end`,`region_strand`,`region_name`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=774063 ;

-- --------------------------------------------------------

--
-- Table structure for table `region_go`
--

CREATE TABLE IF NOT EXISTS `region_go` (
  `region_go_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_id` int(10) unsigned NOT NULL,
  `accession` varchar(20) NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  `source` varchar(20) DEFAULT NULL,
  `go_class` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`region_go_id`),
  UNIQUE KEY `region_id` (`region_id`,`accession`,`description`,`source`,`go_class`),
  KEY `region_id_idxfk` (`region_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `trap`
--

CREATE TABLE IF NOT EXISTS `trap` (
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
  `sequencing_direction` enum('5','3','na','I','R','F','L','R','5S','3S','5SPK','3SPK') DEFAULT NULL,
  `esclone_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`trap_id`),
  UNIQUE KEY `trap_name` (`trap_name`,`gb_id`,`gb_locus`,`gss_id`,`mol_type`,`seq_length`,`seq_length_not_N`,`max_frag_length_N_splitted`,`nrepeat`,`xrepeat`,`paired_tag_id`,`x_percent_masked`,`n_percent_masked`,`sequencing_direction`,`esclone_id`),
  KEY `paired_tag_id_idx` (`paired_tag_id`),
  KEY `esclone_id_idx` (`esclone_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=559007 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapadditional`
--

CREATE TABLE IF NOT EXISTS `trapadditional` (
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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=700175 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapblock`
--

CREATE TABLE IF NOT EXISTS `trapblock` (
  `trapblock_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `trapmap_id` int(10) unsigned NOT NULL DEFAULT '0',
  `start` int(10) NOT NULL DEFAULT '0',
  `end` int(10) NOT NULL DEFAULT '0',
  `strand` enum('1','-1') NOT NULL,
  `perc_ident` float DEFAULT '0',
  `accepted` tinyint(4) DEFAULT NULL,
  PRIMARY KEY (`trapblock_id`),
  UNIQUE KEY `trapmap_id` (`trapmap_id`,`start`,`end`,`strand`,`perc_ident`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=1266998 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapblock_annotation`
--

CREATE TABLE IF NOT EXISTS `trapblock_annotation` (
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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=189078 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapcheck`
--

CREATE TABLE IF NOT EXISTS `trapcheck` (
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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=292545 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapcluster`
--

CREATE TABLE IF NOT EXISTS `trapcluster` (
  `trapcluster_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `accession` varchar(40) NOT NULL DEFAULT '',
  `overlapping` varchar(40) DEFAULT NULL,
  `maxicluster_id` int(10) NOT NULL,
  `sequence` mediumtext,
  `trapclusterblocks` int(10) unsigned DEFAULT '0',
  `traps` int(10) unsigned DEFAULT '0',
  `checked` tinyint(1) DEFAULT '0',
  `hypertrapped` tinyint(1) DEFAULT '0',
  `intronic` tinyint(1) DEFAULT '0',
  `link_to_ensembl` mediumtext,
  `link_to_ucsc` mediumtext,
  `flanking` tinyint(1) DEFAULT '0',
  `isolated` tinyint(1) DEFAULT '0',
  `checked_flanking` tinyint(1) DEFAULT '0',
  `antisense` tinyint(1) DEFAULT '0',
  PRIMARY KEY (`trapcluster_id`),
  UNIQUE KEY `acc` (`accession`),
  KEY `overlapping` (`overlapping`),
  KEY `isolated` (`isolated`),
  KEY `antisense` (`antisense`),
  KEY `checked_flanking` (`checked_flanking`),
  KEY `checked` (`checked`),
  KEY `flanking` (`flanking`),
  KEY `hypertrapped` (`hypertrapped`),
  KEY `traps` (`traps`),
  KEY `trapclusterblocks` (`trapclusterblocks`),
  KEY `intronic` (`intronic`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=6 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapclusteradditional`
--

CREATE TABLE IF NOT EXISTS `trapclusteradditional` (
  `trapclusteradditional_id` int(11) NOT NULL AUTO_INCREMENT,
  `processing` varchar(255) DEFAULT NULL,
  `quality` varchar(255) DEFAULT NULL,
  `comment` varchar(255) DEFAULT NULL,
  `verified` tinyint(1) DEFAULT NULL,
  `label` varchar(255) DEFAULT NULL,
  `tigem_ex_ref` varchar(40) NOT NULL COMMENT 'old tigem trap name',
  `trapcluster_id` int(10) unsigned DEFAULT NULL,
  `user` varchar(255) DEFAULT NULL,
  `note` mediumtext,
  PRIMARY KEY (`trapclusteradditional_id`),
  UNIQUE KEY `pair` (`trapcluster_id`,`label`),
  KEY `trap_id` (`trapcluster_id`),
  KEY `processing` (`processing`),
  KEY `quality` (`quality`),
  KEY `verified` (`verified`),
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=9 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapclusterblock`
--

CREATE TABLE IF NOT EXISTS `trapclusterblock` (
  `trapclusterblock_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `trapclustermap_id` int(10) NOT NULL DEFAULT '0',
  `start` int(10) NOT NULL DEFAULT '0',
  `end` int(10) NOT NULL DEFAULT '0',
  `strand` enum('1','-1') NOT NULL,
  `sequence` mediumtext NOT NULL,
  PRIMARY KEY (`trapclusterblock_id`),
  KEY `t` (`trapclustermap_id`),
  KEY `c` (`start`,`end`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=6 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapclusterblock_annotation`
--

CREATE TABLE IF NOT EXISTS `trapclusterblock_annotation` (
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
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapclustercheck`
--

CREATE TABLE IF NOT EXISTS `trapclustercheck` (
  `checkid` int(11) NOT NULL AUTO_INCREMENT,
  `trapcluster_id` int(10) unsigned DEFAULT NULL,
  `mapped` tinyint(1) NOT NULL DEFAULT '0',
  `annotated` int(11) NOT NULL,
  `checked` tinyint(1) NOT NULL DEFAULT '0',
  `mapping_type` tinyint(1) DEFAULT NULL,
  `splk` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`checkid`),
  UNIQUE KEY `trap_id` (`trapcluster_id`),
  UNIQUE KEY `trap_id_2` (`trapcluster_id`,`mapped`,`checked`,`mapping_type`,`splk`),
  KEY `annotated` (`annotated`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapclustermap`
--

CREATE TABLE IF NOT EXISTS `trapclustermap` (
  `trapclustermap_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `trapcluster_id` int(10) unsigned NOT NULL DEFAULT '0',
  `hit_id` varchar(40) NOT NULL DEFAULT '0',
  `hit_db` varchar(40) NOT NULL DEFAULT '',
  `start` int(10) NOT NULL DEFAULT '0',
  `end` int(10) NOT NULL DEFAULT '0',
  `strand` enum('1','-1') NOT NULL,
  PRIMARY KEY (`trapclustermap_id`),
  UNIQUE KEY `clusthit` (`trapcluster_id`,`start`,`end`,`strand`),
  KEY `hit_id` (`hit_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=6 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapclustermap_region`
--

CREATE TABLE IF NOT EXISTS `trapclustermap_region` (
  `trapclustermap_region_id` int(11) NOT NULL AUTO_INCREMENT,
  `region_id` int(10) unsigned NOT NULL DEFAULT '0',
  `type` varchar(40) NOT NULL,
  `annotation_ambiguous` tinyint(1) NOT NULL DEFAULT '1',
  `total_number_exons` int(11) NOT NULL,
  `number_trapclusterblocks` int(10) unsigned DEFAULT NULL,
  `number_annotated_trapclusterblocks` int(10) unsigned DEFAULT NULL,
  `overlap` int(10) NOT NULL,
  `trapclustermap_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`trapclustermap_region_id`),
  UNIQUE KEY `type` (`type`,`annotation_ambiguous`,`number_trapclusterblocks`,`number_annotated_trapclusterblocks`,`overlap`,`trapclustermap_id`,`region_id`),
  KEY `region_id_idx` (`region_id`),
  KEY `trapmap_id_idx` (`trapclustermap_id`),
  KEY `total_number_exons` (`total_number_exons`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=422 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapmap`
--

CREATE TABLE IF NOT EXISTS `trapmap` (
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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=592203 ;

-- --------------------------------------------------------

--
-- Table structure for table `trapmap_region`
--

CREATE TABLE IF NOT EXISTS `trapmap_region` (
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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=335416 ;

-- --------------------------------------------------------

--
-- Table structure for table `traprepeat`
--

CREATE TABLE IF NOT EXISTS `traprepeat` (
  `traprepeat_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `trap_id` int(10) unsigned DEFAULT NULL,
  `start` int(10) NOT NULL DEFAULT '0',
  `end` int(10) NOT NULL DEFAULT '0',
  `score` float DEFAULT NULL,
  `tag` varchar(25) DEFAULT NULL,
  PRIMARY KEY (`traprepeat_id`),
  UNIQUE KEY `trap_id` (`trap_id`,`start`,`end`,`score`,`tag`),
  KEY `trap_id_idxfk` (`trap_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `trap_error`
--

CREATE TABLE IF NOT EXISTS `trap_error` (
  `trap_error_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `logic_name_id` int(10) unsigned NOT NULL,
  `table_id` int(10) unsigned NOT NULL DEFAULT '0',
  `note` text NOT NULL,
  `table_name` varchar(255) DEFAULT NULL,
  `timestamp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`trap_error_id`),
  UNIQUE KEY `table_id` (`table_id`,`table_name`,`logic_name_id`),
  KEY `table_id_idxfk` (`table_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=2 ;

-- --------------------------------------------------------

--
-- Table structure for table `trap_maxicluster`
--

CREATE TABLE IF NOT EXISTS `trap_maxicluster` (
  `trap_maxicluster_id` int(11) NOT NULL AUTO_INCREMENT,
  `trap_id` int(11) NOT NULL,
  `maxicluster_id` int(11) NOT NULL,
  `trapmap_id` int(11) NOT NULL,
  PRIMARY KEY (`trap_maxicluster_id`),
  KEY `trap_id` (`trap_id`,`trapmap_id`),
  KEY `maxicluter_id` (`maxicluster_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=106 ;

-- --------------------------------------------------------

--
-- Table structure for table `trap_trapcluster`
--

CREATE TABLE IF NOT EXISTS `trap_trapcluster` (
  `trap_trapcluster_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `trap_id` int(10) DEFAULT NULL,
  `trapcluster_id` int(10) DEFAULT NULL,
  `trapmap_id` int(10) NOT NULL DEFAULT '0',
  PRIMARY KEY (`trap_trapcluster_id`),
  UNIQUE KEY `tum` (`trap_id`,`trapcluster_id`,`trapmap_id`),
  KEY `u` (`trap_id`),
  KEY `t` (`trapcluster_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=13 ;

-- --------------------------------------------------------

--
-- Table structure for table `trap_unitrap`
--

CREATE TABLE IF NOT EXISTS `trap_unitrap` (
  `trap_unitrap_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `trap_id` int(10) unsigned DEFAULT NULL,
  `trapblock_id` int(10) NOT NULL DEFAULT '0',
  `unitrap_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`trap_unitrap_id`),
  UNIQUE KEY `tum` (`trap_id`,`unitrap_id`,`trapblock_id`),
  KEY `trap_id_idxfk` (`trap_id`),
  KEY `unitrap_id_idxfk` (`unitrap_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=87655 ;

-- --------------------------------------------------------

--
-- Table structure for table `unitrap`
--

CREATE TABLE IF NOT EXISTS `unitrap` (
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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=24476 ;

-- --------------------------------------------------------

--
-- Table structure for table `unitrap_history`
--

CREATE TABLE IF NOT EXISTS `unitrap_history` (
  `history_id` int(11) NOT NULL AUTO_INCREMENT,
  `old_unitrap_id` varchar(255) NOT NULL,
  `new_unitrap_id` varchar(255) NOT NULL,
  PRIMARY KEY (`history_id`),
  KEY `old_unitrap_id` (`old_unitrap_id`,`new_unitrap_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 COMMENT='historic link between unitrap accession on the web site and ' AUTO_INCREMENT=18843 ;

-- --------------------------------------------------------

--
-- Table structure for table `xref`
--

CREATE TABLE IF NOT EXISTS `xref` (
  `xref_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `region_id` int(10) unsigned NOT NULL,
  `dbname` varchar(40) DEFAULT NULL,
  `accession` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`xref_id`),
  UNIQUE KEY `region_id` (`region_id`,`dbname`,`accession`),
  KEY `region_id_idxfk` (`region_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;
