# ************************************************************
# Sequel Pro SQL dump
# Version 4541
#
# http://www.sequelpro.com/
# https://github.com/sequelpro/sequelpro
#
# Host: 127.0.0.1 (MySQL 5.6.41)
# Database: crop_pal2v2
# Generation Time: 2019-08-06 06:26:56 +0000
# ************************************************************


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;


# Dump of table taxa
# ------------------------------------------------------------

DROP TABLE IF EXISTS `taxa`;

CREATE TABLE `taxa` (
  `taxaid` int(11) unsigned NOT NULL,
  `species_name` varchar(255) DEFAULT NULL,
  `species_name2` varchar(255) DEFAULT NULL,
  `species_name3` varchar(255) DEFAULT NULL,
  `species` varchar(255) DEFAULT NULL,
  `common_name` varchar(255) DEFAULT NULL,
  `comments` text
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

LOCK TABLES `taxa` WRITE;
/*!40000 ALTER TABLE `taxa` DISABLE KEYS */;

INSERT INTO `taxa` (`taxaid`, `species_name`, `species_name2`, `species_name3`, `species`, `common_name`, `comments`)
VALUES
	(3702,'arabidopsis_thaliana','athaliana','a_thaliana','arabidopsis thaliana','thale cress','dicot'),
	(4530,'oryza_sativa','osativa','o_sativa','oryza sativa','rice','monocot'),
	(4577,'zea_mays','zmays','z_mays','zea mays','maize','monocot'),
	(4081,'solanum_lycopersicum','slycopersicum','s_lycopersicum','solanum lycopersicum','tomato','dicot'),
	(4113,'solanum_tuberosum','stuberosum','s_tuberosum','solanum tuberosum','potato','dicot'),
	(29760,'vitis_vinifera','vvinifera','v_vinifera','vitis vinifera','wine grape','dicot'),
	(3708,'brassica_napus','bnapus','b_napus','brassica napus','canola','dicot'),
	(3847,'glycine_max','gmax','g_max','glycine max','soybean','dicot'),
	(4558,'sorghum_bicolor','sbicolor','s_bicolor','sorghum bicolor','sorghum','monocot'),
	(4641,'musa_acuminata','macuminata','m_acuminata','musa acuminata','banana','monocot'),
	(3711,'brassica_rapa','brapa','b_rapa','brassica rapa','field mustard','dicot'),
	(4513,'hordeum_vulgare','hvulgare','h_vulgare','hordeum vulgare','barley','monocot'),
	(4565,'triticum_aestivum','taestivum','t_aestivum','triticum aestivum','bread wheat','monocot');

/*!40000 ALTER TABLE `taxa` ENABLE KEYS */;
UNLOCK TABLES;



/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;
/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
