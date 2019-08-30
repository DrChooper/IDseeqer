from io import StringIO
import requests
from Bio import SeqIO

from lxml.html import parse

from sqlalchemy import create_engine, text
from sqlalchemy import select, MetaData

import sys
import os
import platform

from ftplib import FTP
import subprocess

import settings

import Blast_Primer
import Blast_Prot
import Blast_Nucl


seq_type = ["cdna", "pep"]

def my_download(species_name, species_name2, taxaid, file_name, release, seqtype) :
	
	if(platform.system() == "Windows") :
		mycurl = "curl"
		mygzcat = "zcat"
	else :
		mycurl = "curl"
		mygzcat = "gzcat"
		
	
	if(os.path.isfile(file_name) == False) : # file does not it, download it
		print("Please, wait while {} is downloaded.\n".format(file_name))
		d_string = f'ftp://ftp.ensemblgenomes.org/pub/release-{settings.db_release}/plants/fasta/{species_name}/{seqtype}/{file_name}'
		
		res = subprocess.run([mycurl, "-o", file_name, d_string], universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		res = res.returncode
	else :
		print("File '{}' already exists !\n".format(file_name))
		res = 0
		
	if(res == 0) :
		print("Building BLAST library for '{}'.\n".format(file_name))
		#prefix = species_name.split('_', 1)
		#prefix = prefix[0][0] + prefix[1][0]
		prefix = str(taxaid).strip()
		
		d_string = mygzcat + " " + file_name
		
		#res = subprocess.run([mygzcat, file_name], universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		#res = res.returncode
		
		if(res == 0) :
			if(seqtype == seq_type[0]) : #CDNA
				d_string += " |makeblastdb -in - -out blastdb_" + prefix + "_nucl -input_type=fasta -dbtype nucl -title " + species_name2
				
			
			elif(seqtype == seq_type[1]) : #PEP
				d_string += " |makeblastdb -in - -out bastdb_" + prefix + "_prot -input_type=fasta -dbtype prot -title " + species_name2
			
			print(d_string)
			res2 = os.system(d_string)
			#print("\nres = {}.\n".format(res))
			#res = subprocess.run(["makeblastdb", "-in", "-", "-out", "blastdb"_+prefix+"_prot", "-input_type=fasta", "-dbtype prot", "-title "+species_name2])
			
		
def Retrieval(taxa_chunk) :	

	q = select([settings.species_list.c.taxaid, settings.species_list.c.species_name, settings.species_list.c.species_name2])
	q = q.limit(taxa_chunk)

	species=[[r.species_name, r.species_name2, r.taxaid] for r in settings.engine.execute(q).fetchall()]

	cdna_fasta_files = []
	pep_fasta_files = []

	try :
		ftp = FTP('ftp.ensemblgenomes.org')
		login_res = ftp.login()
	except Exception as e :
		print("CAN NOT LOGIN AT ftp.ensemblgenomes.org !\n", e)
		sys.exit()

	for species_name, species_name2, taxaid in species:
		print(species_name)
		for i in range(len(seq_type)) :
			URL_list = f'/pub/release-{settings.db_release}/plants/fasta/{species_name}/{seq_type[i]}/' 
			ftp.cwd(URL_list)
			l=ftp.nlst()
			next_file = [f for f in l if f.endswith('.all.fa.gz')][0]
			
			if(i == 0) :
				cdna_fasta_files.append([species_name, species_name2, taxaid, next_file])
			elif(i == 1) :
				pep_fasta_files.append([species_name, species_name2, taxaid, next_file])

	ftp.close()

	# Download CDNA Files and Build Blast
	for species_name, species_name2, taxaid, file_name in cdna_fasta_files :
		my_download(species_name, species_name2, taxaid, file_name, settings.db_release, seq_type[0])
		
	# Download PEP Files and Build Blast
	for species_name, species_name2, taxaid, file_name in pep_fasta_files :
		my_download(species_name, species_name2, taxaid, file_name, settings.db_release, seq_type[1])	
		
	
	total = len(species)
	counter = 0
	# Now Retrieve info from BLAST
	for species_name, species_name2, taxaid in species :
		counter += 1
		print("\n>>> BLAST is searching '{}', item {} out of {}.".format(taxaid, counter, total))
		
		print("\n > NOW SEARCHING BLAST_PRIMER FOR TAXA_ID: '{}'".format(taxaid))
		Blast_Primer.Search_Blast(taxaid)
		
		print("\n > NOW SEARCHING BLAST_PROT FOR TAXA_ID: '{}'".format(taxaid))
		Blast_Prot.Search_Blast(taxaid)
		
		print("\n > NOW SEARCHING BLAST_NUCL FOR TAXA_ID: '{}'".format(taxaid))
		Blast_Nucl.Search_Blast(taxaid)
		

    
