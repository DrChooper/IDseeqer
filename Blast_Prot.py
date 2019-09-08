import os
import sys

from sqlalchemy import bindparam
from sqlalchemy import create_engine, text
from sqlalchemy import select, MetaData

import pandas as pd

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import settings

####DEFINE FUNCTION FOR UPDATING:

def updatedb_prot(recs, taxaid):
    u = settings.blast_input.update()
    u = u.values({settings.blast_input.c.locus_seq: bindparam('subject acc.ver')})
    u = u.where(settings.blast_input.c.id == bindparam('query acc.ver'))
    u = u.where(settings.blast_input.c.locus_seq==None)

    proxy = settings.engine.execute(u, *recs)
    return print("updated ",proxy.rowcount," entries for ", taxaid)
	

def Search_Blast(taxaid) : 

	#select species sequences - TA
	q0 = select([settings.blast_input.c.seq_prot,settings.blast_input.c.id])
	q2 = q0.where(settings.blast_input.c.locus_seq==None)
	q3 = q2.where(settings.blast_input.c.seq_prot!=None)
	q1 = q3.where(settings.blast_input.c.taxa==int(taxaid)) #ENTER TAXA FOR PLANT HERE
	qp = q1.limit(settings.blast_chunk)
	rec = list(settings.engine.execute(qp).fetchall())
	
	taxaid = str(taxaid).strip()
	
	records = [SeqRecord(id=str(r.id),seq=Seq(r.seq_prot))for r in rec]
	
	SeqIO.write(records,"query_prot_"+taxaid+".fasta","fasta")# ENTER THE species taxa here and in the blast query below (twice)!!
	#Enter the right blast database
	os.system("blastp -max_hsps 1 -max_target_seqs 5 -outfmt 6 -query query_prot_"+taxaid+".fasta -db blastdb_"+taxaid+"_prot -out prot_"+taxaid+".tsv")

	header='query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score'.split(',')
	header = [h.strip() for h in header]
	df = pd.read_csv("prot_"+taxaid+".tsv",header = 0, sep = "\t", names=header)#ENTER THE taxa here

	try :
		os.remove("query_prot_" + taxaid + ".fasta")
		os.remove("prot_" + taxaid + ".tsv")
	except Exception as e :
		pass

	if len(df) == 0:
		print("Could not find any hits in the blast database for the ",len(rec)," records.")

	else:

		dfs=df.sort_values(by = ["query acc.ver","bit score"], ascending = [True,False])
		idx = df.groupby("query acc.ver")['bit score'].transform(max) == df['bit score']
		dff=df[idx]
		d=dff[['query acc.ver', 'subject acc.ver']]
		recs=d.to_dict(orient='records')

		updatedb_prot(recs, taxaid)
		print("Updated taxa:  {}".format(taxaid)) #Final message with info on updated species taxa
