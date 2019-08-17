from io import StringIO
import requests
from Bio import SeqIO

import requests
import re
from lxml.html import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from sqlalchemy import create_engine, text
from sqlalchemy import select, MetaData

import time

import settings



#find protein sequence from ID
def getuniparc(name, session=None):
    resp = requests.get(
        "https://www.uniprot.org/uniparc/",
         params=dict(query=name)
    )
    r = parse(StringIO(resp.text))
    r = r.getroot()
    t = list(r.xpath('.//table[@id="results"]//td[@class="entryID"]/a/text()'))
    rec = t[0]
   
    return rec 

def getuniparc_d(name, session=None):
    resp = requests.get(
        "https://www.uniprot.org/uniparc/{name}".format(name=name)    
    )
    r = parse(StringIO(resp.text))
    r = r.getroot()
    t = list(r.xpath('.//table[@id="results"]//td[@class="entryID"]/a/text()'))
    rec = t[0]
   
    return rec 



#update the protein sequences in the databae
def updatedb_uniparc1(id_paper1, rec, seq_source):
	u = settings.gfp_input_2019.update()
	u = u.values({settings.gfp_input_2019.c.id_paper2: rec.upper()})
	u = u.values({settings.gfp_input_2019.c.seq_source: seq_source})
	u = u.where(settings.gfp_input_2019.c.id_paper1 == id_paper1)
	u = u.where(settings.gfp_input_2019.c.seq_prot == None)

	proxy = settings.engine.execute(u)
	return proxy.rowcount

def updatedb_uniparc_d(id_paper1, rec):
    u = settings.gfp_input_2019.update()
    u = u.values({settings.gfp_input_2019.c.id_paper1: rec.upper()})
    u = u.where(settings.gfp_input_2019.c.id_paper1 == id_paper1)
    u = u.where(settings.gfp_input_2019.c.seq_prot == None)

    proxy = settings.engine.execute(u)
    return proxy.rowcount

def updateerr(id_paper1, num):
    u = settings.gfp_input_2019.update()
    u = u.values({settings.gfp_input_2019.c.errcol: num})
    u = u.where(settings.gfp_input_2019.c.id_paper1 == id_paper1)

    proxy = settings.engine.execute(u)
    return proxy.rowcount
	


#find protein sequence from ID
def getuniparc2(name, session=None):
    resp = requests.get(
        "https://www.uniprot.org/uniparc/{name}".format(name=name)
    )
    r = parse(StringIO(resp.text))
    r = r.getroot()
    t = list(r.xpath('.//table[@id="sequence"]//pre[@class="sequence"]/text()'))
    S=re.compile(r'[\s\n0-9]+')
    seq=S.sub('',''.join(t))
    rec=SeqRecord(seq)
    return rec 


#UPDATE and Error
def updatedb_uniparc2(id_paper2, rec):
    u = settings.gfp_input_2019.update()
    u = u.values({settings.gfp_input_2019.c.seq_prot: str(rec.seq).upper()})
    u = u.where(settings.gfp_input_2019.c.id_paper2 == id_paper2)
    u = u.where(settings.gfp_input_2019.c.seq_prot == None)

    proxy = settings.engine.execute(u)
    return proxy.rowcount

def updateerr(id_paper2, num):
    u = settings.gfp_input_2019.update()
    u = u.values({settings.gfp_input_2019.c.errcol: num})
    u = u.where(settings.gfp_input_2019.c.id_paper2 == id_paper2)

    proxy = settings.engine.execute(u)
    return proxy.rowcount


def inner_search (id_paper2, uniparc_params) :
	q0 = select([settings.gfp_input_2019.c.id_paper2.distinct()])
	q1 = q0.where(settings.gfp_input_2019.c.id_paper2 == id_paper2)
	#q1 = q0.where(settings.gfp_input_2019.c.id_paper2.like("UPI%"))
	#q2 = q10.where(settings.gfp_input_2019.c.id_paper2!=None)
	qs = q1.where(settings.gfp_input_2019.c.seq_prot==None)
	#q22 = q21.where(settings.gfp_input_2019.c.seq_nucl==None)
	#q20 = q22.where(settings.gfp_input_2019.c.errcol<5) #0,1 for uniprot 0,1,2 for NCBI 0,1,2,3 Gramene
	#qs = q22.limit(uniparc_params.chunk)

	
	res=[r.id_paper2 for r in settings.engine.execute(qs).fetchall()]
	tot = len(res)
	
	temp_count = 0

	for id_paper2 in res:

		temp_count += 1
		
		print("\n>>>>>> In UNIPARC: Searching Item {} out of {}.".format(temp_count, tot))
		try:
			rec = getuniparc2(id_paper2)
			found = updatedb_uniparc2(id_paper2, rec)
			print(' >>>>>> ', id_paper2, 'updated', found, 'rows')
			#time.sleep(retrieval_delay)
		except Exception as e:
			print(' >>>>>> ', 'failed for:',id_paper2, e)
			updateerr(id_paper1, uniparc_err)

	return 0




#Searching Uniparc ids
def search_id (id_paper1, uniparc_params) :
	try:
		rec = getuniparc(id_paper1)
		found = updatedb_uniparc1(id_paper1, rec, uniparc_params.seq_source)
		print(id_paper1, 'updated', found, 'rows')
		inner_search(rec.upper(), uniparc_params)
		#time.sleep(retrieval_delay)
		return True
	except Exception as e:
		print('failed for:',id_paper1, e)
		updateerr(id_paper1, uniparc_params.uniparc_err)

		return False







