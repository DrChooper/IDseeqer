from io import StringIO
import requests
from Bio import SeqIO

import time

import settings



#gramene
def getgramene(pid, session=None):
    resp = (session or requests).get(
        "http://rest.ensembl.org/sequence/id/%s?content-type=text/x-fasta;type=protein"% pid)
    return SeqIO.read(StringIO(resp.text), format = "fasta")
	
#update the protein sequences in the databae
def updatedb_prot(id_paper1, rec, seq_source):
	u = settings.gfp_input_2019.update()
	u = u.values({settings.gfp_input_2019.c.seq_prot: str(rec.seq).upper()})
	u = u.values({settings.gfp_input_2019.c.seq_source: seq_source})
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

def updateresult(id_paper1, source):
	u = settings.gfp_input_2019.update()
	u = u.values({settings.gfp_input_2019.c.id_type: source})
	u = u.where(settings.gfp_input_2019.c.id_paper1 == id_paper1)

	proxy = settings.engine.execute(u)
	return proxy.rowcount


def search_id (id_paper1, gramene_params) :
	try:
		rec = getgramene(id_paper1)
		found = updatedb_prot(id_paper1, rec, gramene_params.seq_source)
		print(id_paper1, 'updated', found, 'rows')
		#updateresult(id_paper1,"ensemble")
		#time.sleep(gramene_params.retrieval_delay)
		return True
	except Exception as e:
		print('failed for:',id_paper1, e)
		updateerr(id_paper1, gramene_params.ensembl_err)
		
		return False

		

	


