from io import StringIO
import requests
from Bio import SeqIO

import time

import settings


#UNIprot
def getuniprot(pid, session=None):
    resp = (session or requests).get("https://www.uniprot.org/uniprot/%s.xml" % pid)
    return SeqIO.read(StringIO(resp.text), format="uniprot-xml")
	
#update the protein sequences in the databae
def updatedb_prot(id_paper1, rec):
    u = settings.gfp_input_2019.update()
    u = u.values({settings.gfp_input_2019.c.seq_prot: str(rec.seq).upper()})
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
	


def search_id (id_paper1, uniprot_params) :

	try:
		rec = getuniprot(id_paper1)
		found = updatedb_prot(id_paper1, rec)
		print(id_paper1, 'updated', found, 'rows')
		#time.sleep(uniprot_params.retrieval_delay)
		return 0
	except Exception as e:
		print('failed for:',id_paper1, e)
		updateerr(id_paper1, uniprot_params.uniprot_err)
		return -1
		
