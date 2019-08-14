from io import StringIO
import requests
from Bio import SeqIO

import time

import settings

#find protein sequence from ID
def getncbi(pid, session=None):
    if isinstance(pid, (list, set)):
        pid = ",".join(str(s) for s in pid)
    resp = (session or requests).get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params=dict(id=pid, db="protein", rettype="gb", retmode="text"),
    )

    return SeqIO.read(StringIO(resp.text), format="gb")

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
	
def search_id (id_paper1, ncbi_params) :
	try:
		rec = getncbi(id_paper1)
		found = updatedb_prot(id_paper1, rec, ncbi_params.seq_source)
		print(id_paper1, 'updated', found, 'rows')
		#time.sleep(ncbi_params.retrieval_delay)
		return 0
	except Exception as e:
		print('failed for:',id_paper1, e)
		updateerr(id_paper1, ncbi_params.ncbi_prot_err)
		return -1

