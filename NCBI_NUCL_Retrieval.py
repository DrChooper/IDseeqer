from io import StringIO
import requests
from Bio import SeqIO

import time

import settings

#find nucleotide sequence from ID (if no entry in protein database)
def getncbi_nuc(pid, session=None):
    if isinstance(pid, (list, set)):
        pid = ",".join(str(s) for s in pid)
    resp = (session or requests).get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params=dict(id=pid, db="nucleotide", rettype="gb", retmode="text"),
    )

    return SeqIO.read(StringIO(resp.text), format="gb")


#UPDATE
#update the nucleotide sequences in the databae
def updatedb_nucl(id_paper1, rec, seq_source):
	u = settings.gfp_input_2019.update()
	u = u.values({settings.gfp_input_2019.c.seq_nucl: str(rec.seq).upper()})
	u = u.values({settings.gfp_input_2019.c.seq_source: seq_source})
	u = u.where(settings.gfp_input_2019.c.id_paper1 == id_paper1)
	u = u.where(settings.gfp_input_2019.c.seq_nucl == None)

	proxy = settings.engine.execute(u)
	return proxy.rowcount


def updateerr(id_paper1, num):
    u = settings.gfp_input_2019.update()
    u = u.values({settings.gfp_input_2019.c.errcol: num})
    u = u.where(settings.gfp_input_2019.c.id_paper1 == id_paper1)

    proxy = settings.engine.execute(u)
    return proxy.rowcount
	

#Searching the NCBI Nucleotide database
def search_id (id_paper1, ncbinucl_params) :
	try:
		rec = getncbi_nuc(id_paper1)
		found = updatedb_nucl(id_paper1, rec, ncbinucl_params.seq_source)
		print(id_paper1, 'updated', found, 'rows')
		#time.sleep(retrieval_delay)
		return True
	except Exception as e:
		print('failed for:',id_paper1, e)
		updateerr(id_paper1, ncbinucl_params.ncbi_nuc_err)

		return False
