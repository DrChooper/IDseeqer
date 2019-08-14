from io import StringIO
import requests
from Bio import SeqIO

from lxml.html import parse

import time

import settings



    
#find protein sequence from ID
def getrap(name, session=None):
    resp = requests.get(
        "https://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1",
         params=dict(name=name)
    )
    r = parse(StringIO(resp.text))
    r = r.getroot()
    t = ''.join(r.xpath('.//pre/text()'))
    rec = list(SeqIO.parse(StringIO(t), format="fasta"))[0]
    return rec 
    
  
#UPDATE FUNCTION
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
	

def search_id (id_paper1, rap_params) :

	try:
		rec = getrap(id_paper1)
		found = updatedb_nucl(id_paper1, rec, rap_params.seq_source)
		print(id_paper1, 'updated', found, 'rows')
		#time.sleep(rap_params.retrieval_delay)
		return 0
	except Exception as e:
		print('failed for:',id_paper1, e)
		updateerr(id_paper1, rap_params.rap_err)

		return -1
		
