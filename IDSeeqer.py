from sqlalchemy import create_engine, text, select, MetaData

import time
import sys

import settings
import Gramene_Retrieval as gramene
import UniProt_Retrieval as uniprot
import NCBI_PROT_Retrieval as ncbiprot
import RAP_Retrieval as rap


def Perform_Search( fun, *args ):
	return fun( *args )

p = 2
r = 5

#print(Perform_Search(settings.add, p , r))
#print(Perform_Search(settings.multiply, p, r))

#sys.exit()


if __name__ == "__main__" : 

	connection_status = settings.init()
	if (connection_status != 0) : # Error connecting to database
		print("\n> CAN NOT CONNECT TO DATABASE : ", connection_status )
		sys.exit()
	

	
	###Define the nucleotide query that reoccurs after every target 
	q0 = select([settings.gfp_input_2019.c.id_paper1.distinct()])
	q1 = q0.where(settings.gfp_input_2019.c.id_paper1!=None)
	q9 = q1.where(settings.gfp_input_2019.c.seq_prot==None)
	q10 = q9.where(settings.gfp_input_2019.c.seq_nucl==None)
	#q20 = q10.where(settings.gfp_input_2019.c.errcol>1 )
	qp = q10.limit(settings.gramene_params.chunk)

	#I have to think about whether the query needs to change or how the user may want to interact with it.

	res=[r.id_paper1 for r in settings.engine.execute(qp).fetchall()]
	

	#sys.exit()

	max_fail = settings.gramene_params.fail_num
	
	counter = 0
	total = len(res)
	
	# Loop over all ids and search
	for id_paper1 in res:
		counter += 1
		
		print("\n>>> Searching Item {} out of {}.".format(counter, total))
	
		print("> Will now search Uniprot")
		solution = uniprot.search_id(id_paper1, settings.uniprot_params)
				
		if (solution == 0) : # Found in Uniprot
			#print("Uniprot Returned ", solution)
			continue
		# Not Found in Uniprot
						
		max_fail -= 1
		if (max_fail <= 0):
			print("Maximum number (%d) of failures reached. "% max_fail)
			break
		#time.sleep(uniprot_params.error_delay)
		
		# Not found in uniprot, now search ncbiprot
		print(">> Not found in Uniprot, will now search NCBI_PROT")
		solution = ncbiprot.search_id(id_paper1, settings.ncbiprot_params)
		
		if (solution == 0) :
			#print("NCBI Returned ", solution)
			continue
		# Not Found in NCBI
		
		max_fail -= 1
		if (max_fail <= 0):
			print("Maximum number (%d) of failures reached. "% max_fail)
			break
		#time.sleep(ncbi_params.error_delay)
		
		# Not found in ncbi, now search gramene
		print(">>> Not found in NCBI_PROT, will now search gramene")
		solution = gramene.search_id(id_paper1, settings.gramene_params)
		
		if (solution == 0) :
			#print("Gramene Returned ", solution)
			continue
		# Not Found in Gramene	
		max_fail -= 1
		if (max_fail <= 0) :
			print("Maximum number (%d) of failures reached. "% max_fail)
			break
		#time.sleep(gramene_params.error_delay)
		
		# Not found in gramene, now search rap
		print(">>> Not found in gramene, will now search RAP")
		solution = rap.search_id(id_paper1, settings.rap_params)
		
		if (solution == 0) :
			#print("RAP Returned ", solution)
			continue
		# Not Found in RAP	
		max_fail -= 1
		if (max_fail <= 0) :
			print("Maximum number (%d) of failures reached. "% max_fail)
			break
		
			
				
	print(' > Program Completed ')