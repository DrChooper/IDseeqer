from sqlalchemy import create_engine, text, select, MetaData

import time
import sys

import settings
import Gramene_Retrieval as gramene
import UniProt_Retrieval as uniprot
import NCBI_Retrieval as ncbi


class Parameters:
	chunk = 100000 #how many ids in the run limit
	retrieval_delay = 2 #time break after successful retrieval in seconds
	error_delay =2 #time break after error in seconds
	max_fail=  200
	fail_num = 200
	ncbi_nuc_err = 1
	uniprot_err = 2
	ncbi_prot_err = 3
	ensembl_err = 4
	err_status = 2 #where to start sending the sequence?
	


if __name__ == "__main__" : 

	connection_status = settings.init()
	if (connection_status != 0) : # Error connecting to database
		print("\n> CAN NOT CONNECT TO DATABASE : ", connection_status )
		sys.exit()
	

	# Define Run Variables
	### define errors
	gramene_params = Parameters()

	uniprot_params = Parameters()
	uniprot_params.chunk = 5000 #how many ids in the run limit
	uniprot_params.retrieval_delay = 3 #time break after successful retrieval in seconds
	uniprot_params.error_delay =3 #time break after error in seconds
	uniprot_params.max_fail=  50 #some number to play around with for what happens if there are too many fails? e.g. move onto the next script
	uniprot_params.fail_num = 50 #as above
	#some system to save where the failed... or maybe change it to which one succeeded in the end.
	uniprot_params.ncbi_nuc_err = 1 
	uniprot_params.uniprot_err = 2
	uniprot_params.ncbi_prot_err = 3
	uniprot_params.ensembl_err = 4
	uniprot_params.err_status = 2 
	
	ncbi_params = Parameters()
	ncbi_params.chunk = 100000 #how many ids in the run limit
	ncbi_params.retrieval_delay = 2 #time break after successful retrieval in seconds
	ncbi_params.error_delay =2 #time break after error in seconds
	ncbi_params.max_fail=  200
	ncbi_params.fail_num = 200
	ncbi_params.ncbi_nuc_err = 1
	ncbi_params.uniprot_err = 2
	ncbi_params.ncbi_prot_err = 3
	ncbi_params.ensembl_err = 4
	ncbi_params.err_status = 2 #where to start sending the sequence?

	###Define the nucleotide query that reoccurs after every target 
	q0 = select([settings.gfp_input_2019.c.id_paper1.distinct()])
	q1 = q0.where(settings.gfp_input_2019.c.id_paper1!=None)
	q9 = q1.where(settings.gfp_input_2019.c.seq_prot==None)
	q10 = q9.where(settings.gfp_input_2019.c.seq_nucl==None)
	q20 = q10.where(settings.gfp_input_2019.c.errcol>1 )
	qp = q20.limit(gramene_params.chunk)

	#I have to think about whether the query needs to change or how the user may want to interact with it.

	res=[r.id_paper1 for r in settings.engine.execute(qp).fetchall()]
	

	#sys.exit()

	max_fail = gramene_params.fail_num
	
	counter = 0
	total = len(res)
	
	# Loop over all ids and search
	for id_paper1 in res:
		counter += 1
		
		print("\n>>> Searching Item {} out of {}.".format(counter, total))
	
		print("> Will now search Uniprot")
		solution = uniprot.search_id(id_paper1, uniprot_params)
				
		if (solution == 0) : # Found in Uniprot
			#print("Uniprot Returned ", solution)
			continue
		# Not Found in Uniprot
						
		max_fail -= 1
		if (max_fail <= 0):
			print("Maximum number (%d) of failures reached. "% max_fail)
			break
		#time.sleep(uniprot_params.error_delay)
		
		# Not found in uniprot, now search ncbi
		print(">> Not found in Uniprot, will now search NCBI")
		solution = ncbi.search_id(id_paper1, ncbi_params)
		
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
		print(">>> Not found in NCBI, will now search gramene")
		solution = gramene.search_id(id_paper1, gramene_params)
		
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
		
			
		
			
		
			
	print(' > Program Completed ')