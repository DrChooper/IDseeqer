from sqlalchemy import create_engine, text, MetaData

# host = "localhost" # or "127.0.0.1"
# database_name = "crop_pal2v2"
# database_user = "root"
# database_password = ""
# table_name = "test_idseeqer_temp"

class Parameters:
	chunk = 100000 #how many ids in the run limit
	retrieval_delay = 2 #time break after successful retrieval in seconds
	error_delay = 2 #time break after error in seconds
	max_fail =  200
	fail_num = 200
	ncbi_nuc_err = 1
	uniprot_err = 2
	ncbi_prot_err = 3
	ensembl_err = 4
	err_status = 2 #where to start sending the sequence?
	rap_err = 20
	seq_source = "gramene"

def init(args) :
	global gfp_input_2019	
	global engine
	global gramene_params
	global uniprot_params
	global ncbiprot_params
	global ncbinucl_params
	global rap_params
	global uniparc_params
	
	
	host = args[0]
	database_name = args[1]
	database_user = args[2]
	database_password = args[3]
	table_name = args[4]
	
		
	#Set up connection to database
	try :
		conn_string = "mysql+pymysql://" + database_user + ":" + database_password + "@" + host + "/" + database_name
		engine = create_engine(conn_string)
		
		# test connection
		engine.execute(text("show variables like \"ver%\"")).fetchall()
		

		#extract the table information
		meta = MetaData(bind=engine)
		meta.reflect()
		meta.tables.keys()


		#define the table
		gfp_input_2019 = meta.tables[table_name]
	
	except Exception as e :
		return e
	
	
	
	### define errors
	gramene_params = Parameters()

	uniprot_params = Parameters()
	uniprot_params.chunk = 5000 #how many ids in the run limit
	uniprot_params.retrieval_delay = 3 #time break after successful retrieval in seconds
	uniprot_params.error_delay = 3 #time break after error in seconds
	uniprot_params.max_fail =  50 #some number to play around with for what happens if there are too many fails? e.g. move onto the next script
	uniprot_params.fail_num = 50 #as above
	#some system to save where the failed... or maybe change it to which one succeeded in the end.
	uniprot_params.ncbi_nuc_err = 1 
	uniprot_params.uniprot_err = 2
	uniprot_params.ncbi_prot_err = 3
	uniprot_params.ensembl_err = 4
	uniprot_params.err_status = 2 
	uniprot_params.seq_source = "uni_prot"
	
	ncbiprot_params = Parameters()
	ncbiprot_params.chunk = 100000 #how many ids in the run limit
	ncbiprot_params.retrieval_delay = 2 #time break after successful retrieval in seconds
	ncbiprot_params.error_delay = 2 #time break after error in seconds
	ncbiprot_params.max_fail =  200
	ncbiprot_params.fail_num = 200
	ncbiprot_params.ncbi_nuc_err = 1
	ncbiprot_params.uniprot_err = 2
	ncbiprot_params.ncbi_prot_err = 3
	ncbiprot_params.ensembl_err = 4
	ncbiprot_params.err_status = 2 #where to start sending the sequence?
	ncbiprot_params.seq_source = "ncbi_prot"
	
	rap_params = Parameters()
	rap_params.chunk = 1000 #how many ids in the run limit
	rap_params.retrieval_delay = 2 #time break after successful retrieval in seconds
	rap_params.error_delay = 10 #time break after error in seconds
	rap_params.max_fail =  500
	rap_params.fail_num = 500
	rap_params.ncbi_nuc_err = 10
	rap_params.uniprot_err = 20
	rap_params.ncbi_prot_err = 2
	rap_params.ensembl_err = 4
	rap_params.err_status = 10 #where to start sending the sequence?
	rap_params.rap_err = 20
	rap_params.seq_source = "rap"
	
	ncbinucl_params = Parameters()
	ncbinucl_params.chunk = 10000 #how many ids in the run limit
	ncbinucl_params.retrieval_delay = 2 #time break after successful retrieval in seconds
	ncbinucl_params.error_delay = 5 #time break after error in seconds
	ncbinucl_params.max_fail=  500
	ncbinucl_params.fail_num = 500
	ncbinucl_params.ncbi_nuc_err = 10
	ncbinucl_params.rap_err = 20
	ncbinucl_params.uniprot_err = 2
	ncbinucl_params.ncbi_prot_err = 3
	ncbinucl_params.ensembl_err = 4
	ncbinucl_params.uniparc_err = 5
	ncbinucl_params.err_status = 10 #where to start sending the sequence?
	ncbinucl_params.seq_source = "ncbi_nucl"
	
	uniparc_params = Parameters()
	uniparc_params.chunk = 10000 #how many ids in the run limit
	uniparc_params.retrieval_delay = 2 #time break after successful retrieval in seconds
	uniparc_params.error_delay = 5 #time break after error in seconds
	uniparc_params.max_fail=  500
	uniparc_params.fail_num = 500
	uniparc_params.ncbi_nuc_err = 10
	uniparc_params.rap_err = 20
	uniparc_params.uniprot_err = 2
	uniparc_params.ncbi_prot_err = 3
	uniparc_params.ensembl_err = 4
	uniparc_params.uniparc_err = 5
	uniparc_params.err_status = 10 #where to start sending the sequence?
	uniparc_params.seq_source = "uniparc"
	

	
	return 0
	
