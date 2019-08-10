from sqlalchemy import create_engine, text, MetaData


def init() :
	global gfp_input_2019	
	global engine
	#global id_paper1
	#Set up connection to database

	#engine = create_engine('mysql+pymysql://root@localhost/crop_pal2v2')
	engine = create_engine('mysql+pymysql://root:PASSWORD@localhost/DATABASE')#Enter your own password and database
	# test connection
	engine.execute(text("show variables like \"ver%\"")).fetchall()

	#extract the table information
	
	meta = MetaData(bind=engine)
	meta.reflect()
	meta.tables.keys()


	#define the table
	gfp_input_2019 = meta.tables['test_idseeqer_temp']