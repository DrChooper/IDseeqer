###################################################################################
#
# The IDSeeqer.py script can be run directly from terminal like below:
#
# Run directly with parameters
# python IDSeeqer.py "['localhost', 'crop_pal2v2', 'user', 'password', 'test_idseeqer_temp', 2]"
#
# Run without parameters to display GUI
# python IDSeeqer.py
#
##################################################################################

from sqlalchemy import create_engine, text, select, MetaData

import ast
import time
import sys

import settings

import Gramene_Retrieval as gramene
import UniProt_Retrieval as uniprot
import NCBI_PROT_Retrieval as ncbiprot
import RAP_Retrieval as rap
import NCBI_NUCL_Retrieval as ncbinucl
import UniParc_Retrieval as uniparc
#import Nuccore_Retrieval as nuccore

import tkinter as tk
import tkinter.font as font

delay = 2

def check_delay(delay) :
	
	if(not isinstance(delay, int)) :
		return False
	elif (delay < 0) :
		return False
	else :
		return True

def Correct_Inputs(root, entries) :

	global args
	global proceed
	global delay
	
	host = entries['Host'].get()
	database_name = entries['Database'].get()
	user_name = entries['User'].get()
	user_password = entries['Password'].get()
	table_name = entries['Table'].get()
	
	
	
	args = [host, database_name, user_name, user_password, table_name]
	
	try :
		delay = int(entries['delay'].get())
		if (not check_delay(delay)) :
			message.config(text = "Please, enter a valid number for retrieval delay !", fg = 'red')
			return -1
	except Exception :
		message.config(text = "Please, enter a valid number for retrieval delay !", fg = 'red')
		return -1
	
	connected = settings.init(args)
	
	if(connected != 0) :
		message.config(text = "Failed to connect to database. Please, check the information provided !", fg = 'red')
		return -1
	else :		
		#message.config(text = "Successful connection to database. Program Running !", fg = 'green')
		proceed = True
	
	root.destroy()
	
	
def Cancel(root) :
	global proceed
	proceed = False
	root.destroy()
	


def makeform(root, fields):
	global message
	entries = {}
	counter = 0
	for field in fields:
		counter += 1
		row = tk.Frame(root)
		lab = tk.Label(row, width=13, text=field+": ", anchor='w')
		ent = tk.Entry(row, width=50)
		
		if (counter == 1) :
			ent.insert(0, "127.0.0.1")
		elif (counter == 2) :
			ent.insert(0, "crop_pal2v2")
		elif (counter == 3) :
			ent.insert(0, "root")
		elif (counter == 4) :
			ent.insert(0, "")
			ent.config(show = "") # Use show = "*" to display * in password field
		elif (counter == 5) :
			ent.insert(0, "test_idseeqer_temp")
			
		row.pack(side=tk.TOP, padx=5, pady=10)
		#row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=15)
		lab.pack(side=tk.LEFT)
		ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
		entries[field] = ent
	
	row = tk.Frame(root)
	lab = tk.Label(row, width=13, text="Retrieval Delay:", anchor = 'w')
	ent = tk.Entry(row, width=50)
	ent.insert(0, 2)
	row.pack(side=tk.TOP, padx=5, pady=15)
	lab.pack(side=tk.LEFT)
	ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
	entries['delay'] = ent
	
	row = tk.Frame(root)	
	arial10b= font.Font(family='Arial', size=10, weight='bold')
	message = tk.Label(row, width=400, height=4, text="", anchor='w', font = arial10b, borderwidth=2, relief="groove")
	row.pack(side=tk.TOP, padx=5, pady=15)
	message.pack(side=tk.LEFT, padx=15, pady=15)
	
	return entries
	
def Get_Inputs() :
	global root
	fields = ('Host', 'Database', 'User', 'Password', 'Table')
	
	root = tk.Tk()
	
	width = 500
	height = 470
	
	# get screen width and height
	screen_width = root.winfo_screenwidth()
	screen_height = root.winfo_screenheight()

	# calculate position x and y coordinates
	x = (screen_width/2) - (width/2)
	y = (screen_height/2) - (height/2)

	root.geometry('%dx%d+%d+%d' % (width, height, x, y)) # Windows width, height and position on screen.
	root.title('IDSeeqer - [Input Parameters]')
	root.resizable(False, False) # Window not resizable in both direction
	
	ents = makeform(root, fields)
	
	arial10b = font.Font(family='Arial', size=10, weight='bold')
	
	b1 = tk.Button(root, text='Process', command=(lambda e=ents: Correct_Inputs(root, e)), height = 2, width = 10, font = arial10b)
	b1.pack(side=tk.LEFT, padx=15, pady=15)
	
	b2 = tk.Button(root, text='Cancel', command=(lambda : Cancel(root)), height = 2, width = 10, font = arial10b)
	b2.pack(side=tk.RIGHT, padx=15, pady=15)
	
	root.bind('<Return>', (lambda e=ents, b1=b1: b1.invoke()))
	
	#root.protocol("WM_DELETE_WINDOW", Cancel(root)) # Catch the close (X) button
	root.mainloop()
	


def Run_Function( fun, *args ):
	return fun( *args )



if __name__ == "__main__" : 

	proceed = False
	
	used_gui = False
	
	
	
	function_names = [uniprot.search_id, ncbiprot.search_id, 
						gramene.search_id, uniparc.search_id, 
						rap.search_id, ncbinucl.search_id]
						
	display_names = {uniprot.search_id : 'UniProt', ncbiprot.search_id : 'NCBI_PROT',
				   gramene.search_id : 'Gramene', uniparc.search_id : 'UniParc', 
				   rap.search_id : 'RAP', ncbinucl.search_id : 'NCBI_NUCL'}
				   
	args = sys.argv
	
	try :
		args = ast.literal_eval(sys.argv[1])
		delay = int(args[5])
		
		if( not(type(args) is list) or (len(args) < 6) or (not check_delay(delay) ) ) :
			#print("Please, provide host, database, user, password and table as input")
			used_gui = True
			Get_Inputs()
						
			#sys.exit()
		else :
			delay = args[5]
			del args[5] 
			proceed = True
	except Exception as e: #(IndexError, SyntaxError) :
		used_gui = True
		Get_Inputs()


	if (proceed == False) :
		sys.exit()

	#### INITIALIZE SETTINGS.PY WHEN GUI IS NOT INVOLVED.
	
	if (not used_gui) : 
		connection_status = settings.init(args)
		if (connection_status != 0) : # Error connecting to database
			print("\n> CAN NOT CONNECT TO DATABASE : ", connection_status )
			sys.exit()
	
	
	# if (connection_status != 0) : # Error connecting to database
		# print("\n> CAN NOT CONNECT TO DATABASE : ", connection_status )
		# sys.exit()
	

	print("\n> SUCCESSFUL CONNECTION TO DATABASE.")
	
	###Define the nucleotide query that reoccurs after every target 
	q0 = select([settings.gfp_input_2019.c.id_paper1.distinct()])
	q1 = q0.where(settings.gfp_input_2019.c.id_paper1!=None)
	q9 = q1.where(settings.gfp_input_2019.c.seq_prot==None)
	q10 = q9.where(settings.gfp_input_2019.c.seq_nucl==None)
	#q20 = q10.where(settings.gfp_input_2019.c.errcol>1 )
	qp = q10.limit(settings.gramene_params.chunk)

	#I have to think about whether the query needs to change or how the user may want to interact with it.

	res=[r.id_paper1 for r in settings.engine.execute(qp).fetchall()]
	
	
	function_params = {uniprot.search_id : settings.uniprot_params, ncbiprot.search_id : settings.ncbiprot_params,
				   gramene.search_id : settings.gramene_params, uniparc.search_id : settings.uniparc_params, 
				   rap.search_id : settings.rap_params, ncbinucl.search_id : settings.ncbinucl_params}
		
	

	max_fail = settings.gramene_params.fail_num
	
	counter = 0
	total = len(res)
	
	# Loop over all ids and search
	for id_paper1 in res:
		counter += 1
		
		#if(counter > 5) :  # The first 5 records, for testing.
			#break
		
		
		print("\n>>> Searching Item {} out of {}.".format(counter, total))
	
		for func in function_names :
			
			print("> Searching '{}' ...".format(display_names[func].upper()))
			if( Run_Function(func, id_paper1, function_params[func]) ) :
				
				# func is the last function in list where search was successful
				# Remove the list
				function_names.remove(func)
				# Now insert it at the beginning of the list.
				function_names.insert(0, func)
				
				time.sleep(delay) # Wait for delay seconds after successful retrieval
				break 
			else : 
				max_fail -= 1
		
		
		
		if(max_fail <= 0) :
			print("\n > Reached ", settings.gramene_params.fail_num, " fails.\nExitting !\n")
			break
		
			
				
	print(' > Program Completed ')