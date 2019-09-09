###################################################################################
#
# The IDSeeqer.py script can be run directly from terminal like below:
#
# Run directly with parameters
# python IDSeeqer.py "['localhost', 'crop_pal2v2', 'user', 'password', 'test_idseeqer_temp', 2, 'taxa', 40, 2, 1000, 'full path to blast folder']"
#
# Run without parameters to display GUI
# python IDSeeqer.py
#
##################################################################################

from sqlalchemy import create_engine, text, select, MetaData

import ast
import time
import sys
import os

import settings

import Gramene_Retrieval as gramene
import UniProt_Retrieval as uniprot
import NCBI_PROT_Retrieval as ncbiprot
import RAP_Retrieval as rap
import NCBI_NUCL_Retrieval as ncbinucl
import UniParc_Retrieval as uniparc
#import Nuccore_Retrieval as nuccore

import Blast_Retrieval

import tkinter as tk
from tkinter import ttk
import tkinter.font as font

from tkinter import filedialog
from tkinter import *

from timeit import default_timer as timer


delay = 2

def Enter_Hit():
	current_tab = tab_parent.select()
	if(current_tab == '.!notebook.!frame') :
		#print('Correct Tab')
		but_process.invoke()
	else :
		pass
	#sys.exit()
	
	

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
	global taxa_table
	global taxa_chunk
	global program_to_run
	
	host = entries['Host'].get()
	database_name = entries['Database'].get()
	user_name = entries['User'].get()
	user_password = entries['Password'].get()
	table_name = entries['Table'].get()
	
	taxa_table = entries['taxa_table'].get()
	db_release = entries['database_release'].get()
	taxa_chunk = int(entries['taxa_chunk'].get())
	blast_chunk = int(entries['blast_chunk'].get())
	
	program_to_run = entries['program'].get()
	
	
	blast_folder = entries['blast_folder'].get()
	
	blast_folder = blast_folder.strip()
	if (not os.path.isdir(blast_folder)) :
		message.config(text="Can not find the Blast folder you provided !", fg = 'red')
		return -1
	
	Save_Blast_Folder(blast_folder)
	
	args = [host, database_name, user_name, user_password, table_name, taxa_table, db_release, blast_chunk, blast_folder]
	
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
	
	return
	
	
def Cancel(root) :
	global proceed
	proceed = False
	root.destroy()
	
def browse_button():
	# Allow user to select a directory and store it in global var
	
	filename = filedialog.askdirectory(initialdir=init_blast_folder, title='Please select the directory of the Blast database.')
	blast_folder.set(filename)
	
def Get_Blast_Folder() :
	if(os.path.isfile("init.dat")) :
		pf = open("init.dat", "r")
		init_blast_folder = pf.readline()
		pf.close()
		if (not os.path.isdir(init_blast_folder)) :
			init_blast = os.getcwd()
			
	else :
		init_blast_folder = os.getcwd()
		
	return init_blast_folder

def Save_Blast_Folder(the_line) :
	try :
		pf = open("init.dat", "w")
		pf.write(the_line)
		pf.close()
	except Exception :
		pass
	
	

def makeform(root, tab_parent, tab_connection, tab_blast, fields):
	global message
	global blast_folder
	global init_blast_folder
	entries = {}
	counter = 0
	
	blast_folder = StringVar()
	
	init_blast_folder = Get_Blast_Folder()
		
	brow = tk.Frame(tab_blast)
	lab = tk.Label(brow, width=10, text="Blast Folder:", anchor="w")
	ent = tk.Entry(brow, width=50, textvariable=blast_folder)
	but = tk.Button(tab_blast, text="Browse", command=browse_button)
	 
	
	brow.pack(side=tk.TOP, padx=5, pady=10)
	lab.pack(side=tk.LEFT)
	ent.pack(side=tk.LEFT, expand=tk.YES, fill=tk.X)
	ent.insert(0, init_blast_folder)
	but.place(in_=tab_blast, bordermode=OUTSIDE, height=30, width=50, x=445, y=5)
	entries['blast_folder'] = ent
	
	
	row = tk.Frame(tab_blast)
	lab = tk.Label(row, width = 10, text="Taxa Table:", anchor="w")
	ent = tk.Entry(row, width = 50)
	ent.insert(0, "taxa")
	row.pack(side=tk.TOP, padx=5, pady=10)
	lab.pack(side=tk.LEFT)
	ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
	entries['taxa_table'] = ent
	
	row = tk.Frame(tab_blast)
	lab = tk.Label(row, width = 10, text="DB Release:", anchor="w")
	ent = tk.Entry(row, width = 50)
	ent.insert(0, "40")
	row.pack(side=tk.TOP, padx=5, pady=10)
	lab.pack(side=tk.LEFT)
	ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
	entries['database_release'] = ent
	
	row = tk.Frame(tab_blast)
	lab = tk.Label(row, width = 10, text="Taxa Chunk:", anchor="w")
	ent = tk.Entry(row, width = 50)
	ent.insert(0, 1000)
	row.pack(side=tk.TOP, padx=5, pady=10)
	lab.pack(side=tk.LEFT)
	ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
	entries['taxa_chunk'] = ent
	
	row = tk.Frame(tab_blast)
	lab = tk.Label(row, width = 10, text="Blast Chunk:", anchor="w")
	ent = tk.Entry(row, width = 50)
	ent.insert(0, 1000)
	row.pack(side=tk.TOP, padx=5, pady=10)
	lab.pack(side=tk.LEFT)
	ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
	entries['blast_chunk'] = ent
	
	row = tk.Frame(tab_connection)
	lab = tk.Label(row, width=13, text="Program To Run:", anchor = "w")
	cb  = ttk.Combobox(row, state="readonly", values=["Sequence Retrieval", "Blast Search"],width=47)
	row.pack(side=tk.TOP, padx=5, pady=10)
	lab.pack(side=tk.LEFT)
	cb.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
	cb.current(0)
	
	entries['program'] = cb	
	
	
	for field in fields:
		counter += 1
		row = tk.Frame(tab_connection)
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
			ent.focus()
		elif (counter == 5) :
			ent.insert(0, "test_idseeqer_temp")
			
		row.pack(side=tk.TOP, padx=5, pady=10)
		#row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=10)
		lab.pack(side=tk.LEFT)
		ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
		entries[field] = ent
	
	row = tk.Frame(tab_connection)
	lab = tk.Label(row, width=13, text="Retrieval Delay:", anchor = 'w')
	ent = tk.Entry(row, width=50)
	ent.insert(0, "2")
	row.pack(side=tk.TOP, padx=5, pady=15)
	lab.pack(side=tk.LEFT)
	ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
	entries['delay'] = ent
	
	row = tk.Frame(tab_connection)	
	arial10b= font.Font(family='Arial', size=10, weight='bold')
	message = tk.Label(row, width=400, height=4, text="", anchor='w', font = arial10b, borderwidth=2, relief="groove")
	row.pack(side=tk.TOP, padx=5, pady=15)
	message.pack(side=tk.LEFT, padx=15, pady=15)
	
	return entries
	
def Get_Inputs() :
	global root, tab_parent, tab_connection, tab_query, tab_parameters, but_process
	fields = ('Host', 'Database', 'User', 'Password', 'Table')
	
	root = tk.Tk()
	
	width = 500
	height = 520
	
	# get screen width and height
	screen_width = root.winfo_screenwidth()
	screen_height = root.winfo_screenheight()

	# calculate position x and y coordinates
	x = ( (screen_width - width) /2)
	y = ( (screen_height - height) /2)

	root.geometry('%dx%d+%d+%d' % (width, height, x, y)) # Windows width, height and position on screen.
	root.title('IDSeeqer - [Input Parameters]')
	root.resizable(False, False) # Window not resizable in both direction
	
	tab_parent = ttk.Notebook(root)
	tab_connection = ttk.Frame(tab_parent)
	tab_blast = ttk.Frame(tab_parent)
	tab_parameters = ttk.Frame(tab_parent)
	
	ents = makeform(root, tab_parent, tab_connection, tab_blast, fields)
	
	arial10b = font.Font(family='Arial', size=10, weight='bold')
	
	but_process = tk.Button(tab_connection, default='active', text='Process', command=(lambda e=ents: Correct_Inputs(root, e)), height = 2, width = 10, font = arial10b)
	but_process.pack(side=tk.LEFT, padx=15, pady=15)
	
	but_cancel = tk.Button(tab_connection, text='Cancel', command=(lambda : Cancel(root)), height = 2, width = 10, font = arial10b)
	but_cancel.pack(side=tk.RIGHT, padx=15, pady=15)
	
	#root.bind('<Return>', (lambda e=ents, but_process=but_process: but_process.invoke()))
	root.bind("<Return>", (lambda e=ents: Enter_Hit()))
	
	tab_parent.add(tab_connection, text = " Connection ")
	tab_parent.add(tab_blast,      text = "   Blast    ")
	tab_parent.add(tab_parameters, text = " Parameters ")

	tab_parent.pack(expand = 1, fill = 'both')
	
	#root.protocol("WM_DELETE_WINDOW", Cancel(root)) # Catch the close (X) button
	
	
	try : 
		dirpath = os.getcwd()
		imgicon = tk.PhotoImage(file= dirpath + '/img/icon2.gif')
		root.tk.call('wm', 'iconphoto', root._w, imgicon)
	except Exception :
		pass
		
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
		blast_folder = args[10].strip()
		
		if( not(type(args) is list) or (len(args) < 11) or (not check_delay(delay) ) or \
			(not os.path.isdir(blast_folder)) )  :
			#print("Please, provide host, database, user, password and table as input")
			used_gui = True
			Get_Inputs()
						
			#sys.exit()
		else :
			delay = args[5]
			taxa_table = args[6]
			db_release = int(arg[7])
			taxa_chunk = int(args[8])
			blast_chunk = int(args[9])
			
			del args[5] # Remove delay
			del args[7] # Remove taxa_chunk
			
			
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
	
	if (program_to_run == "Sequence Retrieval") :
		counter = 0
		total = len(res)
		
		t0 = 0
		acc_time = 0
		
		# Loop over all ids and search
		for id_paper1 in res:
			counter += 1
			
			#if(counter > 5) :  # The first 5 records, for testing.
				#break
			
			t0 = timer()
			
			if( acc_time != 0) :
				rem = (total - (counter - 1)) * (acc_time/(counter - 1))
				m_rem, s_rem = divmod(int(rem), 60)
				h_rem, m_rem = divmod(m_rem, 60)
				
				print("\n### APPROX. TIME REMAINING: {:d}hour {:02d}min {:02d}s".format(h_rem, m_rem, s_rem))
				
					
			print("\n>>> Searching Item {} out of {}.".format(counter, total))
		
			for fun in function_names :
				
				print("> Searching '{}' ...".format(display_names[fun].upper()))
				
				if( Run_Function(fun, id_paper1, function_params[fun]) ) :
					
					# func is the last function in list where search was successful
					# Remove from the list
					function_names.remove(fun)
					# Now insert it at the beginning of the list.
					function_names.insert(0, fun)
					
					time.sleep(delay) # Wait for delay seconds after successful retrieval
					break 
				else : 
					max_fail -= 1
			
			temp_time = timer() - t0
			acc_time += temp_time
			
			if(max_fail <= 0) :
				print("\n > Reached ", settings.gramene_params.fail_num, " fails.\nExitting !\n")
				break
			
	elif (program_to_run == "Blast Search") :
		# DO THE BLAST RETRIEVAL
		Blast_Retrieval.Retrieval(taxa_chunk)
	
				
	print(' > Program Completed ')