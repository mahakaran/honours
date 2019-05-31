import multiprocessing as mp
import MDAnalysis

###########################################
##
## This function will extract each frame
## of a MD and then will quit the binary
## part
##
###########################################
MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False

def parallelExtract(frames):
	
	#now we will create a universe and we will loop over timesteps

	universe=MDAnalysis.Universe(gro, xtc)
	selection=None
	try:
		selection=universe.selectAtoms(sel) #in v<0.11 we need to use select_atoms
	except:
		selection=universe.select_atoms(sel)
	for ID, frame in frames:
		universe.trajectory[frame] #will read actual frame
		writerPDB=MDAnalysis.Writer(outputFolder+"/.temp_RIP-MD/frame_"+str(ID)+".pdb")
		writerPDB.write(selection)
	return


#This code is cooked. Essentially, the function below generates a symmetric list of only half of the 
#trajectory, because it skips every nth frame given by flag -sf. This is super stupid because it means that the 
#edge IDs that are printed out don't match the actual frame number. Ths should be fixed. I would fix it but I 
#don't want to mess up the entire code. 

def xtc(GRO, XTC, folder, frameSeparation,frameStart,frameEnd, nproc, Selection):
	global xtc, gro, outputFolder, sel
	sel=Selection
	xtc=XTC
	gro=GRO
	outputFolder=folder
	pdbFrameName={}
	
	u=MDAnalysis.Universe(gro, xtc)
	#########################################
	##conditions for frameStart
	#########################################
	start=None
	if int(frameStart)==0:
		start=0
	else:
		if int(frameStart)>len(u.trajectory):
			print "Error: Start frame is bigger than the number of frames that compose MD trajectory, exiting..."
			exit()
			
		else:
			start=int(frameStart)


	#########################################
	##conditions for frameEnd
	#########################################
	end=None
	if int(frameEnd)==-1:
		end= len(u.trajectory)
	else:
		if int(frameEnd)>len(u.trajectory):
			print "Error: End frame is bigger than the number of snapshots that compose MD trajectory ("+str(len(u.trajectory))+" frames)... exiting"
			exit()
		else:
			end= int(frameEnd)
	
	###############################
	## comparing start number
	## with end number
	###############################
	if start>=end:
			print "Error: Start frame is bigger or equal than end frame... exiting"
			exit()

	i=start
	if int(frameSeparation)>=len(u.trajectory):
		print "Error: Separation frame value is equal or higher than the trajectory length ("+str(len(u.trajectory))+" frames)... exiting"
		exit()
		
	j=0
	framesToExtract=[]
	while i<end:
		framesToExtract.append([j,i]) #j will be the frame id and i is the frame to extract
		i+=int(frameSeparation)+1
		j+=1
	
	print framesToExtract	


	######################################
	## defining subgroups to extract
	## framesin a parallel way
	#####################################
	listToWork=[]
	for i in range(int(nproc)):
		listToWork.append([])
	j=0
	for i in range(len(framesToExtract)):
		listToWork[j].append(framesToExtract[i])
		j+=1
		if j==int(nproc):
			j=0

	if len(framesToExtract)==1 or len(framesToExtract)==0:
		print "There is not enough frames to compute a consensus of interactions... exiting"
		exit()

	pool=mp.Pool(processes=int(nproc)) #for multiprocessing
	listOfListOfNames=pool.map(parallelExtract,(listToWork))
	for vector in framesToExtract:
		pdbFrameName["frame_"+str(vector[0])]=outputFolder+"/.temp_RIP-MD/frame_"+str(vector[0])+".pdb"
	return pdbFrameName
