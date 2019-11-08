import os, logging, subprocess
import PSPFunc
from ete3 import EvolTree

def bppSite(bppFile, bppMixed, alnFile, alnFormat, treeFile, lModels, outDir, baseName, logger):	
	### SITE ANALYSIS: BIO++

	logger.info("Bio++ Site Analysis")
	logger.info("Models to be run: {:s}".format(", ".join(model for model in lModels)))
	logger.info("Bppml parameter file: {:s}".format(bppFile))
	logger.info("Bppmixedlilkelihood parameter file: {:s}".format(bppMixed))
	
	nodes = PSPFunc.nbNode(treeFile, logger)
	## Bppml
	""" 
	Optimize tree and model using bppml
	Variables to include are
		INPUTFILE - alignement file
		FORMAT - format of the aln file (here, phyx)
		TREEFILE - tree file for the analyzed aln
		MODEL - choose which model you want run on the data YNGKP_M0 through 8, same models as PAML
		NODES - number of nodes in the tree file
		IGNORE - parameters to ignore for optimization, for example if one is fixated (ex: omegas in M8a)
		OUTTREE - name of the optimized output tree
		OUTPARAMS - name of the output file summarizing parameters
		BACKUP - name of log file
	"""

	# Bppml output file names - dictionaries that associate model number with output file name for the model
	outSite = outDir+"bpp_site/"
	if not os.path.exists(outSite):
		subprocess.Popen("mkdir "+outSite, shell =  True).wait()

	outFileName = outSite+baseName
	dModelTrees = {model:outFileName+"_"+model+".dnd" for model in lModels}
	dModelParams = {model:outFileName+"_"+model+".params" for model in lModels}
	dModelLog = {model:outFileName+"_optimization"+model for model in lModels}
	dModelSyntax = {model:"YNGP_"+model+"(frequencies=F3X4(initFreqs=observed))" for model in lModels}		# dictionary model number - MODEL argument for bppml
	dLogLlh = {}		# dictionary(model:logllh)
	
	for model in lModels:
		print(model)
		# take into account the specificities of each model (number of classes n for example)
		if int(model[1]) > 2:
			dModelSyntax[model] = dModelSyntax[model].replace(model+"(", model+"(n=8, ")
		if len(model) > 2:
			dModelSyntax[model] = dModelSyntax[model].replace(model, model[:-1]).replace(",", ", omegas=1, ")
			
		# if M0 optimization in models, use tree optimized in M0 for subsequent model optimizations
		if "M0" in lModels:
			if not model == "M0":
				treeFile = dModelTrees["M0"]
			if model == "M8a":
				ignore = "BrLen,YNGP_M8.omegas*"
			elif model == "M0":
				ignore = ""
			else:
				ignore = "BrLen"
		else:
			if model == "M8a":
				ignore = "YNGP_M8.omegas*"
			else:
				ignore = ""

		# create dictionary with all elements of the two argument lists to build commands
		dBppCmd = {"INPUTFILE":alnFile, "FORMAT":alnFormat, "TREEFILE":treeFile, "MODEL":dModelSyntax[model], "NODES":nodes, "IGNORE":ignore, "OUTTREE":dModelTrees[model], "OUTPARAMS":dModelParams[model], "BACKUP":dModelLog[model], "param":bppFile}
		
		# running bppml
		logger.info("Running {:s} optimization".format(model))
		
		# join each couple of the cmd dictionary so that it reads "k1 = v1" "k2 = v2" etc...
		argsMx = "\""+"\" \"".join([k+"="+v for k, v in dBppCmd.items()])+"\""
		logger.debug("bppml "+argsMx)
		runMx = subprocess.Popen("bppml "+argsMx, shell=True, stdout=subprocess.PIPE).wait()
		logger.debug(subprocess.PIPE)
		
		# fill dictionary with loglikelihoods of each model
		if os.path.exists(dModelParams[model]):
			with open(dModelParams[model], "r") as params:
				dLogLlh[model] = float(params.readline().strip().split("= ")[-1])
				logger.info("Log Likelihood = {}".format(dLogLlh[model]))
		else:
			logger.info("Possible failed optimization, likelihood has not been calculated.")
		
	
	# perform LRT
	# M1 vs M2
	if "M1" and "M2" in lModels:
		if "M1" and "M2" in dLogLlh:
			LR12, p12 = PSPFunc.LRT(dLogLlh["M1"], dLogLlh["M2"], 2)
			logger.info("LRT of M1 vs M2: {}".format(p12))
		else:
			logger.info("Possible failed optimization, likelihoods of M1 and M2 have not been computed.")
	if "M7" and "M8" in lModels:
		if "M7" and "M8" in dLogLlh:
			LR78, p78 = PSPFunc.LRT(dLogLlh["M7"], dLogLlh["M8"], 2)
			logger.info("LRT of M7 vs M8: {}".format(p78))
		else:
			logger.info("Possible failed optimization, likelihoods of M7 and M8 have not been computed.")
	if "M8" and "M8a" in lModels:
		if "M8" and "M8a" in dLogLlh:
			LR88a, p88a = PSPFunc.LRT(dLogLlh["M8a"], dLogLlh["M8"], 1)
			ts88a = 0.5*p88a + 0.5
			logger.info("LRT of M8 vs M8a: {} (Treshold: {})".format(p88a, ts88a))
		else:
			logger.info("Possible failed optimization, likelihoods have not been computed.")
			
	# Bppmixedlikelihoods
	""" 
	Optimize tree and model using bppml
	Variables to include are
		INPUTFILE - alignement file
		FORMAT - format of the aln file (here, phyx)
		TREEFILE - tree file for the analyzed aln
		PARAMS - .params file from model optimization (bppml)
		OUTINFO - name of the results file (info about sites etc.)
	"""
	
	for model in lModels:
		# use tree optimized in M0 for each model
		if "M0" in lModels:
			treeFile = dModelTrees["M0"]
		else:
			treeFile = dModelTrees[model]
		
		if model == "M0":
			logger.info("M0 is not a mixed model, skipping.")
			continue
			
		# dictionary(model:results file name)
		dModelResults = {model:outDir+baseName+"_results"+model+".log" for model in lModels}

		dMixCmd = {"INPUTFILE":alnFile, "FORMAT":alnFormat, "TREEFILE":treeFile, "PARAMS":dModelParams[model], "OUTINFO":dModelResults[model], "param":bppMixed}
		
		logger.info("Running mixed likelihoods with model {:s}".format(model))
		argsMx = "\""+"\" \"".join([k+"="+v for k, v in dMixCmd.items()])+"\""
		logger.debug("bppmixedlikelihoods "+argsMx)
		runMx = subprocess.Popen("bppmixedlikelihoods "+argsMx, shell=True, stdout=subprocess.PIPE).wait()
		logger.debug(subprocess.PIPE)
			

def pamlSite(alnFile, treeFile, lModels, pamlParams, outDir, baseName, logger):
	
	tree = EvolTree(treeFile)
	os.mkdir(outDir+"paml_site/")
	tree.workdir = outDir+"paml_site/"
	tree.link_to_alignment(alnFile, "Fasta")
	logger.info("PAML codeml")
	
	dModelRun = {}
	if "M8a" in lModels:
		lModels.remove("M8a")
	for model in lModels:
		logger.info("Running {:s}".format(model))
		dModelRun[model] = tree.run_model(model)
	
	if "M1" and "M2" in dModelRun:
		p12 = tree.get_most_likely("M2", "M1")
		logger.info("LRT of M1 vs M2 = {}".format(p12))
	if "M7" and "M8" in dModelRun:
		p78 = tree.get_most_likely("M8", "M7")
		logger.info("LRT of M7 vs M8 = {}".format(p78))
	"""
	if "M8a" and "M8" in dModelRun:
		p88a = tree.get_most_likely("M8a", "M8")
		logger.info("LRT of M8 vs M8a = {}".format(p88a))
	"""
	