import PSPFunc, GeneAnalysis, SiteAnalysis, BranchAnalysis, os
from time import localtime, strftime
import logging
import subprocess, re
import Init
from Logging import setup_logger

def pspAnalysis(params, rule):
    """
    procedure which execute functions for psp step

    @param1 data: basicData object

        @return Output directory name
    """
    logger=logging.getLogger(".".join(["main","positiveSelection",rule]))

    tree = os.path.join(params["outdir"],params["queryName"]+"_tree.dnd")
    aln = os.path.join(params["outdir"],params["queryName"]+"_align.fasta")
    
    timeStamp = strftime("%Y%m%d%H%M", localtime())

    ### set up new directory
    
    outDir = os.path.join(params["outdir"],params["queryName"] + "_positive_selection")
    os.makedirs(outDir, exist_ok=True)

    params["outDir"]=outDir
    
    ### Terminal output for user
    logger.info("Output directory: {:s}".format(outDir))
    logger.info("Alignement: {:s}".format(aln))
    logger.info("Tree: {:s}".format(tree))

    ### Run the different analysis as determined by control file
    logger.info("Starting positive selection analyses.")
    logger.info("POSITIVE SELECTION ANALYSIS: ")
    logger.info("Analysis to be run:")

    ###########################################################
    #### HYPHY

    if rule in ["busted","meme"]:
      cladoFile =  PSPFunc.supBoot(params)
      if rule == "busted":
        try:		
          GeneAnalysis.hyphyBusted(aln, cladoFile, outDir, logger)
        except Exception:
          logger.info("BUSTED encountered an unexpected error, skipping.")

      if rule=="meme":
        try:
          BranchAnalysis.memeBranchSite(aln, cladoFile, outDir, logger)
        except Exception:
          logger.error("MEME encountered an unexpected error, skipping.")


    ###########################################################
    #### BPP

    if rule in ["bppml","opb","paml"]:
      lModels = list(map(str.strip,re.compile(r"[,;]").split(params["models"])))

    if rule in ["bppml","opb"]:
      outBPP = outDir+"/bpp_site"
      os.makedirs(outBPP, exist_ok=True)
      
      PSPFunc.pspFileCreation(outBPP+"/base.bpp","bppml")
      
      if params["mixedlikelihood"]:
        params["mixedlikelihood"]=outBPP+"/base_mll.bpp"
        PSPFunc.pspFileCreation(params["mixedlikelihood"],"bppmixedlikelihood")
      else:
        params["mixedlikelihood"]=outBPP+"/base.bpp"

      SiteAnalysis.bppSite(aln, 
                           tree, 
                           outDir, 
                           outBPP+"/base.bpp",
                           params["mixedlikelihood"], 
                           lModels, 
                           logger)
    
      lPSNodes = []
      if rule == "opb":
        PSPFunc.pspFileCreation(outBPP+"/base_opb.bpp","opb")
        try:
          params = BranchAnalysis.bppBranch(aln, 
                                            tree, 
                                            outDir, 
                                            outBPP+"/base_opb.bpp", 
                                            logger)	
        except Exception:
          logger.error("Bio++ Branch Analysis encountered an unexpected error, skipping.")
    
      # if params["opb"] and params["gnh"] and len(lPSNodes) > 1:
      # params["gnh"]=outBPP+"/base_gnh.bpp"
      # PSPFunc.pspFileCreation(params["gnh"],"gnh")
      # try:
      #       BranchAnalysis.bppBranchSite(aln, tree, outDir, params["gnh"], lPSNodes, logger)
      # except Exception:
      #  	logger.error("Bio++ Pseudo Branch-Site Analysis encountered an unexpected error, skipping.")


    ###########################################################
    #### PAML
        
    if rule=="paml":
      SiteAnalysis.pamlSite(aln, 
                            tree, 
                            outDir, 
                            params["paml"], 
                            lModels,
                            logger)

      # """try:
      # SiteAnalysis.pamlSite(aln, tree, lModels, dCtrls["paml"], outDir, data.baseName, logger)
      # except Exception:
      # logger.info("PAML (codeml) Site encountered an unexpected error, skipping.")"""

    logger.info("Finished positive selection analyses.")


if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    # Logging
    logger = logging.getLogger("main.positiveSelection."+snakemake.rule)
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)

    aln = os.path.join(config["outdir"], config["queryName"] + "_align.fasta")

    # Run step
    parameters = Init.paramDef(config)
    pspAnalysis(parameters,snakemake.rule)

