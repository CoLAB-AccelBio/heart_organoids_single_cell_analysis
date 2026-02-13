################################################################################
#
#   File name: trajectory_inference.R
#
#   Authors: Jacek Marzec ( jacek.marzec@accelbio.pt )
#
#   Biocant Park,
#   Parque Tecnológico de Cantanhede,
#   3060-197 Cantanhede
#
################################################################################

################################################################################
#
#	  Description: Script collecting user-defined parameters for the corresponding trajectory_inference_heart_organoid_S1_3_endothelial_cells.Rmd markdown script performing ....
#
#	  Command line use example: Rscript trajectory_inference.R  ...
#
#   inputFolder:    Full path with name of the folder with input data
#   start_clust:    Starting cluster(s) from which lineages will be drawn. Default is NULL
#   clust_res:      Resolution from the clustering step to be used for the analysis. Default is 0.1
#   omega:          Logical to introduce an artificial cluster omega. Default is TRUE
#   genes:          List of genes for exploration
#   run_slingshot   Logical to run Slingshot. Default is TRUE
#   run_tscan       Logical to run TSCAN. Default is TRUE
#   run_monocle3    Logical to run Monocle 3. Default is TRUE
#   outFolder:      Desired location for the results and report
#   report_name:    Desired core name for the results
#   seed_init:      Set up a seed for random number generation
#   knitr_verbose:  Logical to run knitr in verbose mode to get more detailed output about what is being processed. Default is FALSE
#   hide_code_btn:  Hide the "Code" button allowing to show/hide code chunks in the final HTML report. Default is TRUE
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Load libraries
#===============================================================================

suppressMessages(library(optparse))

#===============================================================================
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option("--inFolder", action="store", default=NA, type='character',
              help="Full path with name of the folder with input data"),
  make_option("--start_clust", action="store", default=NULL, type='character',
              help="Starting cluster(s) from which lineages will be drawn"),
  make_option("--clust_res", action="store", default=0.1, type='numeric',
              help="Resolution from the clustering step to be used for the analysis"),
  make_option("--omega", action="store", default=TRUE, type='logical',
              help="Logical to introduce an artificial cluster omegat"),
  make_option("--genes", action="store", default=NA, type='character',
              help="List of genes to explore"),
  make_option("--run_slingshot", action="store", default=TRUE, type='logical',
              help="Logical to run Slingshot"),
  make_option("--run_tscan", action="store", default=TRUE, type='logical',
              help="Logical to run TSCAN"),
  make_option("--run_monocle3", action="store", default=TRUE, type='logical',
              help="Logical to run Monocle 3"),
  make_option("--outFolder", action="store", default=NA, type='character',
              help="Desired location for the results and report"),
  make_option("--report_name", action="store", default=NA, type='character',
              help="Desired core name for the results"),
  make_option("--seed_init", action="store", default=99999999, type='numeric',
              help="Set up a seed for random number generation"),
  make_option("--knitr_verbose", action="store", default=FALSE, type='logical',
              help="Logical to run knitr in verbose mode to get more detailed output about what is being processed"),
  make_option("--hide_code_btn", action="store", default=FALSE, type='logical',
              help="Hide the \"Code\" button allowing to show/hide code chunks in the final HTML report")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Read in argument from command line and check if all required arguments were provide by the user
if ( is.na(opt$omega) ) {
  
  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript trajectory_inference.R ...\n\n")
  q()
  
}

##### Create user-defined directory for the report
if ( !file.exists("~/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_endothelial_cells/trajectory_inference") ) {
  dir.create("~/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_endothelial_cells/trajectory_inference", recursive=TRUE)
}

##### Collect parameters
param_list <- list(inFolder = opt$inFolder,
                   start_clust = opt$start_clust,
                   clust_res = opt$clust_res,
                   omega = opt$omega,
                   genes = opt$genes,
                   run_slingshot = opt$run_slingshot,
                   run_tscan = opt$run_tscan,
                   run_monocle3 = opt$run_monocle3,
                   outFolder = opt$outFolder,
                   report_name = opt$report_name,
                   seed_init = opt$seed_init,
                   knitr_verbose = opt$knitr_verbose,
                   hide_code_btn = opt$hide_code_btn)

##### Pass the user-defined arguments to the "trajectory_inference.Rmd" R markdown script and generate the report
rmarkdown::render(input = "trajectory_inference.Rmd",
                  output_file = paste0(opt$output_dir, "/", opt$report_name, ".html"),
                  output_dir = opt$outFolder,
                  params = param_list )

##### Remove the assocaited MD file and the redundant folder with plots that are embedded in the HTML report
unlink(paste0(opt$outFolder, "/trajectory_inference.md"), recursive = TRUE)
unlink(paste0(opt$outFolder, "/trajectory_inference_files"), recursive = TRUE)

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
