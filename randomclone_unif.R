#' Method that does a simple random number of clusters and assigns SNVs to closest cluster
#' 
#' Rscript randomclone_unif.R 1e27cc8a-5394-4958-9af6-5ece1fe24516 1e27cc8a-5394-4958-9af6-5ece1fe24516_allDirichletProcessInfo.txt 0.77 GBM-US output/unif/
#' 
# set.seed(123) 
source("~/repo/randomclone/util.R")
MIN_CLUSTERS = 1
MAX_CLUSTERS = 5

args = commandArgs(T)
samplename = args[1]
dpclustinput_infile = args[2]
purity = as.numeric(args[3])
project = args[4]
outdir = args[5]

#' Load the data
dat = parse_data(dpclustinput_infile)
if (nrow(dat) < 2) {
  print("Not enough SNVs")
  q(save="no")
}

res = randomclone_unif(dat)

#' Write the output
write_output_calibration_format(samplename, dat, res$structure, res$assignments, purity, outdir)
write_output_summary_table(res$structure, outdir, samplename, project, purity)

