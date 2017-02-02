#' Method that takes mutations and assigns all of them to 1 cluster and takes the median CCF of the SNVs as the cluster CCF
#' 
#' Rscript randomclone_single.R 1e27cc8a-5394-4958-9af6-5ece1fe24516 1e27cc8a-5394-4958-9af6-5ece1fe24516_allDirichletProcessInfo.txt 0.77 GBM-US output/single/
source("~/repo/randomclone/util.R")

args = commandArgs(T)
samplename = args[1]
dpclustinput_infile = args[2]
purity = as.numeric(args[3])
project = args[4]
outdir = args[5]

#' Load the data
dat = parse_data(dpclustinput_infile)
cluster_ccf = median(dat$subclonal.fraction, na.rm=T)
cluster_cp = cluster_ccf * purity

structure_df = data.frame(cluster=1, 
                          n_ssms=nrow(dat), 
                          proportion=cluster_cp,
                          ccf=cluster_ccf)

write_output_calibration_format(samplename, dat, structure_df, rep(1, nrow(dat)), purity, outdir)
write_output_summary_table(structure_df, outdir, samplename, project, purity)


