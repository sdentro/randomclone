#' Method that takes mutations and assigns all of them to 1 cluster and takes the median CCF of the SNVs as the cluster CCF
#' 
#' Rscript randomclone_single.R 1e27cc8a-5394-4958-9af6-5ece1fe24516 1e27cc8a-5394-4958-9af6-5ece1fe24516_allDirichletProcessInfo.txt 0.77 GBM-US output/single/
FORCE_CLONE = T

args = commandArgs(T)
libpath = args[1]
samplename = args[2]
dpclustinput_infile = args[3]
purity = as.numeric(args[4])
project = args[5]
outdir = args[6]

source(file.path(libpath,"util.R"))
#' Load the data
dat = parse_data(dpclustinput_infile)
if (nrow(dat) < 2) {
  print("Not enough SNVs")
  q(save="no")
}

if (!FORCE_CLONE) {
  cluster_ccf = median(dat$subclonal.fraction, na.rm=T)
} else {
  cluster_ccf = 1
}
cluster_cp = cluster_ccf * purity

structure_df = data.frame(cluster=1, 
                          n_ssms=nrow(dat), 
                          proportion=cluster_cp,
                          ccf=cluster_ccf)

write_output_calibration_format(samplename, dat, structure_df, rep(1, nrow(dat)), purity, outdir)
write_output_summary_table(structure_df, outdir, samplename, project, purity)


