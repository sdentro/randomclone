#' Method that does a simple random number of clusters and assigns SNVs to closest cluster
#' 
#' Rscript randomclone_unif.R 1e27cc8a-5394-4958-9af6-5ece1fe24516 1e27cc8a-5394-4958-9af6-5ece1fe24516_allDirichletProcessInfo.txt 0.77 GBM-US output/unif/
#' 
set.seed(123) 
source("util.R")
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

#' Draw number of clusters
n_clusters = sample(MIN_CLUSTERS:MAX_CLUSTERS, 1)

#' Get area where 95% of the data lives and set those boundaries as min and max CCF for clusters
ccf_boundaries = quantile(dat$subclonal.fraction, probs=c(.025,.975))

ccf_min = ccf_boundaries[[1]]
ccf_max = ccf_boundaries[[2]]

#' Draw clusters from uniform distribution
cluster_locations = runif(n_clusters, min=ccf_min, max=ccf_max)

#cluster_locs = c(0.75, 1.05)

#' Assign SNVs
cluster_assignments = sapply(dat$subclonal.fraction, function(x) { which.min(abs(x-cluster_locations)) })

#' Re-adjust the cluster locations based on assignments
for (i in 1:n_clusters) {
  cluster_locations[i] = median(dat$subclonal.fraction[cluster_assignments==i], na.rm=T)
}

structure_df = data.frame(table(cluster_assignments), 
                          proportion=cluster_locations * purity,
                          ccf=cluster_locations)
colnames(structure_df)[1] = "cluster"
colnames(structure_df)[2] = "n_ssms"

#' Write the output
write_output_calibration_format(samplename, dat, structure_df, cluster_assignments, purity, outdir)
write_output_summary_table(structure_df, outdir, samplename, project, purity)

