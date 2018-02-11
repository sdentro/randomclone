#' Method that does a simple random number of clusters and assigns SNVs by drawing random breaks in the CCF sorted data
#' 
#' Rscript randomclone_stick.R 1e27cc8a-5394-4958-9af6-5ece1fe24516 1e27cc8a-5394-4958-9af6-5ece1fe24516_allDirichletProcessInfo.txt 0.77 GBM-US output/stick/
#' 
# set.seed(123)
MIN_CLUSTERS = 1
MAX_CLUSTERS = 5

args = commandArgs(T)
libpath = args[1]
samplename = args[2]
dpclustinput_infile = args[3]
purity = as.numeric(args[4])
project = args[5]
outdir = args[6]

source(file.path(libpath, "util.R"))

#' Load the data
dat = parse_data(dpclustinput_infile)
if (nrow(dat) < 2) {
  print("Not enough SNVs")
  q(save="no")
}

# #' Draw number of clusters
# n_clusters = sample(MIN_CLUSTERS:MAX_CLUSTERS, 1)
# 
# #' Sort SNVs by their CCF
# snv_order = order(dat$subclonal.fraction)
# assignments = rep(NA, nrow(dat))
# if (n_clusters > 1) {
#   #' Put n-1 breaks
#   breaks = sample(1:nrow(dat), n_clusters-1)
#   
#   #' for the n clusters take median CCF as the cluster locations
#   for (i in 1:length(breaks)) {
#     cluster_break = breaks[i]
#     assignments[snv_order < cluster_break] = i
#   }
# } else {
#   cluster_break = 0
#   i = 0
# }
# 
# #' Assign SNVs of the final cluster
# assignments[snv_order > cluster_break] = i+1
# 
# #' Determine cluster locations
# cluster_locations = rep(NA, n_clusters)
# for (i in 1:n_clusters) {
#   cluster_locations[i] = median(dat$subclonal.fraction[assignments==i], na.rm=T)
# }
# cluster_locations = cluster_locations[!is.na(cluster_locations)]
# 
# structure_df = data.frame(table(assignments), 
#                           proportion=cluster_locations * purity,
#                           ccf=cluster_locations)
# colnames(structure_df)[1] = "cluster"
# colnames(structure_df)[2] = "n_ssms"

res = randomclone_stick(dat)

#' Write the output
write_output_calibration_format(samplename, dat, res$structure, res$assignments, purity, outdir)
write_output_summary_table(res$structure, outdir, samplename, project, purity)





