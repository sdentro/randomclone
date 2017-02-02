#' Method that takes mutations and assigns all of them to 1 cluster and takes the median CCF of the SNVs as the cluster CCF
#' 
#' Rscript single_cluster.R 1e27cc8a-5394-4958-9af6-5ece1fe24516 1e27cc8a-5394-4958-9af6-5ece1fe24516_allDirichletProcessInfo.txt 0.77 GBM-US $PWD/



args = commandArgs(T)
samplename = args[1]
dpclustinput_infile = args[2]
purity = as.numeric(args[3])
project = args[4]
outdir = args[5]

# samplename = "1e27cc8a-5394-4958-9af6-5ece1fe24516"
# dpclustinput_infile = "1e27cc8a-5394-4958-9af6-5ece1fe24516_allDirichletProcessInfo.txt"
# purity = 0.77

dat = read.table(dpclustinput_infile, header=T, stringsAsFactors=F)
dat = dat[!is.na(dat$subclonal.fraction),]
cluster_ccf = median(dat$subclonal.fraction, na.rm=T)
cluster_cp = cluster_ccf * purity

purity_df = data.frame(sample=samplename, 
                       purity=purity, 
                       ploidy=NA)
assignments_df = data.frame(chr=dat$chr, 
                            pos=dat$end, 
                            cluster=rep(1, nrow(dat)))
structure_df = data.frame(cluster=1, 
                          n_ssms=nrow(assignments_df), 
                          proportion=cluster_cp,
                          ccf=cluster_ccf)
num_clusters_df = data.frame(sample=samplename, 
                             clusters=nrow(structure_df))
multiplicity_df = data.frame(chr=dat$chr, 
                             pos=dat$end, 
                             tumour_copynumber=(dat$nMaj1+dat$nMin1)*dat$frac1 + ifelse(!is.na(dat$frac2), (dat$nMaj2+dat$nMin2)*dat$frac2, 0), 
                             multiplicity=dat$no.chrs.bearing.mut, 
                             multiplicity_options=NA, 
                             probabilities=NA)

#' Write out the calibration format data
write.table(purity_df, file=file.path(outdir, paste0(samplename, "_purity_ploidy.txt")), row.names=F, sep="\t", quote=F)
write.table(assignments_df, file=file.path(outdir, paste0(samplename, "_mutation_assignments.txt")), row.names=F, sep="\t", quote=F)
write.table(structure_df, file=file.path(outdir, paste0(samplename, "_subclonal_structure.txt")), row.names=F, sep="\t", quote=F)
write.table(num_clusters_df, file=file.path(outdir, paste0(samplename, "_number_of_clusters.txt")), row.names=F, sep="\t", quote=F)
write.table(multiplicity_df, file=file.path(outdir, paste0(samplename, "_multiplicity.txt")), row.names=F, sep="\t", quote=F)

#' Determine which cluster is clonal
is_clonal = structure_df$ccf > 0.9
if (any(is_clonal)) {
  if (sum(is_clonal, na.rm=T) > 1) {
    clone_id = which.min(abs(structure_df$ccf[is_clonal]-1))
    is_clonal = rep(F, length(is_clonal))
    is_clonal[clone_id] = T
  }
  clonal_cluster = structure_df$cluster[is_clonal]
} else {
  clonal_cluster = NA
}

if (!is.na(clonal_cluster)) {
  num_clonal = structure_df$n_ssms[structure_df$cluster==clonal_cluster]
  num_subclonal = ifelse(nrow(structure_df) > 1, structure_df$n_ssms[structure_df$cluster!=clonal_cluster], 0)
} else {
  num_clonal = 0
  num_subclonal = sum(structure_df$n_ssms)
}


#' Write out the summary table entry
summary_table = data.frame(cancer_type=project, 
                           samplename=samplename, 
                           num_subclones=nrow(structure_df)-1, 
                           purity=purity, 
                           ploidy=NA,
                           num_clonal=num_clonal,
                           num_subclonal=num_subclonal,
                           frac_clonal=num_clonal / (num_clonal+num_subclonal),
                           noCNA=NA,
                           clonal=NA,
                           subclonal=NA)
write.table(summary_table, file=file.path(outdir, paste0(samplename, "_summary_table.txt")), row.names=F, sep="\t", quote=F)




