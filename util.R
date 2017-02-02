
parse_data = function(dpclustinput_infile) {
  dat = read.table(dpclustinput_infile, header=T, stringsAsFactors=F)
  dat = dat[!is.na(dat$subclonal.fraction),]
  return(dat)
}

write_output_calibration_format = function(samplename, dat, structure_df, assignments, purity, outdir) {
  
  assignments_df = data.frame(chr=dat$chr, 
                              pos=dat$end, 
                              cluster=assignments)
  num_clusters_df = data.frame(sample=samplename, 
                               clusters=nrow(structure_df))
  multiplicity_df = data.frame(chr=dat$chr, 
                               pos=dat$end, 
                               tumour_copynumber=(dat$nMaj1+dat$nMin1)*dat$frac1 + ifelse(!is.na(dat$frac2), (dat$nMaj2+dat$nMin2)*dat$frac2, 0), 
                               multiplicity=dat$no.chrs.bearing.mut, 
                               multiplicity_options=NA, 
                               probabilities=NA)
  purity_df = data.frame(sample=samplename, 
                         purity=purity, 
                         ploidy=NA)
  
  #' Write out the calibration format data
  write.table(purity_df, file=file.path(outdir, paste0(samplename, "_purity_ploidy.txt")), row.names=F, sep="\t", quote=F)
  write.table(assignments_df, file=file.path(outdir, paste0(samplename, "_mutation_assignments.txt")), row.names=F, sep="\t", quote=F)
  write.table(structure_df, file=file.path(outdir, paste0(samplename, "_subclonal_structure.txt")), row.names=F, sep="\t", quote=F)
  write.table(num_clusters_df, file=file.path(outdir, paste0(samplename, "_number_of_clusters.txt")), row.names=F, sep="\t", quote=F)
  write.table(multiplicity_df, file=file.path(outdir, paste0(samplename, "_multiplicity.txt")), row.names=F, sep="\t", quote=F)
}

write_output_summary_table = function(structure_df, outdir, samplename, project, purity) {
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
}

mutationCopyNumberToMutationBurden<-function(copyNumber,totalCopyNumber,cellularity,normalCopyNumber = rep(2,length(copyNumber))){
  #burden = copyNumber*cellularity/(cellularity*(totalCopyNumber-burden)+2*(1-cellularity))
  burden = copyNumber*cellularity/(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
  #burden[is.nan(burden)]=0
  burden[is.nan(burden)|(burden<0.000001)]=0.000001
  #190412
  #burden[burden>1]=1
  burden[burden>0.999999]=0.999999
  return(burden)  
}

calc.new.likelihood = function(y, n, kappa, thetas) {
  lfoy = log.f.of.y(y,n,kappa,thetas)
  new.likelihood = sum(lfoy)
  return(new.likelihood)
}

log.f.of.y <- function(y1, n1, kappa1, x) {
  #x=1 and kappa=1 causes problems
  x[x>0.999 & kappa1==1] = 0.999
  #allow kappa = 0, for mutations on deleted chromosomes
  if (class(kappa1) == 'numeric') {
    # Case input consists of vectors
    no.kappa.nonzero = length(which(kappa1!=0))
    no.subsamples = length(y1)
  } else {
    # Case input consists of matrices
    no.kappa.nonzero = rowSums(kappa1!=0)
    no.subsamples = ncol(y1)
  }
  kappa1[kappa1==0] = NA
  res = lchoose(n1, y1) + y1 * log(kappa1*x) + (n1-y1) * log(1-kappa1*x)
  if (class(res) == "numeric") { res = matrix(res, nrow=1) }
  resSums = rowSums(res,na.rm=T) * no.subsamples/no.kappa.nonzero
  return(resSums)
}

aic <- function(likelihood, num.samples, num.nodes) {
  return(- 2 * likelihood + 2 * num.samples * num.nodes)
}

bic <- function(likelihood, num.samples, num.nodes, num.muts) {
  return(-2 * likelihood + 2 * num.samples * num.nodes * log(num.muts))
}

binom_ll = function(cluster_location, mutcount, wtcount, tumourCopyNumber, copyNumberAdjustment, purity, normalCopyNumber) {
  mutBurdens = mutationCopyNumberToMutationBurden(cluster_location * copyNumberAdjustment, tumourCopyNumber, purity, normalCopyNumber)
  assignment_ll = sapply(1:length(mutcount), function(k, mc, wt, mb) {  mc[k]*log(mb[k]) + wt[k]*log(1-mb[k]) }, mc=mutcount, wt=wtcount, mb=mutBurdens)
  return(sum(assignment_ll))
}

calc_all_metrics = function(dat, purity, res) {
  kappa = mutationCopyNumberToMutationBurden(1, dat$subclonal.CN, purity) * dat$no.chrs.bearing.mut
  num_muts = nrow(dat)
  num_samples = 1
  likelihoods = rep(NA, ITERATIONS)
  aics = rep(NA, ITERATIONS)
  bics = rep(NA, ITERATIONS)
  binom_ll = rep(NA, ITERATIONS)
  for (i in 1:ITERATIONS) {
    structure_df = res[[i]]$structure
    likelihoods[i] = calc.new.likelihood(dat$mut.count, (dat$mut.count+dat$WT.count), kappa, structure_df$ccf)
    aics[i] = aic(likelihoods[i], num_samples, nrow(structure_df))
    bics[i] = bic(likelihoods[i], num_samples, nrow(structure_df), log(num_muts))
    
    binom_ll_sample = 0
    for (j in 1:nrow(structure_df)) {
      binom_ll_sample = binom_ll(structure_df$ccf[j], dat$mut.count, dat$WT.count, dat$subclonal.CN, dat$no.chrs.bearing.mut, purity, rep(2, nrow(dat)))
    }
    binom_ll[i] = binom_ll_sample
  }
  all_metrics = data.frame(likelihood=likelihoods, aic=aics, bic=bics, binom_ll=binom_ll)
  return(all_metrics)
}

randomclone_unif = function(dat) {
  #' Draw number of clusters
  n_clusters = sample(MIN_CLUSTERS:MAX_CLUSTERS, 1)
  
  #' Get area where 95% of the data lives and set those boundaries as min and max CCF for clusters
  ccf_boundaries = quantile(dat$subclonal.fraction, probs=c(.025,.975))
  
  ccf_min = ccf_boundaries[[1]]
  ccf_max = ccf_boundaries[[2]]
  
  #' Draw clusters from uniform distribution
  cluster_locations = runif(n_clusters, min=ccf_min, max=ccf_max)
  
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
  return(list(structure=structure_df, assignments=cluster_assignments))
}


