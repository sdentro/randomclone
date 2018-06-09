
parse_data = function(dpclustinput_infile) {
  dat = read.table(dpclustinput_infile, header=T, stringsAsFactors=F)
  dat = dat[!is.na(dat$subclonal.fraction),]
  return(dat)
}

write_output_calibration_format = function(samplename, dat, structure_df, assignments, purity, outdir) {
  
  # Catch case where simulated data does not contain these columns
  if (!"frac2" %in% colnames(dat)) {
    dat$nMaj1 = dat$subclonal.CN
    dat$nMin1 = 0
    dat$frac1 = 1
    dat$frac2 = NA
  }
  
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

run_mtimer = function(libpath, mtimerpath, clusters, vcf_snv, bb_file, purity, ploidy, sex, is_wgd, q=0.05, min_read_diff=2, rho_snv=0.01, deltaFreq=0.00, round_subclonal_cn=F, remove_subclonal_cn=F, xmin=3) {
  source(file.path(mtimerpath, "MutationTime.R"))
  #source("~/repo/MutationTime.R/MutationTime.R")
  source(file.path(libpath, "util.R"))
  
  #' reset cluster numbers
  clusters = clusters[with(clusters, order(proportion, decreasing=T)),]
  clusters$cluster = 1:nrow(clusters)
  
  #' Load copy number and variants
  bb <- loadBB(bb_file, round_subclones=round_subclonal_cn, remove_subclones=remove_subclonal_cn)
  bb$clonal_frequency = 1
  vcf_snv <- readVcf(vcf_snv, genome="GRCh37")
  #' Merge too close clusters
  if (nrow(clusters) > 1) { clusters = mergeClustersByMutreadDiff(clusters, purity, ploidy, vcf_snv, min_read_diff) }
  #' Calc assignment probs
  MCN <- computeMutCn(vcf_snv, bb, clusters, purity, gender=sex, isWgd=is_wgd, rho=rho_snv, n.boot=0, xmin=xmin, deltaFreq=deltaFreq)
  #' Calc tail probabilities  
  qq_snv <- mean(MCN$D$pMutCNTail < q/2 | MCN$D$pMutCNTail > 1-q/2, na.rm=T)
  # p_snv = pbinom(sum(MCN$D$pMutCNTail < q/2 | MCN$D$pMutCNTail > 1-q/2, na.rm=T), nrow(MCN$D), 0.05, lower.tail=TRUE)
  return(qq_snv)
}

binom_ll_diff = function(structure_df, assignments, mutcount, wtcount, tumourCopyNumber, copyNumberAdjustment, purity, normalCopyNumber) {
  # assignment_ll = array(NA, c(length(mutcount), nrow(structure_df)))
  # for (j in 1:nrow(structure_df)) {
  #   mutBurdens = mutationCopyNumberToMutationBurden(structure_df$ccf[j] * copyNumberAdjustment, tumourCopyNumber, purity, normalCopyNumber)
  #   # assignment_ll = sapply(1:length(mutcount), function(k, mc, wt, mb) {  mc[k]*log(mb[k]) + wt[k]*log(1-mb[k]) }, mc=mutcount, wt=wtcount, mb=mutBurdens)
  #   assignment_ll[,j] = sapply(1:length(mutcount), function(k, mc, wt, mb) {  mutcount[k]*log(mutBurdens[k]) + wtcount[k]*log(1-mutBurdens[k]) }, mc=mutcount, wt=wtcount, mb=mutBurdens)
  # }
  # 
  # # calc likelihood of assigned cluster minus not assigned cluster
  # ll_diff = rep(NA, nrow(dat))
  # for (j in 1:nrow(dat)) {
  #   assigned = assignments[j]
  #   ll_diff[j] = assignment_ll[j, assigned]
  #   for (c in (1:nrow(structure_df))[-assigned]) {
  #     ll_diff[j] = ll_diff[j] - assignment_ll[j, c]
  #   }
  # }
  # 
  # return(sum(ll_diff))
  return(NA)
}

binom_ll_2 = function(structure_df, mutcount, wtcount, tumourCopyNumber, copyNumberAdjustment, purity, normalCopyNumber) {
  binom_ll_sample = 0
  for (j in 1:nrow(structure_df)) {
    mutBurdens = mutationCopyNumberToMutationBurden(structure_df$ccf[j] * copyNumberAdjustment, tumourCopyNumber, purity, normalCopyNumber)
    binom_ll_sample = binom_ll_sample + sapply(1:length(mutcount), function(k, mc, wt, mb) {  mc[k]*log(mb[k]) + wt[k]*log(1-mb[k]) }, mc=mutcount, wt=wtcount, mb=mutBurdens)
  }
  return(sum(binom_ll_sample))
}



calc_all_metrics = function(util_libpath, mtimer_libpath, dat, purity, res, vcf_snv, bb_file, ploidy, sex, is_wgd, q=0.05, min_read_diff=2, rho_snv=0.01, deltaFreq=0.00, round_subclonal_cn=F, remove_subclonal_cn=F, xmin=3) {

  kappa = mutationCopyNumberToMutationBurden(1, dat$subclonal.CN, purity) * dat$no.chrs.bearing.mut
  num_muts = nrow(dat)
  num_samples = 1
  # likelihoods = rep(NA, ITERATIONS)
  # aics = rep(NA, ITERATIONS)
  # bics = rep(NA, ITERATIONS)
  # binom_lls = rep(NA, ITERATIONS)
  # binom_ll_2s = rep(NA, ITERATIONS)
  # binom_ll_diffs = rep(NA, ITERATIONS)
  # mtimer_ll = rep(NA, ITERATIONS)
  # for (i in 1:ITERATIONS) {
  #	  structure_df = res[[i]]$structure
  #	      assignments = res[[i]]$assignments
  #	  mtimer_ll = run_mtimer(mtimer_libpath, structure_df, vcf_snv, bb_file, purity, ploidy, sex, is_wgd, q=q, min_read_diff=min_read_diff, rho_snv=rho_snv, deltaFreq=deltaFreq, xmin=xmin)
  #}


  scores = mclapply(1:ITERATIONS, function(i) {
    structure_df = res[[i]]$structure
    assignments = res[[i]]$assignments
    likelihoods = calc.new.likelihood(dat$mut.count, (dat$mut.count+dat$WT.count), kappa, structure_df$ccf)
    aics = aic(likelihoods, num_samples, nrow(structure_df))
    bics = bic(likelihoods, num_samples, nrow(structure_df), log(num_muts))
    
    # Original binom_ll with bug
    binom_ll_sample = 0
    for (j in 1:nrow(structure_df)) {
      binom_ll_sample = binom_ll(structure_df$ccf[j], dat$mut.count, dat$WT.count, dat$subclonal.CN, dat$no.chrs.bearing.mut, purity, rep(2, nrow(dat)))
    }
    binom_lls = binom_ll_sample
    
    binom_ll_2s = binom_ll_2(structure_df, dat$mut.count, dat$WT.count, dat$subclonal.CN, dat$no.chrs.bearing.mut, purity, rep(2, nrow(dat)))
    binom_ll_diffs = binom_ll_diff(structure_df, assignments, dat$mut.count, dat$WT.count, dat$subclonal.CN, dat$no.chrs.bearing.mut, purity, rep(2, nrow(dat)))
    mtimer_ll = run_mtimer(util_libpath, mtimer_libpath, structure_df, vcf_snv, bb_file, purity, ploidy, sex, is_wgd, q=q, min_read_diff=min_read_diff, rho_snv=rho_snv, deltaFreq=deltaFreq, xmin=xmin)
    return(data.frame(likelihood=likelihoods, aic=aics, bic=bics, binom_ll=binom_lls, binom_ll_2=binom_ll_2s, binom_ll_diff=binom_ll_diffs, mtimer_ll=mtimer_ll))
  }, mc.cores=MAXCORES)
  all_metrics = do.call(rbind, scores)
  return(all_metrics)
}

merge_superclones = function(cluster_locations, cluster_assignments, min_ccf_clone=0.95, max_ccf_clone=1.05) {
  to_remove = which(cluster_locations > max_ccf_clone)
  if (any(cluster_locations > min_ccf_clone & cluster_locations < max_ccf_clone)) {
    clonal_cluster = which(cluster_locations > min_ccf_clone & cluster_locations < max_ccf_clone)
  } else if (sum(cluster_locations > max_ccf_clone) > 1) {
    # no clone, but multiple superclones
    clonal_cluster = to_remove[which.min(abs(1-cluster_locations[cluster_locations > max_ccf_clone]))]
    to_remove = to_remove[!to_remove %in% clonal_cluster]
  } else {
    # no clone and the single superclone, that is the clone, leave it in place
    return(list(cluster_locations=cluster_locations, cluster_assignments=cluster_assignments))
  } 

  #' Pick highest CCF if multiple clonal candidates
  if (length(clonal_cluster) > 1) {
    clonal_cluster = clonal_cluster[which.max(cluster_locations[clonal_cluster])]
  }

  cluster_assignments[cluster_assignments %in% to_remove] = clonal_cluster
  cluster_locations = cluster_locations[-to_remove]
  
  #' Reindex the clusters - if needed
  cluster_ids = unique(cluster_assignments)
  if (max(cluster_ids, na.rm=T) != length(cluster_ids)) {
    new_assignments = cluster_assignments
    
    for (i in 1:length(unique(new_assignments))) {
      cluster_id = cluster_ids[i]
      new_assignments[cluster_assignments==cluster_id] = i
    }
    cluster_assignments = new_assignments
  }
  return(list(cluster_locations=cluster_locations, cluster_assignments=cluster_assignments))
}

randomclone_unif = function(dat, min_bound_data, max_bound_data, force_clone=F) {
  #' Draw number of clusters
  n_clusters = sample(MIN_CLUSTERS:MAX_CLUSTERS, 1)
  
  #' Get area where 95% of the data lives and set those boundaries as min and max CCF for clusters
  ccf_boundaries = quantile(dat$subclonal.fraction, probs=c(min_bound_data, max_bound_data))
  
  ccf_min = ccf_boundaries[[1]]
  ccf_max = ccf_boundaries[[2]]
  
  #' Draw clusters from uniform distribution
  cluster_locations = sort(runif(n_clusters, min=ccf_min, max=ccf_max), decreasing=T)
  if (force_clone) {
    cluster_locations[sample(1:n_clusters, 1)] = 1
  }
  
  #' Assign SNVs
  cluster_assignments = sapply(dat$subclonal.fraction, function(x) { which.min(abs(x-cluster_locations)) })
  
  #' Re-adjust the cluster locations based on assignments
  for (i in 1:n_clusters) {
    # Update only when there are mutations assigned
    if (sum(cluster_assignments==i, na.rm=T) > 0) {
      cluster_locations[i] = median(dat$subclonal.fraction[cluster_assignments==i], na.rm=T)
    }
  }
  
  #' Check if there is a superclone found and there is a clone - then merge
  if (any(cluster_locations > 0.90 & cluster_locations < 1.10) & any(cluster_locations > 1.10)) {
    res = merge_superclones(cluster_locations, cluster_assignments, min_ccf_clone=0.90, max_ccf_clone=1.10)
    cluster_locations = res$cluster_locations
    cluster_assignments = res$cluster_assignments
  }
  
  if (force_clone) {
    cluster_locations[which.min(abs(1-cluster_locations))] = 1
  }
  
  # #' Merge clusters too close
  # cluster_pair_distance = cluster_locations[1:(length(cluster_locations)-1)] - cluster_locations[2:length(cluster_locations)]
  # print(cluster_pair_distance)
  # if (any((cluster_pair_distance < 0.05) & cluster_locations[1:(length(cluster_locations)-1)] < 0.5)) {
  #   
  #   min_index = which.min(cluster_pair_distance)
  #   # merge cluster_locations[min_index] and cluster_locations[min_index+1]
  #   print(paste0("MERGING ", cluster_locations[min_index], " with ", cluster_locations[min_index+1]))
  #   
  #   
  #   cluster_assignments[cluster_assignments==min_index+1] = min_index
  #   cluster_pair_distance = cluster_pair_distance[-(min_index+1)]
  # }
  cluster_assignments = factor(cluster_assignments, levels=1:length(cluster_locations)) 
  structure_df = data.frame(table(cluster_assignments), 
                            proportion=cluster_locations * purity,
                            ccf=cluster_locations)
  colnames(structure_df)[1] = "cluster"
  colnames(structure_df)[2] = "n_ssms"
  
  return(list(structure=structure_df, assignments=cluster_assignments))
}

randomclone_stick = function(dat, force_clone=F) {
  #' Draw number of clusters
  n_clusters = sample(MIN_CLUSTERS:MAX_CLUSTERS, 1)
  
  #' Sort SNVs by their CCF
  snv_order = order(dat$subclonal.fraction)
  assignments = rep(NA, nrow(dat))
  if (n_clusters > 1) {
    #' Put n-1 breaks
    breaks = sort(sample(2:nrow(dat), n_clusters-1), decreasing=F)
    
    #' for the n clusters take median CCF as the cluster locations
    previous = 0
    for (i in 1:length(breaks)) {
      cluster_break = breaks[i]
      assignments[snv_order >= previous & snv_order < cluster_break] = i
      previous = cluster_break
    }

  } else {
    cluster_break = 0
    i = 0
  }
  
  #' Assign SNVs of the final cluster
  assignments[snv_order >= cluster_break] = i+1
  
  #' Determine cluster locations
  cluster_locations = rep(NA, n_clusters)
  for (i in 1:n_clusters) {
    cluster_locations[i] = median(dat$subclonal.fraction[assignments==i], na.rm=T)
  }
  cluster_locations = cluster_locations[!is.na(cluster_locations)]

  #' Check if there is a superclone found and there is a clone - then merge
  if (any(cluster_locations > 0.90 & cluster_locations < 1.10) & any(cluster_locations > 1.10)) {
    res = merge_superclones(cluster_locations, assignments, min_ccf_clone=0.90, max_ccf_clone=1.10)
    cluster_locations = res$cluster_locations
    assignments = res$cluster_assignments
  }
  
  #' Shift one cluster to be at exactly 1
  if (force_clone) {
    cluster_locations[which.min(abs(1-cluster_locations))] = 1
  }

  structure_df = data.frame(table(assignments), 
                            proportion=cluster_locations * purity,
                            ccf=cluster_locations)
  colnames(structure_df)[1] = "cluster"
  colnames(structure_df)[2] = "n_ssms"

  # remove emoty clusters
  structure_df = structure_df[structure_df$n_ssms > 0,]
  return(list(structure=structure_df, assignments=assignments))
}

calc_exp_mutreads_ccf = function(ccf, purity, ploidy, mean_depth) {
  return(ccf / (purity*ploidy + (1-purity)*2) * mean_depth)
}

mergeClustersByMutreadDiff = function(clusters, purity, ploidy, vcf_snv, min_read_diff) {
  clusters_new = clusters
  exp_reads = sapply(clusters$ccf, calc_exp_mutreads_ccf, purity=purity, ploidy=ploidy, mean_depth=mean(getTumorDepth(vcf_snv), na.rm=T))
  ccf_diff = exp_reads[1:(length(exp_reads)-1)] - exp_reads[2:length(exp_reads)]

  if (any(ccf_diff < min_read_diff)) {

    #' Iteratively merge a pair of clusters untill no more pairs within distance can be found
    merged = T
    while(merged) {
      merged = F

      exp_reads = sapply(clusters_new$ccf, calc_exp_mutreads_ccf, purity=purity, ploidy=ploidy, mean_depth=mean(getTumorDepth(vcf_snv), na.rm=T))
      ccf_diff = exp_reads[1:(length(exp_reads)-1)] - exp_reads[2:length(exp_reads)]
      to_merge = which(ccf_diff < min_read_diff)

      if (length(to_merge)==0) {
        merged = F
        break
      } else {
        i = to_merge[1]
        clusters_new$ccf[i] = sum(clusters_new$ccf[c(i, i+1)]*clusters_new$n_ssms[c(i, i+1)]) / sum(clusters_new$n_ssms[c(i, i+1)])
        clusters_new$n_ssms[i] = sum(clusters_new$n_ssms[c(i, i+1)])
        clusters_new = clusters_new[-(i+1),]
        merged = T
      }

      if (nrow(clusters_new)==1) {
        merged = F
        break
      }
    }
  }
  clusters_new$proportion = clusters_new$ccf * purity
  
  # sorting and renumbering clusters
  clusters_new = clusters_new[with(clusters_new, order(proportion, decreasing=T)),]
  clusters_new$cluster = 1:nrow(clusters_new)
  return(clusters_new)
}

assign_mtimer = function(MCN, clusters, purity) {
  best_cluster = sapply(MCN$D$CNF, function(x) if (is.na(x)) NA else which.min(abs(x-clusters$proportion)))
  cluster_counts = table(factor(best_cluster, levels=clusters$cluster))
  clusters_new_2 = data.frame(clusters$cluster, sapply(clusters$cluster, function(x) cluster_counts[[as.character(x)]]), clusters$proportion, clusters$ccf)
  colnames(clusters_new_2) = colnames(clusters)

  mcn = dpclust3p::mutationBurdenToMutationCopyNumber(burden=MCN$D$altCount / (MCN$D$altCount + MCN$D$wtCount), cellularity=purity, normalCopyNumber=rep(2, nrow(MCN$D)), totalCopyNumber=MCN$D$MajCN + MCN$D$MinCN)
  if (all(is.na(best_cluster))) {
    mcn = rep(NA, length(best_cluster))
  }
  ccf = MCN$D$CNF / purity
  plot_data_2 = data.frame(mcn=mcn, ccf=ccf, cluster=factor(best_cluster, levels=rev(unique(sort(clusters$cluster)))))
  return(list(plot_data=plot_data_2, clusters_new=clusters_new_2))
}

pcawg11_output = function(snv_mtimer, MCN, consensus_vcf_file) {
  # Cluster locations
  final_clusters = snv_mtimer$clusters

  # Assignments
  snv_assignments = data.frame(chr=as.character(seqnames(vcf_snv)), pos=as.numeric(start(vcf_snv)), cluster=snv_mtimer$plot_data$cluster)

  # Multiplicities
  snv_mult = data.frame(chr=snv_assignments$chr,
                        pos=snv_assignments$pos,
                        tumour_copynumber=MCN$D$MajCN+MCN$D$MinCN,
                        multiplicity=MCN$D$MutCN, multiplicity_options=NA, probabilities=NA)

  get_probs = function(final_clusters, MCN, vcf_snv) {
    n_subclones = nrow(final_clusters)-1
    if (n_subclones==0) {
      r = t(t(sapply(MCN$D$pAllSubclones, function(x) 0)))
    } else if (n_subclones==1) {
      #r = t(t(sapply(MCN$D$pAllSubclones, function(x) if(length(x)!=0) x else rep(NA, n_subclones))))
      r = matrix(unlist(sapply(MCN$D$pAllSubclones, function(x) if(length(x)!=0) x else rep(NA, n_subclones))))
    } else {
      # r = t(sapply(MCN$D$pAllSubclones, function(x) if(length(x)!=0) x else rep(1, n_subclones)))
      # r = matrix(unlist(sapply(MCN$D$pAllSubclones, function(x) if(length(x)!=0) x else rep(1, n_subclones))), ncol=n_subclones, byrow=T)
      r = matrix(unlist(lapply(MCN$D$pAllSubclones, function(x) if (is.null(x)) { rep(NA, n_subclones) } else { x })), ncol=n_subclones, byrow=T)
    }
    snv_assignments_prob = data.frame(chr=as.character(seqnames(vcf_snv)),
                                      pos=as.numeric(start(vcf_snv)),
                                      clone=1-rowSums(r),
                                      r, stringsAsFactors=F)

    if (n_subclones==0) {
      snv_assignments_prob = snv_assignments_prob[,1:3]
    }

    # set cluster number in the header
    colnames(snv_assignments_prob) = c("chr", "pos", paste0("cluster_", final_clusters$cluster))
    return(snv_assignments_prob)
  }

  # Obtain probabilities of assignments - snv
  snv_assignments_prob = get_probs(final_clusters, MCN, vcf_snv)

  # Recalculate the size of the clusters
  final_clusters$n_snvs = colSums(snv_assignments_prob[, grepl("cluster", colnames(snv_assignments_prob)), drop=F], na.rm=T)

  return(list(final_clusters=final_clusters,
        snv_assignments=snv_assignments,
        snv_mult=snv_mult,
        snv_assignments_prob=snv_assignments_prob))
}
