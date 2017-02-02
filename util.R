
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

