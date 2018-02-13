#
# This script merges two summary table entry output files from the randomclone informed procedure into one
#
args = commandArgs(T)
summary_table_mtimer_file = args[1]
summary_table_randomclone_file = args[2]
subclonal_structure_file = args[3]
samplename = args[4]
outdir = args[5]

summ = readr::read_tsv(summary_table_mtimer_file)
anno = readr::read_tsv(summary_table_randomclone_file)
unlink(c(summary_table_mtimer_file, summary_table_randomclone_file), recursive=F)

structure = readr::read_tsv(subclonal_structure_file)

anno$num_clonal = summ$num_clonal
anno$num_subclonal = summ$num_subclonal
anno$frac_clonal = summ$frac_clonal
anno$cluster_sizes = paste(structure$n_snvs, collapse=";")
anno = anno[, (!grepl("samplename.1", colnames(anno)))]
write.table(anno, file.path(outdir, paste0(samplename, "_summary_table.txt")), quote=F, sep="\t", row.names=F)
