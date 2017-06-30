#' Method that does a simple random number of clusters and assigns SNVs to closest cluster
#' 
#' Rscript randomclone_unif.R 1e27cc8a-5394-4958-9af6-5ece1fe24516 1e27cc8a-5394-4958-9af6-5ece1fe24516_allDirichletProcessInfo.txt 0.77 GBM-US output/unif/
#' 
# set.seed(123) 
source("~/repo/randomclone/util.R")
MIN_CLUSTERS = 1
MAX_CLUSTERS = 5
ITERATIONS = 10
FORCE_CLONE = T
MIN_BOUND_DATA = .025
MAX_BOUND_DATA = .975

usemethod = "stick"
run_assessment = T

args = commandArgs(T)
samplename = args[1]
dpclustinput_infile = args[2]
purity = as.numeric(args[3])
project = args[4]
outdir = args[5]
if (length(args) > 5) {
  bb_file = args[6]
  vcf_snv = args[7]
  ploidy = as.numeric(args[8])
  sex = args[9]
  is_wgd = as.logical(args[10])
}

# samplename = "Sim_500_3"
# dpclustinput_infile = "test_data/dirichlet_input/Sim_500_003_allDirichletProcessInfo.txt"
# purity = 0.89
# project = "Bladder-TCC"

if (run_assessment) {
  # outdir = file.path("sandbox", samplename)
  outdir_bic = file.path(outdir, "bic")
  outdir_aic = file.path(outdir, "aic")
  outdir_binom = file.path(outdir, "binom")
  outdir_binom_2 = file.path(outdir, "binom_2")
  outdir_binom_diff = file.path(outdir, "binom_diff")
  outdir_ll = file.path(outdir, "ll")
  outdir_mtimer = file.path(outdir, "mtimer")
  if (!file.exists(outdir)) { dir.create(outdir) }
  if (!file.exists(outdir_bic)) { dir.create(outdir_bic) }
  if (!file.exists(outdir_binom)) { dir.create(outdir_binom) }
  if (!file.exists(outdir_binom_2)) { dir.create(outdir_binom_2) }
  if (!file.exists(outdir_binom_diff)) { dir.create(outdir_binom_diff) }
  if (!file.exists(outdir_aic)) { dir.create(outdir_aic) }
  if (!file.exists(outdir_ll)) { dir.create(outdir_ll) }
  if (!file.exists(outdir_mtimer)) { dir.create(outdir_mtimer) }
}

#' Load the data
dat = parse_data(dpclustinput_infile)
if (nrow(dat) < 2) {
  print("Not enough SNVs")
  q(save="no")
}

#' Run the method
res = list()
for (i in 1:ITERATIONS) {
  if (usemethod=="unif") {
    res[[i]] = randomclone_unif(dat, min_bound_data=MIN_BOUND_DATA, max_bound_data=MAX_BOUND_DATA, force_clone=FORCE_CLONE)
  } else if (usemethod=="stick") {
    res[[i]] = randomclone_stick(dat, force_clone=FORCE_CLONE)
  } else {
    print(paste0("Uknown method ", usemethod))
    q(save="no", status=1)
  }
}

# for (i in ITERATIONS:(ITERATIONS*2)) {
#   res[[i]] = randomclone_unif(dat, min_bound_data=MIN_BOUND_DATA, max_bound_data=MAX_BOUND_DATA, force_clone=!FORCE_CLONE)
# }

#' Calc overall likelihoods for every solution
save.image("testing.RData")
all_metrics2 = calc_all_metrics(dat, purity, res, vcf_snv, bb_file, ploidy, sex, is_wgd, q=0.05, min_read_diff=2, rho_snv=0.01, deltaFreq=0.00) 

if (run_assessment) {
  #' pick the best solution
  best_bic = which.min(all_metrics2$bic)
  best_aic = which.min(all_metrics2$aic)
  best_ll = which.min(all_metrics2$likelihood)
  best_binom = which.min(all_metrics2$binom_ll)
  best_binom_2 = which.min(all_metrics2$binom_ll_2)
  # best_binom_diff = which.min(all_metrics2$binom_ll_diff)
  best_mtimer = which.min(all_metrics2$mtimer_ll)
  
  #' Write the output
  save_output = function(res, best_index, outdir) {
    structure_df = res[[best_index]]$structure
    assignments = res[[best_index]]$assignments
    write_output_calibration_format(samplename, dat, structure_df, assignments, purity, outdir)
    write_output_summary_table(structure_df, outdir, samplename, project, purity)
  }
  save_output(res, best_ll, outdir_ll)
  save_output(res, best_bic, outdir_bic)
  save_output(res, best_aic, outdir_aic)
  save_output(res, best_binom, outdir_binom)
  save_output(res, best_binom_2, outdir_binom_2)
  save_output(res, best_mtimer, outdir_mtimer)
  # save_output(res, best_binom_diff, outdir_binom_diff)
} else {
  best_binom = which.min(all_metrics2$binom_ll)
}

if (run_assessment) {
  #' Make postprocess figure
  library(ggplot2)
  library(gridExtra)

  make_base_plot = function(dat) {
    p = ggplot() + geom_histogram(data=dat, mapping=aes(x=subclonal.fraction, y=..count..), fill="grey", colour="black", binwidth=0.025) +
      xlim(0, 1.5) + xlab("CCF") + ylab("Count")
    pb = ggplot_build(p)
    max_height = max(pb$data[[1]]$count)
    return(list(baseplot=p, max_height=max_height))
  }

  make_plot = function(cluster_locs, baseplot, plot_title, line_colour="red") {
    p = baseplot$baseplot
    max_height = baseplot$max_height
    p = p + geom_segment(data=cluster_locs, mapping=aes(x=ccf, xend=ccf, y=0, yend=max_height + 10), colour=line_colour)
    p = p + ggtitle(plot_title)
    return(p)
  }

  baseplot = make_base_plot(dat)

  truth = read.table(paste0("test_data/truth/Subclonal_Structure/", samplename, ".subclonal_structure.txt"), header=T, stringsAsFactors=F)
  truth_plot = make_plot(truth, baseplot, "truth", line_colour="green")

  aic_struct = read.table(file.path(outdir_aic, paste0(samplename, "_subclonal_structure.txt")), header=T, stringsAsFactors=F)
  aic_plot = make_plot(aic_struct, baseplot, "aic")

  bic_struct = read.table(file.path(outdir_bic, paste0(samplename, "_subclonal_structure.txt")), header=T, stringsAsFactors=F)
  bic_plot = make_plot(bic_struct, baseplot, "bic")

  ll_struct = read.table(file.path(outdir_ll, paste0(samplename, "_subclonal_structure.txt")), header=T, stringsAsFactors=F)
  ll_plot = make_plot(ll_struct, baseplot, "log-likelihood")

  binom_struct = read.table(file.path(outdir_binom, paste0(samplename, "_subclonal_structure.txt")), header=T, stringsAsFactors=F)
  binom_plot = make_plot(binom_struct, baseplot, "binomial")

  binom_2_struct = read.table(file.path(outdir_binom_2, paste0(samplename, "_subclonal_structure.txt")), header=T, stringsAsFactors=F)
  binom_2_plot = make_plot(binom_2_struct, baseplot, "binomial_2")
  
  mtimer_struct = read.table(file.path(outdir_mtimer, paste0(samplename, "_subclonal_structure.txt")), header=T, stringsAsFactors=F)
  mtimer_plot = make_plot(mtimer_struct, baseplot, "mtimer")

  # binom_diff_struct = read.table(file.path(outdir_binom_diff, paste0(samplename, "_subclonal_structure.txt")), header=T, stringsAsFactors=F)
  # binom_diff_plot = make_plot(binom_diff_struct, baseplot, "binomial_diff")

  png(file.path(outdir, paste0(samplename, "_solutions.png")), height=700, width=1000)
  grid.arrange(truth_plot, binom_plot, binom_2_plot, mtimer_plot, bic_plot, ll_plot, ncol=2,
               top=paste0(samplename, " - min_lust=", MIN_CLUSTERS, " max_lust=", MAX_CLUSTERS, " iters=", ITERATIONS, " force_clone=", FORCE_CLONE, " min_data=", MIN_BOUND_DATA, " max_data=", MAX_BOUND_DATA))
  dev.off()
}

#' Write the output - taking best_binom as best solution
write_output_calibration_format(samplename, dat, res[[best_binom]]$structure, res[[best_binom]]$assignments, purity, outdir)
write_output_summary_table(res[[best_binom]]$structure, outdir, samplename, project, purity)


