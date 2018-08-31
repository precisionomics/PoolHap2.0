# Set the directory and get the data. 
setwd("C:/Users/Lauren Mak/Dropbox/University of Calgary/Project Pathogen WHEvolution/PoolHap_Testing/Pipeline_MultiTest")
stats_data = read.table(file = 'PHX_multitest_summary.txt', sep = '\t', header = TRUE,  as.is = numeric())

# GC, rjMCMC and PHX_* inter_freq_diff ~ ver_var_pos + num_haps + prop_error + prop_quasi.
# NOTE: Sometimes, inter_freq_diff is greater than 1. In that case, remove whole row. 
moduleTypes = c("GC", "rjMCMC", "PHX_0", "PHX_1", "PHX_2")
for (m in 1:length(moduleTypes)) {
	mod = moduleTypes[m]
	tmp_haps = as.vector(t(stats_data[which(stats_data$Module==mod & stats_data$Summary_Stat=="num_haps"),-c(1,2)]))
	tmp_err = as.vector(t(stats_data[which(stats_data$Module==mod & stats_data$Summary_Stat=="prop_error"),-c(1,2)]))	
	tmp_quasi = as.vector(t(stats_data[which(stats_data$Module==mod & stats_data$Summary_Stat=="prop_quasi"),-c(1,2)]))	
	tmp_ifd = as.vector(t(stats_data[which(stats_data$Module==mod & stats_data$Summary_Stat=="inter_freq_diff"),-c(1,2)]))
	tmp = setNames(data.frame(t(stats_data[3,-c(1)])), t(stats_data[3,-c(1)])[1,])
	tmp = data.frame(tmp[-1,])
	colnames(tmp) <- c("ver_var_pos")
	tmp$num_haps = tmp_haps
	tmp$prop_error = tmp_err
	tmp$prop_quasi = tmp_quasi
	tmp$inter_freq_diff = tmp_ifd
	tmp_df1 = transform(tmp, ver_var_pos = as.numeric(as.character(ver_var_pos)))
	tmp_df2 = tmp_df1[(tmp_df1$inter_freq_diff <= 1), ]
	tmp_df = tmp_df2[!(is.na(tmp_df2$inter_freq_diff)), ]
	tmp_lm = lm(inter_freq_diff ~ ver_var_pos + prop_error, data = tmp_df)
	print(mod)
	print(summary(tmp_lm))	
}
