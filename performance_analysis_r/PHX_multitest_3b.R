# Set the directory and get the data. 
setwd("C:/Users/Lauren Mak/Dropbox/University of Calgary/Project Pathogen WHEvolution/PoolHap_Testing/Pipeline_MultiTest")
stats_data = read.table(file = 'PHX_multitest_summary.txt', sep = '\t', header = TRUE,  as.is = numeric())

# GC, rjMCMC and PHX prop_quasi ~ ver_var_pos + num_haps + prop_error
moduleTypes = c("PHX_0", "PHX_1", "PHX_2")
tmp_all = data.frame()
for (m in 1:length(moduleTypes)) {
	mod = moduleTypes[m]
	tmp_haps = as.vector(t(stats_data[which(stats_data$Module==mod & stats_data$Summary_Stat=="num_haps"),-c(1,2)]))
	tmp_md = as.vector(t(stats_data[which(stats_data$Module==mod & stats_data$Summary_Stat=="max_diff"),-c(1,2)]))	
	tmp_err = as.vector(t(stats_data[which(stats_data$Module==mod & stats_data$Summary_Stat=="prop_error"),-c(1,2)]))	
	tmp_quasi = as.vector(t(stats_data[which(stats_data$Module==mod & stats_data$Summary_Stat=="prop_quasi"),-c(1,2)]))	
	tmp = setNames(data.frame(t(stats_data[3,-c(1)])), t(stats_data[3,-c(1)])[1,])
	tmp = data.frame(tmp[-1,])
	colnames(tmp) <- c("ver_var_pos")
	tmp$num_haps = tmp_haps
	tmp$max_diff = tmp_md
	tmp$prop_error = tmp_err
	tmp$prop_quasi = tmp_quasi
	tmp_df = transform(tmp, ver_var_pos = as.numeric(as.character(ver_var_pos)))
	tmp_all = rbind(tmp_all,tmp_df)
}
tmp_lm = lm(prop_quasi ~ max_diff + prop_error, data = tmp_all)
print(summary(tmp_lm))