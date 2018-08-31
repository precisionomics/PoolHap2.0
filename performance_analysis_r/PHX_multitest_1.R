# Set the directory and get the data. 
setwd("C:/Users/Lauren Mak/Dropbox/University of Calgary/Project Pathogen WHEvolution/PoolHap_Testing/Pipeline_MultiTest")
stats_data = read.table(file = 'PHX_multitest_summary.txt', sep = '\t', header = TRUE,  as.is = numeric())

# GC, rjMCMC and PHX num_haps ~ ver_var_pos
moduleTypes = c("GC", "rjMCMC", "PHX_2")
for (m in 1:length(moduleTypes)) {
	mod = moduleTypes[m]
	tmp_haps = as.vector(t(stats_data[which(stats_data$Module==mod & stats_data$Summary_Stat=="num_haps"),-c(1,2)]))
	tmp = setNames(data.frame(t(stats_data[c(1,3,5),-c(1)])), t(stats_data[c(1,3,5),-c(1)])[1,])[-1,]
	tmp$num_haps = tmp_haps
	tmp_df1 = transform(tmp, orig_haps = as.numeric(as.character(orig_haps)), ver_var_pos = as.numeric(as.character(ver_var_pos)), prop_false = as.numeric(as.character(prop_false)))
	tmp_df = tmp_df1[(tmp_df1$num_haps != 0), ]	
	tmp_lm = lm(num_haps ~ ver_var_pos, data = tmp_df)
	print(mod)
	print(summary(tmp_lm))
}