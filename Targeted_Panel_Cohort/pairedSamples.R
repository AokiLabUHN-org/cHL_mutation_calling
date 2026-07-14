library(stringr)
library(data.table)
mergedFunc = function(X) {
	collapsed = paste0(X, collapse = "_")
	return(gsub(" ", "", collapsed))
}

sample_pairs = as.data.frame(fread("/cluster/projects/aokigroup/Stueckmann/input/CHIP/sample_pairs.txt", header = F))
colnames(sample_pairs) = c("cHL", "BM")

global_DF = NULL
#loop through dataframe of paired samples
for (i in 1:nrow(sample_pairs) ) {
	BM = as.data.frame(fread(paste0("/cluster/projects/aokigroup/Stueckmann/out/sequenceData/mutations_cd/parsed/", sample_pairs[i, "BM"], ".parsed.tsv")))
	cHL = as.data.frame(fread(paste0("/cluster/projects/aokigroup/Stueckmann/out/sequenceData/mutations_cd/parsed/", sample_pairs[i, "cHL"], ".parsed.tsv")))
	#filter candidate somatic variants in blood on VAF to remove germline calls
	BM_filt = BM[which(BM$VAF > 0.02 & BM$VAF < 0.35),]
	muts_BM = apply(BM_filt[,1:6], 1, mergedFunc)
	rownames(BM_filt) = muts_BM
	#filter cHL variants to only those existing as somatic calls in matched blood
	muts_cHL = apply(cHL[,1:6], 1, mergedFunc)
	rownames(cHL) = muts_cHL
	muts_int = muts_cHL[which(muts_cHL %in% muts_BM)]

	#subset both BM and cHL dataframes to intersecting mutations
	cHL_filt = cHL[muts_int,]
	BM_filt = BM_filt[muts_int,]

	#filter based on consequence and evidence filter
	ix = which(cHL_filt$consequence %in% c("missense", "nonsense", "nonstop", "stop_lost", "frameshift", "inframe_indel", "splice_site", "start_codon") & cHL_filt$FILTER == "PASS")

	#see if any mutations remain
	if(length(ix) > 0) {
		tmpDF = data.frame(mutation = rownames(cHL_filt)[ix], cHL_vaf = cHL_filt$VAF[ix], BM_vaf = BM_filt$VAF[ix], cHL_sample = rep(sample_pairs[i, "cHL"], length(ix)), BM_sample = rep(sample_pairs[i, "BM"], length(ix)), cHL_support = cHL_filt$total_depth[ix], BM_support = BM_filt$total_depth[ix])
		if (is.null(global_DF)) {
			global_DF = tmpDF
		} else {
			global_DF = rbind(global_DF, tmpDF)
		}
	}


}

#mutations of interest
ix_m = which(global_DF$BM_support > 10 & global_DF$cHL_support > 10 & global_DF$cHL_vaf > global_DF$BM_vaf)
moi = global_DF[ix_m,]
saveRDS(moi, "/cluster/projects/aokigroup/Stueckmann/out/sequenceData/mutations_cd/paired_mutations_of_interest.rds")

#part two: analysis of paired HRS and ME samples
setwd("/cluster/projects/aokigroup/Stueckmann/out/sequenceData/mutations_cd/parsed/")
me_samples = list.files()[str_which(list.files(), "ME")]
paired_samples = gsub("ME", "", me_samples)

global_DF = NULL
for(i in 1:length(paired_samples)) {
	#only sample with single digit ID that does not prepend a zero
	if (paired_samples[i] == "cHL-LMD-3.parsed.tsv") {
		rhs_name = paired_samples[i]
	#fix incorrect ID for sample #50
	} else if (paired_samples[i] == "cHL-sort50-.parsed.tsv") {
		rhs_name = "cHL-sort-50.parsed.tsv"
	} else if(nchar(paired_samples[i]) == 20) {
		rhs_name = gsub("LMD-", "LMD-0", paired_samples[i])
	} else {
		rhs_name = paired_samples[i]
	}
	
	BM = as.data.frame(fread(me_samples[i]))
	cHL = as.data.frame(fread(rhs_name))

 	BM_filt = BM[which(BM$VAF > 0.02 & BM$VAF < 0.35),]
        muts_BM = apply(BM_filt[,1:6], 1, mergedFunc)
        rownames(BM_filt) = muts_BM
        #filter cHL variants to only those existing as somatic calls in matched blood
        muts_cHL = apply(cHL[,1:6], 1, mergedFunc)
        rownames(cHL) = muts_cHL
        muts_int = muts_cHL[which(muts_cHL %in% muts_BM)]

        #subset both BM and cHL dataframes to intersecting mutations
        cHL_filt = cHL[muts_int,]
        BM_filt = BM_filt[muts_int,]

        #filter based on consequence and evidence filter
        ix = which(cHL_filt$consequence %in% c("missense", "nonsense", "nonstop", "stop_lost", "frameshift", "inframe_indel", "splice_site", "start_codon") & cHL_filt$FILTER == "PASS")
 #see if any mutations remain
        if(length(ix) > 0) {
                tmpDF = data.frame(mutation = rownames(cHL_filt)[ix], cHL_vaf = cHL_filt$VAF[ix], BM_vaf = BM_filt$VAF[ix], cHL_sample = rep(rhs_name, length(ix)), BM_sample = rep(me_samples[i], length(ix)), cHL_support = cHL_filt$total_depth[ix], BM_support = BM_filt$total_depth[ix])
                if (is.null(global_DF)) {
                        global_DF = tmpDF
                } else {
                        global_DF = rbind(global_DF, tmpDF)
                }
        }

}
#mutations of interest
ix_m = which(global_DF$BM_support > 10 & global_DF$cHL_support > 10 & global_DF$cHL_vaf > global_DF$BM_vaf)
moi = global_DF[ix_m,]
saveRDS(moi, "/cluster/projects/aokigroup/Stueckmann/out/sequenceData/mutations_cd/paired_mutations_of_interest_ME.rds")



