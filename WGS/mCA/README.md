# cHL_mutation_calling/WGS/mCA
Pipeline adapted from https://github.com/tangdavid/mCAs_WGS to call mCAs in WGS data from cHL patients. The following steps are requied to run:

1. Recompile the binary from source using the updated Makefile (cHL_mutation_calling/WGS/mCA/Makefile) which has been successfully tested on H4H

2. Run the generate_phased_snps pipeline to create a per-sample bcf file which can be used as the input for any downstream chromosomal analysis

3. Run the cram_depth and call_mCAs pipelines as described on the original GitHub
