# plot the size and chromosome of mutations seen in the pedigree
#
require(matrixStats)
require(ggplot2)
require(ggbeeswarm)
require(dplyr)
require(tidyr)
setwd("/Users/Jordan/Documents/Research/Pangenomics/Centromeres/working")

#tab_raw = read.table("mutations_in_child.tsv", header = T)
tab_raw = read.table("aggregated_mutations_in_child.tsv", header = T)

# flag some assemblies has having major errors in assembly or alignment
remove_chrom = data_frame(chrom = c("chr13", "chr7", "chr15"), hap = c(2, 1, 2))
keep_vec = rep.int(T, nrow(tab_raw))
for (i in 1:nrow(remove_chrom)) {
    keep_vec = keep_vec & !(tab_raw$chrom == remove_chrom$chrom[i] & tab_raw$child_hap == remove_chrom$hap[i])
}

tab = tab_raw[keep_vec,]

max_len = max(tab$length)
lims = c(0, max_len)
# lims = c(0, 15000)

plt = (tab %>% 
    select(Length=length, Chromosome=chrom, Mutation=type, Parent=parent) %>%
    ggplot(aes(y = Length, x = Chromosome, color = Mutation, shape = Parent)) + 
    geom_beeswarm(priority = "density", cex = .75, alpha = 1) + 
    ylim(lims)
)
plt
