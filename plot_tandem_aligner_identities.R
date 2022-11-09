require(matrixStats)
require(ggplot2)
require(ggbeeswarm)

setwd("/Users/Jordan/Documents/Research/Pangenomics/Centromeres/working")

raw = read.table("active_array_identities.txt", header = T)

dat1 = raw[2:5]
dat2 = raw[6:9]

cols = c("parent1_haplotype1", "parent1_haplotype2", "parent2_haplotype1", "parent2_haplotype2")
colnames(dat1) = cols
colnames(dat2) = cols

dat = rbind(dat1, dat2)
best_match_val = rowMaxs(as.matrix(dat), na.rm = T)

flattened = unlist(dat)

plot_dat = data.frame(flattened, flattened == rep.int(best_match_val, 4))
colnames(plot_dat) = c("Identity", "Best_alignment_of_child")

plt = (ggplot(plot_dat, aes(y = Identity, x = 1, color = Best_alignment_of_child)) 
       + geom_beeswarm(priority = "density", cex = 1) 
       + scale_color_brewer(palette = "Set1")
       + theme(axis.text.x=element_blank(), 
               axis.ticks.x=element_blank(),
               axis.title.x=element_blank())
)
plt
