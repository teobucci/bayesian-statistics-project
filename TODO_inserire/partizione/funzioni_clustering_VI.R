# Old functions -----------------------------------------------------------

# Get labels for each iterations for each data point
#GDFMM$Partition is a (n_iter x n_data) matrix
part_matrix <- GDFMM$Partition
## ---> this is the object where you saved all the partitions


# Compute similarity matrix
sim_matrix <- psm(part_matrix)

# Heatmap, to understand the uncertainty in the similarity matrix
heatmap(sim_matrix)

# Compute final partition by minimazing the Binder loss or the VI loss functions
binder = minbinder(sim_matrix)
VI     = minVI(sim_matrix)

# Compare partitions just in terms of number of clusters and cluster cardinalities
table(binder$cl)
table(VI$cl)
table(real_partition)

# Compute Rand Index
arandi(binder$cl,real_partition, adjust = F)
arandi(VI$cl,real_partition, adjust = T)


# new functions -----------------------------------------------------------

library("Rcpp")
library("RcppArmadillo")
sourceCpp("wade.cpp")

# Get labels for each iterations for each data point
#GDFMM$Partition is a (n_iter x n_data) matrix
part_matrix <- GDFMM$Partition
## ---> this is the object where you saved all the partitions

# Compute similarity matrix
sim_matrix <- psm(part_matrix)

# Heatmap, to understand the uncertainty in the similarity matrix
heatmap(sim_matrix)

# compute VI loss for all visited partitions
dists <- VI_LB(part_matrix, psm_mat = sim_matrix)

# select best partition (among the visited ones)
final_partition_VI <- part_matrix[which.min(dists),]
table(final_partition_VI)

# Compute Rand Index
arandi(final_partition_VI,real_partition)
