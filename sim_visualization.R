suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(ACutils)))
suppressWarnings(suppressPackageStartupMessages(library(mvtnorm)))
suppressWarnings(suppressPackageStartupMessages(library(salso)))
suppressWarnings(suppressPackageStartupMessages(library(FGM)))
suppressWarnings(suppressPackageStartupMessages(library(gmp)))
suppressWarnings(suppressPackageStartupMessages(library(mcclust)))
suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))
suppressWarnings(suppressPackageStartupMessages(library(logr)))
suppressWarnings(suppressPackageStartupMessages(library(tidygraph)))
suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(pbapply)))
##### options(warn=1)

paths = c(
    "src/utility_functions.R",
    "src/bulky_functions.R",
    "src/data_generation.R"
)

for(p in paths){
    path = file.path(p)
    if(file.exists(path)){
        source(path)
    } else {
        cat("File",path,"was not found in directory, please check.")
    }
}



posterior_analysis <- function(i){
    library("Rcpp")
    library("RcppArmadillo")
    sourceCpp("src/wade.cpp")
    
    simulation_id = grid[i,]$simulation_id
    dir_name = paste("simulation_", simulation_id, sep = "")
    output_path = file.path("output", "data", dir_name)
    dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
    
    sim <- read_rds(file.path("output","data",paste("simulation_",simulation_id,".rds",sep="")))
    rho_true <- grid[i,]$rho_true
    z_true = rho_to_z(rho_true)
    r_true = z_to_r(z_true)
    p = length(z_true)
    num_clusters_true = length(rho_true)
    
    Prec_true = sim$true_precision # precision matrix
    Graph_true = sim$true_graph # true graph
    Graph_true[col(Graph_true)==row(Graph_true)] = 0 # remove self loops
    graph_density = sum(sim$Graph) / (p*(p-1))
    
    r = do.call(rbind, lapply(sim$rho, rho_to_r))
    z = do.call(rbind, lapply(sim$rho, rho_to_z))
    r_true = rho_to_r(sim$true_rho)
    z_true = rho_to_z(sim$true_rho)
    num_clusters = do.call(rbind, lapply(sim$rho, length))
    num_clusters = as.vector(num_clusters)
    
    mean(sim$accepted) # acceptance frequency
    # computing rand index for each iteration
    rand_index = apply(z, 1, mcclust::arandi, z_true)
    
    # compute VI loss for all visited partitions
    # compute similarity matrix 
    sim_matrix <- salso::psm(z)
    # adding names for the heatmap
    rownames(sim_matrix) = 1:length(z_true)
    colnames(sim_matrix) = 1:length(z_true)
    dists <- VI_LB(z, psm_mat = sim_matrix)
    
    # select best partition (among the visited ones)
    final_partition_VI <- z[which.min(dists),]
    unname(table(final_partition_VI))
    
    # compute Rand Index
    mcclust::arandi(final_partition_VI,z_true)
    
    # kl distance
    kl_dist = do.call(rbind, lapply(sim$K, function(k) {
        ACutils::KL_dist(sim$true_precision, k)
    }))
    last = round(tail(kl_dist, n=1), 3)
    
    
    last_G = sim$G[[length(sim$G)]]
    bfdr_select = BFDR_selection(last_G, tol = seq(0.1, 1, by = 0.001))
    bfdr_select$best_treshold
    final_graph = bfdr_select$best_truncated_graph
    
    # create graph for visualization
    g1 <- graph.adjacency(Graph_true)
    edges1 <- get.edgelist(g1)
    edges1 <- cbind(edges1,rep("true",nrow(edges1)))
    g2 <- graph.adjacency(final_graph)
    edges2 <- get.edgelist(g2)
    edges2 <- cbind(edges2,rep("estimated",nrow(edges2)))
    edges <- as.data.frame(rbind(edges1,edges2))
    names(edges) = c("from","to","graph")
    nodes = data.frame(vertices=1:p,true_clust=as.factor(z_true), estim_clust=as.factor(final_partition_VI))
    nodes
    g = graph_from_data_frame(edges, directed = FALSE, nodes)
    lay = create_layout(g, layout = "fr")
    
    
    if(options$true_graph_and_K){
        pdf(file=file.path(output_path,"true_G_K.pdf"), width = 16)
        par(mfrow=c(1,2))
        ACutils::ACheatmap(
            sim$true_graph,
            use_x11_device = F,
            horizontal = F,
            main = "Graph",
            center_value = NULL,
            col.upper = "black",
            col.center = "grey50",
            col.lower = "white"
        )
        
        ACutils::ACheatmap(
            sim$true_precision,
            use_x11_device = F,
            horizontal = F,
            main = "Precision matrix K",
            center_value = NULL,
            col.upper = "black",
            col.center = "grey50",
            col.lower = "white"
        )
        dev.off()
    }
    
    
    if(options$estimated_graph_and_K){
        pdf(file=file.path(output_path,"estimated_G_K.pdf"), width = 16)
        par(mfrow=c(1,2))
        ACutils::ACheatmap(
            tail(sim$G,n=1)[[1]],
            use_x11_device = F,
            horizontal = F,
            main = "\nEstimated plinks matrix",
            center_value = NULL,
            col.upper = "black",
            col.center = "grey50",
            col.lower = "white"
        )
        
        ACutils::ACheatmap(
            tail(sim$K,n=1)[[1]],
            use_x11_device = F,
            horizontal = F,
            main = "\nEstimated Precision matrix K",
            center_value = NULL,
            col.upper = "black",
            col.center = "grey50",
            col.lower = "white"
        )
        dev.off()
        
    }
    
    if(options$plot_graph){
        par(mfrow=c(1,1))
        ggraph(lay) + 
            geom_edge_link(edge_colour = "grey") + 
            geom_node_point(aes(color=true_clust, shape = estim_clust), size = 4) +
            geom_node_text(aes(label = name), repel=TRUE) +
            facet_edges(~graph)
        ggsave(file.path(output_path,"graph.pdf"), device = "pdf", width=12.2, height=8.21)
    }
    
    
    if(options$changepoint_kl){
        pdf(file=file.path(output_path,"indexes.pdf"), width = 16, height = 9)
        par(mfrow=c(1,3))
        ### Barplot of changepoints
        bar_heights = colSums(r)
        cp_true = which(r_true==1)
        color <- ifelse(seq_along(bar_heights) %in% c(cp_true), "red", "gray")
        
        barplot(
            bar_heights,
            names = seq_along(bar_heights),
            border = "NA",
            space = 0,
            yaxt = "n",
            main = "\nChangepoint\n frequency distribution",
            #col = color,
            cex.names=.6,
            las=2
        )
        
        abline(v=cp_true-0.5, col="red", lwd=2)
        legend("topright", legend=c("True"), col=c("red"),
               bty = "n",
               lty = 1,
               cex = 0.6)
        #title(paste("Acceptance frequency:",round(mean(sim$accepted),2)), line = -20)
        
        plot(
            x = 1:length(rand_index),
            y = rand_index,
            type = "n",
            xlab = "Iterations",
            ylab = "Rand Index",
            main = paste(
                "\nRand Index - Traceplot\n",
                "Last:",
                round(tail(rand_index, n=1), 3),
                "- Mean:",
                round(mean(rand_index), 2)
            )
        )
        lines(x = 1:length(rand_index), y = rand_index)
        abline(h = 1, col = "red", lwd = 4)
        
        
        plot(
            x = 1:length(kl_dist),
            y = kl_dist,
            type = "n",
            xlab = "Iterations",
            ylab = "K-L distance",
            main = paste("\nKullback Leibler distance\nLast value:", last)
        )
        lines(x = 1:length(kl_dist), y = kl_dist)
        dev.off()
        
    }
    

    if(options$theta_sigma){
        pdf(file=file.path(output_path,"theta_sigma_prior.pdf"), width = 16, height = 9)
        par(mfrow=c(2,1))
        plot(
            x = 1:length(sim$sigma),
            y = sim$sigma,
            type = "n",
            xlab = "Iterations",
            ylab = TeX(r'($\sigma$ prior)'),
            main = TeX(r'($\sigma$ prior - Traceplot)')
        )
        lines(x = 1:length(sim$sigma), y = sim$sigma, lwd = 0.8)
        
        plot(
            x = 1:length(sim$theta),
            y = sim$theta,
            type = "n",
            xlab = "Iterations",
            ylab = TeX(r'($\theta$ prior)'),
            main = TeX(r'($\theta$ prior - Traceplot)')
        )
        lines(x = 1:length(sim$theta), y = sim$theta, lwd = 0.8)
        dev.off()

    }
    
    
    
    if(options$groups_number){
        par(mfrow=c(1,2))
        histogram = hist(
            num_clusters,
            xlab = "Number of groups",
            ylab = "Frequency",
            border = "darkgray",
            main = paste(
                "Number of groups - Frequency\n",
                "Last:",
                tail(num_clusters, n = 1),
                "- Mean:",
                round(mean(num_clusters), 2),
                "- True:",
                num_clusters_true
            ),
            xaxt = "n"
        )
        
        axis(
            side = 1,
            at = histogram$mids,
            labels = seq(from = floor(histogram$mids[1]),
                         length.out = length(histogram$mids))
        ) 
        
        abline(v = num_clusters_true+0.25,
               col = "red",
               lwd = 4)
        
        
        
        plot(
            x = 1:length(num_clusters),
            y = num_clusters,
            type = "n",
            xlab = "Iterations",
            ylab = "Number of groups",
            main = "Number of groups - Traceplot"
        )
        lines(x = 1:length(num_clusters), y = num_clusters)
        abline(h = length(z_to_rho(z_true)),
               col = "red",
               lwd = 4)
        legend("topleft", legend=c("True"), col=c("red"),
               lty = 1,
               cex = 1)
    }
    
}
filename_data = "output/simulation_table.rds"
grid = readRDS(file = filename_data)
options <- list(true_graph_and_K = TRUE,
                estimated_graph_and_K = TRUE,
                plot_graph = TRUE,
                changepoint_kl = TRUE,
                groups_number = FALSE,
                theta_sigma = TRUE)
pbsapply(1:nrow(grid), posterior_analysis)
