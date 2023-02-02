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
suppressWarnings(suppressPackageStartupMessages(library(latex2exp)))

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

library("Rcpp")
library("RcppArmadillo")
sourceCpp("src/wade.cpp")


posterior_analysis <- function(i){

    simulation_id = grid[i,]$simulation_id
    dir_name = paste("simulation_", simulation_id, sep = "")
    output_path = file.path("output", "data", dir_name)
    dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
    
    sim <- read_rds(file.path("output","data",paste("simulation_",simulation_id,".rds",sep="")))
    
    # compute other partition forms
    rho_true = sim$true_rho
    r_true = rho_to_r(rho_true)
    z_true = rho_to_z(rho_true)
    p = length(z_true)
    num_clusters_true = length(rho_true)
    rho = sim$rho
    r = do.call(rbind, lapply(sim$rho, rho_to_r))
    z = do.call(rbind, lapply(sim$rho, rho_to_z))
    num_clusters = do.call(rbind, lapply(sim$rho, length))
    num_clusters = as.vector(num_clusters)

    # graph related quantities
    last_plinks = tail(sim$G, n=1)[[1]]
    bfdr_select = BFDR_selection(last_plinks, tol = seq(0.1, 1, by = 0.001))
    G_est = bfdr_select$best_truncated_graph # estimated adjacency
    
    # --------------------------------------------------------------------------
    
    # kl distance
    kl_dist = do.call(rbind, lapply(sim$K, function(k) {
        ACutils::KL_dist(sim$true_precision, k)
    }))
    
    # computing rand index for each iteration
    rand_index = apply(z, 1, mcclust::arandi, z_true)
    
    # compute VI
    sim_matrix <- salso::psm(z)
    dists <- VI_LB(z, psm_mat = sim_matrix)
    
    # select best partition (among the visited ones)
    best_partition_index = which.min(dists)
    rho_est = rho[[best_partition_index]]
    z_est = z[best_partition_index,]

    if(options$comparison_G_true_est){
        pdf(file=file.path(output_path,"comparison_G_true_est.pdf"), width = 16, height = 9)
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
            G_est,
            use_x11_device = F,
            horizontal = F,
            main = "\nEstimated Graph",
            center_value = NULL,
            col.upper = "black",
            col.center = "grey50",
            col.lower = "white"
        )
        dev.off()
    }
    
    if(options$comparison_K_true_est){
        pdf(file=file.path(output_path,"comparison_K_true_est.pdf"), width = 16, height = 9)
        par(mfrow=c(1,2))
        ACutils::ACheatmap(
            sim$true_precision,
            use_x11_device = F,
            horizontal = F,
            main = "Precision matrix",
            center_value = NULL,
            col.upper = "black",
            col.center = "grey50",
            col.lower = "white"
        )

        ACutils::ACheatmap(
            tail(sim$K,n=1)[[1]],
            use_x11_device = F,
            horizontal = F,
            main = "\nEstimated Precision matrix",
            center_value = NULL,
            col.upper = "black",
            col.center = "grey50",
            col.lower = "white"
        )
        dev.off()
    }
    
    if(options$plot_graph){
        par(mfrow=c(1,1))
        # create graph for visualization
        g1 <- graph.adjacency(sim$true_graph)
        edges1 <- get.edgelist(g1)
        edges1 <- cbind(edges1, rep("true", nrow(edges1)))
        g2 <- graph.adjacency(G_est)
        edges2 <- get.edgelist(g2)
        edges2 <- cbind(edges2, rep("estimated", nrow(edges2)))
        edges <- as.data.frame(rbind(edges1, edges2))
        names(edges) = c("from", "to", "graph")
        nodes = data.frame(
            vertices = 1:p,
            clust_true = as.factor(z_true),
            clust_est = as.factor(z_est)
        )
        # nodes
        g = graph_from_data_frame(edges, directed = FALSE, nodes)
        lay = create_layout(g, layout = "linear", circular = TRUE)

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
            cex.names=.6,
            las=2
        )
        
        abline(v=cp_true-0.5, col="red", lwd=2)
        legend("topright", legend=c("True"), col=c("red"),
               bty = "n",
               lty = 1,
               cex = 0.6)
        
        # --------------------------------------------

        plot(
            x = seq_along(rand_index),
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
        lines(x = seq_along(rand_index), y = rand_index)
        abline(h = 1, col = "red", lwd = 4)
        
        # --------------------------------------------
        
        plot(
            x = seq_along(kl_dist),
            y = kl_dist,
            type = "n",
            xlab = "Iterations",
            ylab = "K-L distance",
            main = paste(
                "\nKullback-Leibler distance\nLast:",
                round(tail(kl_dist, n=1), 3)
                )
        )
        lines(x = seq_along(kl_dist), y = kl_dist)
        dev.off()
    }
    

    if(options$theta_sigma){
        pdf(file=file.path(output_path,"theta_sigma_prior.pdf"), width = 16, height = 9)
        par(mfrow=c(2,1))
        plot(
            x = seq_along(sim$sigma),
            y = sim$sigma,
            type = "n",
            xlab = "Iterations",
            ylab = TeX(r'($\sigma$ prior)'),
            main = TeX(r'($\sigma$ prior - Traceplot)')
        )
        lines(x = seq_along(sim$sigma), y = sim$sigma, lwd = 0.8)
        
        # --------------------------------------------

        plot(
            x = seq_along(sim$theta),
            y = sim$theta,
            type = "n",
            xlab = "Iterations",
            ylab = TeX(r'($\theta$ prior)'),
            main = TeX(r'($\theta$ prior - Traceplot)')
        )
        lines(x = seq_along(sim$theta), y = sim$theta, lwd = 0.8)
        dev.off()

    }
    
    if(options$numgroups_frequency){
        pdf(file=file.path(output_path,"numgroups_frequency.pdf"), width = 16, height = 9)
        barplot(
            table(num_clusters),
            xlab = "Number of groups",
            ylab = "Frequency",
            main = paste(
                "Number of groups - Frequency\n",
                "Last:",
                tail(num_clusters, n = 1),
                "- Mean:",
                round(mean(num_clusters), 2),
                "- True:",
                num_clusters_true
            )
        )
        dev.off()
    }
    
    if(options$numgroups_traceplot){
        pdf(file=file.path(output_path,"numgroups_traceplot.pdf"), width = 16, height = 9)
        plot(
            x = seq_along(num_clusters),
            y = num_clusters,
            type = "n",
            xlab = "Iterations",
            ylab = "Number of groups",
            main = "Number of groups - Traceplot"
        )
        lines(x = seq_along(num_clusters), y = num_clusters)
        abline(h = length(z_to_rho(z_true)),
               col = "red",
               lwd = 4)
        legend("topleft", legend=c("True"), col=c("red"),
               lty = 1,
               cex = 1)
        dev.off()
    }
    
}

filename_data = "output/simulation_table.rds"
grid = readRDS(file = filename_data)
options <- list(comparison_G_true_est = TRUE,
                comparison_K_true_est = TRUE,
                plot_graph = TRUE,
                changepoint_kl = TRUE,
                numgroups_frequency = TRUE,
                numgroups_traceplot = FALSE,
                theta_sigma = TRUE)
pbsapply(1:nrow(grid), posterior_analysis)
