rm(list=ls())
library(doSNOW)
library(foreach)
library(doParallel)

# Import gridsearch table
source("02a_simulation_grid_generation.R")
grid = make_grid()

# Defining function for cycling on grid (in parallel)
Gibbs <- function(i, grid){
    suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
    suppressWarnings(suppressPackageStartupMessages(library(ACutils)))
    suppressWarnings(suppressPackageStartupMessages(library(mvtnorm)))
    suppressWarnings(suppressPackageStartupMessages(library(salso)))
    suppressWarnings(suppressPackageStartupMessages(library(FGM)))
    suppressWarnings(suppressPackageStartupMessages(library(gmp)))
    suppressWarnings(suppressPackageStartupMessages(library(mcclust)))
    suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))
    suppressWarnings(suppressPackageStartupMessages(library(logr)))
    
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
    ## Generate data
    
    # Define true clustering
    rho_true = grid[i,]$rho_true # c(8,4,8,5)
    
    # Set seed for data generation
    seed_data_gen = grid[i,]$seed_data_gen
    
    # Define number of observations
    n = grid[i,]$n
    
    # Define variance of the Beta
    beta_sig2 = grid[i,]$beta_sig2
    
    # --------------------------------------------------------
    
    z_true = rho_to_z(rho_true)
    r_true = z_to_r(z_true)
    p = length(z_true)
    
    # Generate data
    sim = Generate_BlockDiagonal(n = n, z_true = z_true)
    
    if(grid[i,]$type_data_gen == "B"){
        sim = Generate_Block(
            n=n,
            z_true=z_true,
            p_block_diag = 1,
            p_block_extra_diag = 0,
            p_inside_block = 0.95,
            p_outside_block = 0.05,
            elem_out = 5,
            min_eigenval_correction = 3,
            seed = 1
        )
    }
    
    graph_density = sum(sim$Graph) / (p*(p-1))
    
    ## Set options for the simulation
    
    options = set_options(sigma_prior_0=0.5,
                          sigma_prior_parameters=list("a"=1,"b"=1,"c"=1,"d"=1),
                          theta_prior_0=1,
                          theta_prior_parameters=list("c"=1,"d"=1),
                          rho0=p, # start with one group
                          weights_a0=rep(1,p-1),
                          weights_d0=rep(1,p-1),
                          alpha_target=0.234,
                          beta_mu=graph_density,
                          beta_sig2=beta_sig2,
                          d=3,
                          alpha_add=0.5,
                          adaptation_step=1/(p*1000),
                          update_sigma_prior=TRUE,
                          update_theta_prior=TRUE,
                          update_weights=TRUE,
                          update_partition=TRUE,
                          update_graph=TRUE,
                          perform_shuffle=TRUE)
    
    dir.create(file.path("output", "data"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path("output", "log"), showWarnings = FALSE, recursive = TRUE)
    
    # unique_ID = uuid::UUIDgenerate(use.time = TRUE, n = 1, output = c("string"))
    # unique_ID = dittodb::hash(unique_ID, n = 8)
    # cat("This simulation has been assigned ID:", unique_ID)
    unique_ID = sprintf("%02d", i)
    grid[i,]$simulation_id = unique_ID
    
    filename_data = paste("output/data/simulation_", unique_ID, ".rds", sep = "")
    filename_log = paste("output/log/simulation_", unique_ID, ".log", sep = "")
    
    #log_open(file_name = filename_log, show_notes=FALSE, logdir = FALSE)
    res <- Gibbs_sampler(
        data = sim$data,
        niter = 8000,
        nburn = 2000,
        thin = 1,
        options = options,
        seed = 123456,
        print = TRUE
    )
    #log_close()
    
    res$true_rho = rho_true
    res$true_precision = sim$Prec
    res$true_graph = sim$Graph
    
    # save an object to a file
    saveRDS(res, file = filename_data)
    return(unique_ID)
}


#Setup backend to use many processors
totalCores = detectCores()
cluster <- makeCluster(totalCores[1])
registerDoParallel(cluster)
registerDoSNOW(cluster)
iterations <- nrow(grid)

# Progress bar
pb <- txtProgressBar(min = 1, max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Run for loop in Parallel
results <- foreach(i = 1:nrow(grid), .combine=rbind, .options.snow = opts) %dopar% {
    IDs <- Gibbs(i, grid)
    setTxtProgressBar(pb, i)
    return(IDs)
}

#Stop cluster
stopCluster(cluster)

# Append IDs to the grid
for(i in 1:nrow(grid)){
    grid[i,]$simulation_id <- results[i]
}

# Save the grid including IDs as .rds
saveRDS(grid, file = "output/simulation_table.rds")