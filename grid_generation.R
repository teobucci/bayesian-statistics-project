##### CREATING DIFFERENT GRIDS FOR TESTING #####
rho_true <- c(8,4,8,5)

#### GRID 1: VARYING SEED ####
l <- 8 # Grid Length
simulation_id <- rep("_",l)
n <- rep(500,l)
p <- rep(25,l)
type_data_gen <- rep("BD",l)
seed_data_gen <-  c(22111996,
                    31051999,
                    27051999,
                    29061999,
                    12091997,
                    27091999,
                    27121996)
beta_sig2 <- rep(0.20,l)
rho_true_vec <- rep(list(rho_true),l)
initial_partition <- rep(list(25),l)
grid1 <- cbind(simulation_id = simulation_id, 
               n = n, p = p, 
               type_data_gen = type_data_gen, 
               seed_data_gen = seed_data_gen,
               initial_partition = initial_partition,
               beta_sig2 = beta_sig2, 
               rho_true = rho_true_vec)



#### GRID 2: VARYING BETA_SIG_2 ####
l <- 10 # Grid Length
simulation_id <- rep("_",l)
n <- rep(500,l)
p <- rep(25,l)
type_data_gen <- rep("BD",l)
seed_data_gen <- rep(27121996,l)
beta_sig2 <- seq(1/16,0.2,length.out=l) # Here there should be 0.25 instead of 0.2
rho_true_vec <- rep(list(rho_true),l)
initial_partition = rep(list(25),l)
grid2 <- cbind(simulation_id = simulation_id, 
               n = n, p = p, 
               type_data_gen = type_data_gen, 
               seed_data_gen = seed_data_gen,
               initial_partition = initial_partition,
               beta_sig2 = beta_sig2, 
               rho_true = rho_true_vec)



#### GRID 3: CHANGE INITIAL PARTITION ####
l <- 2 # Grid Length
simulation_id <- rep("_",l)
n <- rep(500,l)
p <- rep(25,l)
type_data_gen <- rep("BD",l)
seed_data_gen <- rep(27121996,l)
beta_sig2 <- rep(0.20,l)
rho_true_vec <- rep(list(rho_true),l)
initial_partition = list(25, c(rep(1,25)))
grid3 <- cbind(simulation_id = simulation_id, 
               n = n, p = p, 
               type_data_gen = type_data_gen, 
               seed_data_gen = seed_data_gen,
               initial_partition = initial_partition,
               beta_sig2 = beta_sig2, 
               rho_true = rho_true_vec)



#### GRID 4: ROBUSTNESS WRT TO TRUE PARTITIONS ####
l <- 4# Grid Length
simulation_id <- rep("_",l)
n <- rep(500,l)
p <- rep(25,l)
type_data_gen <- rep("BD",l)
seed_data_gen <- rep(27121996,l)
beta_sig2 <- rep(0.20,l)
beta_sig2[3] <- 0.02 
rho_true_vec <- list(c(8,4,8,5),
                   c(1,10,2,9,3),
                   c(1,3,2,4,2,3,3,4,3),
                   c(12,13))
initial_partition = list(25, 
                        c(12,13), 
                        c(8,8,9), 
                        c(6,6,6,7))
grid4 <- cbind(simulation_id = simulation_id, 
               n = n, p = p, 
               type_data_gen = type_data_gen, 
               seed_data_gen = seed_data_gen,
               initial_partition = initial_partition,
               beta_sig2 = beta_sig2, 
               rho_true = rho_true_vec)



#### GRID 5: VARYING N (NUMBER OF OBSERVATIONS) ####
l <- 7 # Grid Length
simulation_id <- rep("_",l)
n <- c(500,400,300,200,100,50,20)
p <- rep(25,l)
type_data_gen <- rep("BD",l)
seed_data_gen <-  rep(27121996, l)
beta_sig2 <- rep(0.10,l)
rho_true_vec <- rep(list(rho_true),l)
initial_partition <- rep(list(25),l)
grid5 <- cbind(simulation_id = simulation_id, 
               n = n, p = p, 
               type_data_gen = type_data_gen, 
               seed_data_gen = seed_data_gen,
               initial_partition = initial_partition,
               beta_sig2 = beta_sig2, 
               rho_true = rho_true_vec)


#### GRID 6: VARYING P ####
l <- 10 # Grid Length
simulation_id <- rep("_",l)
n <- rep(500,l)
p <- c(5,10,15,20,25,30,35,40,45,50)
type_data_gen <- rep("BD",l)
seed_data_gen <-  rep(27121996, l)
beta_sig2 <- rep(0.0625,l)
rho_true_vec <- c(
    list(c(3,2)),
    list(c(5,5)),
    list(c(5,5,5)),
    list(c(5,5,5,5)),
    list(c(5,5,5,5,5)),
    list(c(5,5,5,5,5,5)),
    list(c(5,5,5,5,5,5,5)),
    list(c(5,5,5,5,5,5,5,5)),
    list(c(5,5,5,5,5,5,5,5,5)),
    list(c(5,5,5,5,5,5,5,5,5,5))
)

initial_partition <- sapply(p,list)

grid6 <- cbind(simulation_id = simulation_id, 
               n = n, p = p, 
               type_data_gen = type_data_gen, 
               seed_data_gen = seed_data_gen,
               initial_partition = initial_partition,
               beta_sig2 = beta_sig2, 
               rho_true = rho_true_vec)


#### GRID 7: NOISE ON THE GRAPH ####
l <- 1 # Grid Length
simulation_id <- rep("_",l)
n <- rep(500,l)
p <- rep(25,l)
type_data_gen <- rep("B",l)
seed_data_gen <-  rep(27121996, l)
beta_sig2 <- rep(0.2,l)
rho_true_vec <- rep(list(rho_true),l)
initial_partition <- rep(list(25),l)

grid7 <- cbind(simulation_id = simulation_id, 
               n = n, p = p, 
               type_data_gen = type_data_gen, 
               seed_data_gen = seed_data_gen,
               initial_partition = initial_partition,
               beta_sig2 = beta_sig2, 
               rho_true = rho_true_vec)

###########################################
# Creating final grid

grid <- rbind(grid1,grid2,grid3,grid4,grid5,grid6,grid7)
rm(grid1)
rm(grid2)
rm(grid3)
rm(grid4)
rm(grid5)
rm(grid6)
rm(grid7)

