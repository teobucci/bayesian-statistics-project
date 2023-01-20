## Source script to execute samplers on different datasets in parallel ##
## Once set, from R console type source("<path/to/this_file>") ##

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# The path to the script to execute in parallel
scriptpath = "./simulation_file.R"
SimNo <- 1:4

# Parallel execution (4 cores)
idx <- SimNo[1]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[2]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[3]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[4]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
