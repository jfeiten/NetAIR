#
# This is a Plumber API. You can run the API by clicking
# the 'Run API' button above.
#
# Find out more about building APIs with Plumber here:
#
#    https://www.rplumber.io/
#

library(plumber)
source("src/func.R")

#* @apiTitle Network Analysis API

#* Open index page
#* @get /start
#* @html
function(){
  require(readr)
  read_file("index.html")
}


#* Show work directory
#* @get /get
function(){
  getwd()
}

#* Read databases and do permutations between networks
#* @param path1 The path file data1
#* @param path2 The path file data2
#* @param it Number of permutations
#* @get /createnetworks
#* @html
function(file = "", desfecho_input = "", it = 2){
  
  require(readr)
  # copy .rmd ----
  folder_name <<- gsub(":", "", format(Sys.time(), "%X-%d%b%Y"))
  
  dir.create(paste0("cache/", folder_name))
  
  # identify the folders
  current.folder <- "src"
  new.folder <- paste0("cache/", folder_name)
  new.folder
  
  # find the files that you want
  list.of.files <- list.files(current.folder, "report_networks.Rmd", full = TRUE)
  list.of.files
  
  # copy the files to the new folder
  file.copy(list.of.files, new.folder)
  
  # create networks ----
  # define seed
  seed <- 2501
  set.seed(seed)
  n_perm <- as.numeric(it)
  
  
  
  path <- paste0("data/", file)
  desfecho <- desfecho_input
  
  # Baixa o banco sem missings
  dt_merged <- read.csv(path, header = TRUE, stringsAsFactors = FALSE)
  
  vars_type <- dt_merged[1, ]
  vars_type <- as.character(vars_type)
  vars_type <- vars_type[-which(colnames(dt_merged) == desfecho)]
  
  dt_merged <- dt_merged[-1, ]
  
  div_obj <- divideData(as.data.frame(dt_merged), desfecho = desfecho)
  
  data1 <- div_obj$data1
  data2 <- div_obj$data2
  
  set.seed(seed)
  permutations_results <- simulatePermutations(data1, data2, 10, vars_type)
  #permutations_results$near_zero
  
  m <- match(permutations_results$near_zero, colnames(data1))
  if (length(m) != 0) {
    data1 <- data1[, -m]
    data2 <- data2[, -m]
    vars_type <- vars_type[-m]
    
  }
  
  vars_levels <- getVarsLevels(data1)
  vars_levels
  
  data1 <- map_df(data1, function(x){as.numeric(as.character(x))})
  data1 <- as.matrix(data1)
  
  data2 <- map_df(data2, function(x){as.numeric(as.character(x))})
  data2 <- as.matrix(data2)
  
  netmeasures <- calculateNetworkMeasures(data1, data2, n_perm, vars_type, vars_levels)
  
  #save.image(file = paste0("cache/", folder_name, "/network_session.Rdata"))
  
  objs <- list(seed, n_perm, desfecho, dt_merged, vars_type, vars_levels, div_obj, 
               data1, data2, permutations_results, netmeasures)
  names(objs) <- list("seed", "n_perm", "desfecho", "dt_merged", "vars_type", "vars_levels", "div_obj", 
                      "data1", "data2", "permutations_results", "netmeasures")
  
  saveRDS(objs, file = paste0("cache/", folder_name, "/network_session.rds"))
  
  # read finished page ----
  read_file("finished_create_report.html")
}

#* Create results report
#* @get /createreport
#* @html
function(){
  require(readr)
  require(rmarkdown)
  
  file_path <- paste0("cache/", folder_name, "/report_networks.Rmd")
  rmarkdown::render(file_path, encoding = "UTF-8")
  
  file_path2 <- paste0("cache/", folder_name, "/report_networks.html")
  read_file(file_path2)
}

#* Open page create one network
#* @get /startone
#* @html
function(){
  require(readr)
  read_file("start_one.html")
}

#* Read databases and create one network
#* @param file Dataset file in .csv
#* @get /createonenetwork
#* @html
function(file = ""){
  
  require(readr)
  require(mgm)
  require(qgraph)
  
  # copy .rmd ----
  folder_name <<- gsub(":", "", format(Sys.time(), "%X-%d%b%Y"))
  
  dir.create(paste0("cache/", folder_name))
  
  # identify the folders
  current.folder <- "src"
  new.folder <- paste0("cache/", folder_name)
  new.folder
  
  # find the files that you want
  list.of.files <- list.files(current.folder, "report_one_network.Rmd", full = TRUE)
  list.of.files
  
  # copy the files to the new folder
  file.copy(list.of.files, new.folder)
  
  # create networks ----
  # define seed
  seed <- 2501
  set.seed(seed)
  
  path <- paste0("data/", file)
  
  # Baixa o banco sem missings
  dt_merged <- read.csv(path, header = TRUE, stringsAsFactors = FALSE)
  
  vars_type <- dt_merged[1, ]
  vars_type <- as.character(vars_type)
  
  vars_levels <- getVarsLevels(dt_merged)
  vars_levels
  
  dt_merged <- dt_merged[-1, ]
  dt <- map_df(dt_merged, function(x){as.numeric(as.character(x))})
  
  fit_mgm <- mgm(data = dt, type = vars_type, 
                  levels = vars_levels, k = 2, lambdaSel = "CV", 
                  lambdaFolds = 10, ruleReg = "AND")
  
  net <- qgraph(fit_mgm$pairwise$wadj, edge.color = fit_mgm$pairwise$edgecolor,
               nodeNames = colnames(dt), legend = TRUE, legend.cex = 0.1, label.cex = 3, vsize = 3, layout = "spring")
  
  
  objs <- list(seed, dt_merged, vars_type, vars_levels, fit_mgm, net)
  names(objs) <- list("seed", "dt_merged", "vars_type", "vars_levels", "fit_mgm", "net")
  
  saveRDS(objs, file = paste0("cache/", folder_name, "/one_network_session.rds"))
  
  # read finished page ----
  read_file("finished_create_one_network_report.html")
}

#* Create results report
#* @get /createonenetworkreport
#* @html
function(){
  require(readr)
  require(rmarkdown)
  
  file_path <- paste0("cache/", folder_name, "/report_one_network.Rmd")
  rmarkdown::render(file_path, encoding = "UTF-8")
  
  file_path2 <- paste0("cache/", folder_name, "/report_one_network.html")
  read_file(file_path2)
}
