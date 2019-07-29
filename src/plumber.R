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
  print(m)
  
  data1 <- data1[, -m]
  data2 <- data2[, -m]
  
  vars_type <- vars_type[-m]
  
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
  read_file("src/report_networks.html")
}


#* Plot a histogram
#* @png
#* @get /plot
function() {
    rand <- rnorm(100)
    hist(rand)
}

#* Return the sum of two numbers
#* @param a The first number to add
#* @param b The second number to add
#* @post /sum
function(a, b) {
    as.numeric(a) + as.numeric(b)
}
