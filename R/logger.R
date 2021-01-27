# https://www.datatechnotes.com/2017/09/logging-into-file-in-r-scripts.html

# Create log file name
get_log_file_name <- function(file_name="r_log") {
  file_name <- paste(file_name,format(Sys.time(), "%Y%m%d"),sep="_")
  file_name <- paste(log_path, paste(file_name,"log", sep="."), sep="/")
  return(file_name)
}

# logging function
logit <- function(msg, ...) {
   cat(format(Sys.time(), "%Y-%m-%d %X"), ":", paste(msg, ...), "\n",
    append = TRUE, file = logger)
}
