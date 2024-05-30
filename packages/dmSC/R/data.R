# DATA
#
# Files for installing and managing data for DeMayo.
#
# Help for managing global variables within a package came from this article:
# https://www.r-bloggers.com/2017/09/global-variables-in-r-packages/

#' Set the search locations for dmSC datasets.
#' 
#' @param ... An array of paths to search for dmSC datasets.
#' 
#' @export
set_datapath <- function(...) {
  arguments <- unlist(list(...))
  unique_arguments <- unique(arguments)
  if (length(unique_arguments) < length(arguments)) {
    warning("Some duplicate entries were not added to the datapath.")
  }
  data_globals$path <- unique_arguments
}

#' Add search locations for dmSC datasets.
#' @export
append_datapath <- function(...) {
  set_datapath(data_globals$path, unlist(list(...)))
}

#' View search locations for dmSC datasets.
#' 
#' @return A vector describing the current search path.
#' 
#' @export
get_datapath <- function() {
  return(data_globals$path)
}

#' Create the data folder for holding dmSC datasets.
create_data_folder <- function(verbose = FALSE) {
  # Ensure data folder exists
  if (!dir.exists(data_globals$data_folder)) {
    if (verbose) {
      print("Data folder isn't created, creating it now.")
    }
    dir.create(path = data_globals$data_folder)
  }
}

#' Check if dmSC datasets are already installed. Not exported.
#' 
#' @return Whether a data name exists in the available search paths.
check_if_data_exists <- function(data_nm) {
  # Check if nm exists in data.path.
  installed_files <- list.files(path = data_globals$data_folder,
                                pattern = data_globals$valid_pattern,
                                ignore.case = FALSE,
                                recursive = FALSE)
  installed_files_no_ext <- tools::file_path_sans_ext(installed_files)
  return(data_nm %in% installed_files_no_ext)
}

#' List all installable dmSC datasets.
#' 
#' @return A vector of all installable file names. 
#'
#' @export
list_data <- function() {
  # List all data that are currently installed.
  # print("Note: Data versioning to come in the future.")
  print(c("Path currently includes: ", data_globals$path))
  all_files <- character(0)
  for (path in data_globals$path) {
    new_files <- list.files(path = path,
                            pattern = data_globals$valid_pattern,
                            recursive = FALSE)
    all_files <- c(all_files, new_files)
    new_indices <- (length(all_files) - length(new_files)):length(all_files)
    names(all_files)[new_indices] <- path
  }

  # Get rid of duplicates, accrue duplicate sources.
  all_files_unique <- split(seq_along(all_files), all_files)
  all_sources <- names(all_files)
  all_files <- names(all_files_unique)
  all_sources <- sapply(all_files_unique,
                        function(value) {
                          return(paste0(all_sources[value], collapse = ", "))
                        })

  # Get a list of files that are installed and return them.
  installed_files <- list.files(path = data_globals$data_folder,
                                pattern = data_globals$valid_pattern,
                                recursive = FALSE)
  install_queries <- sapply(all_files,
                            function(filename) {
                              return(filename %in% installed_files)
                            })
  file_information <- data.frame(file = tools::file_path_sans_ext(all_files),
                                 installed = install_queries,
                                 sources = all_sources)
  rownames(file_information) <- NULL

  return(file_information)
}

#' Install dmSC datasets from local locations or online locations.
#' 
#' @param nm Name of the dataset to install. Do not include extension.
#' @param ... Names of additional datasets to install. Optional.
#' @param reinstall Whether to re-copy the dataset to your installed
#'                  data files if the requested dataset is installed.
#'                  Default FALSE.
#' @param verbose Whether to log messages. Default FALSE.
#' 
#' @export
install_data <- function(nm, ..., reinstall = FALSE, verbose = FALSE) {
  # Create the data folder if not present.
  create_data_folder(verbose = verbose)

  # Announce the folder we are checking.
  if (verbose) {
    print(sprintf("Checking folder %s for existing data.",
                  data_globals$data_folder))
  }

  # For each of the datasets we are requested to install, install it to
  # the local repository.
  nms <- c(nm, as.character(...))
  for (nm in nms) {
    if (!reinstall) {
      if (check_if_data_exists(nm)) {
        if (verbose) {
          print(sprintf(paste0("Data '%s' already installed. Set 'reinstall = ",
                               "TRUE' if you wish to overwrite the existing ",
                               "file."), nm))
        }
        next
      }
    }

    # nm does not exist or we don't care if it's overwritten. Install it.
    found_nm <- FALSE
    for (path in data_globals$path) {
      is_url <- grepl("www.|http:|https:|ftp:", path)
      if (!is_url) {
        searched_files <- list.files(path = path,
                                     pattern = data_globals$valid_pattern,
                                     ignore.case = FALSE,
                                     recursive = FALSE)
        searched_files_no_ext <- tools::file_path_sans_ext(searched_files)
        if (nm %in% searched_files_no_ext) {
          found_nm <- TRUE
          tgt_file <- searched_files[searched_files_no_ext == nm][1]
          src <- file.path(path, tgt_file)
          dest <- data_globals$data_folder
          if (verbose) {
            print(sprintf("Found local data '%s', copying from %s to %s...",
                          nm, src, dest))
          }
          file.copy(from = src, to = dest)
          break
        }
      } else {
        warning(sprintf(paste0("Online dmSC data_path entry '%s' is not yet ",
                               "supported. Skipping."), path))
      }
    }
    if (!found_nm) {
      stop(sprintf("Could not find data '%s' on dmSC's data_path.", nm))
    }
  }
}

# Make the new environment and make the data folder if not present.
data_globals <- new.env()
data_globals$data_folder <- "./data"
create_data_folder(verbose = FALSE)
data_globals$data_folder <- normalizePath("./data")
data_globals$path <- character()
data_globals$valid_pattern <- paste("^.*\\.R$|^.*\\.r$|^.*\\.rda$|",
                                    "^.*\\.RData$|^.*\\.txt$|^.*\\.TXT$|",
                                    "^.*\\.tab$|^.*csv$|^.*\\.CSV$",
                                    collapse = "")