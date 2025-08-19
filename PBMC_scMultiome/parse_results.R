data_tab <- read.csv("/gstore/project/epigen/benchmark/resource_use.txt")
get_file_paths <- function(df) paste0(df[["path"]],"/",rep("slurm-",nrow(df)),df[["job_id"]],rep(".out", nrow(df)))

output_files <- get_file_paths(data_tab)

get_exit_status <- function(file_path){
    if(!file.exists(file_path)) return(1)
    connection <- file(file_path)
    text <- readLines(connection)
    if(length(text)<3) return(1)
    if(text[length(text)-2]=="Execution halted") return(1)
    if(grepl("error", text[length(text)-2])) return(1)
    if(grepl("error", text[length(text)])) return(1)
    if(grepl("command not found", text[length(text)-2])) return(1)
    if(grepl("Nothing to be done", text[length(text)-2])) return(1)
    if(grepl("No such file or directory", text[length(text)])) return(1)
    last_but_one_line = text[length(text)-1]
    if (grepl("JobId",last_but_one_line)) {
        exit_code = gsub(".*ExitCode=(.{1,3}) .*","\\1",last_but_one_line)
        if(exit_code=="0:0") return(0)
    }
    return(1)
}

output_files_zero_status <- output_files[!as.logical(sapply(output_files, get_exit_status))]

get_measurement_data <- function(path){
    connection <- file(path)
    text <- readLines(connection)
    last_line = text[length(text)]
    job_name = gsub(".*[0-9]{8}\\|(.*?)\\|.*", "\\1", last_line)
    node_ID = gsub(".*NodeList=(.*) BatchHost.*", "\\1", text[length(text)-1])
    memory = gsub(".*,mem=([0-9]{1,4}G),node.*", "\\1", text[length(text)-1])
    new_job_name = gsub(" ", "_", job_name)
    last_line = gsub(job_name, new_job_name, last_line)
    entries <- strsplit(last_line, "\\|")[[1]]
    entries <- unlist(lapply(entries, function(x) {if(x=="") return("-"); return(x)}))
    entries <- unlist(lapply(entries, function(x) strsplit(x, " ")))
    colnames <- entries[1:12]
    # row_data <- setNames(entries[13:24], colnames)
    # df <- data.frame()
    # df <- rbind(df, data.frame(as.list(c(row_data, node = node_ID))))
    # df <- rbind(df, data.frame(as.list(c(setNames(entries[25:36], colnames), node = node_ID))))
    # df <- rbind(df, data.frame(as.list(c(setNames(entries[37:48], colnames), node = node_ID))))
    # df
      data.frame(as.list(c(setNames(c(entries[13:14], entries[27:36]), colnames), node = node_ID, requested_memory = memory)))
}

dfs <- lapply(output_files_zero_status, get_measurement_data)
df <- do.call(rbind, dfs)
df <- df[,c("JobID","JobName","requested_memory","AllocCPUS","Elapsed","AveCPU","CPUTime","MaxVMSizeNode","MaxVMSize","AveVMSize","ConsumedEnergy","MaxRSS","ExitCode","node")]


