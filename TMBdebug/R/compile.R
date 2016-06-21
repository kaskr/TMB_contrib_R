compile <-
function(...){
    args <- list(...)
    if (.Platform$OS.type == "windows") {
        args$flags <- "-O1 -g"
        args$DLLFLAGS <- ""
    } else {
        args$flags <- "-O0 -g"
    }
    do.call(TMB::compile, args)
}
