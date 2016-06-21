MakeADFun <-
function(..., DLL = TMB:::getUserDLL()){
    ## Set local working directory
    orig_dir <- getwd()
    setwd( tempdir() )
    on.exit( setwd(orig_dir) )
    ## Save inputs
    All_inputs <- list(..., DLL=DLL)
    save( All_inputs, file="All_inputs.RData")
    ## Copy DLL to tempdir
    DLL <- All_inputs$DLL
    DLLfull <- paste0(orig_dir,"/",DLL)
    ## Write file to source
    txt <- c("library( TMB )",
             paste0("dyn.load(dynlib('",DLLfull,"'))"),
             "load( 'All_inputs.RData' )",
             "Obj <- do.call(TMB::MakeADFun, All_inputs)"
             )
    writeLines(txt, paste0(DLL,".R"))
    ## Try running
    Bdg_output <- gdbsource(paste0(DLL, ".R"))
    # Sort out outcomes
    if( length(grep("#0",Bdg_output)) > 0 ){
        message("Model has errors")
        print(Bdg_output)
        stop()
    }
    ## OK ==> Safe to run in main session:
    TMB::MakeADFun(...,DLL=DLL)
}
