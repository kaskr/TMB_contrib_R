\name{MakeADFun}
\alias{MakeADFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Overloaded definition of MakeADFun to allow debugging on windows machines
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MakeADFun(..., DLL = TMB:::getUserDLL())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{DLL}{
%%     ~~Describe \code{DLL} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (..., DLL = TMB:::getUserDLL()) 
{
    orig_dir <- getwd()
    setwd(tempdir())
    on.exit(setwd(orig_dir))
    All_inputs <- list(..., DLL = DLL)
    save(All_inputs, file = "All_inputs.RData")
    DLL <- All_inputs$DLL
    DLLfull <- paste0(orig_dir, "/", DLL)
    txt <- c("library( TMB )", paste0("dyn.load(dynlib('", DLLfull, 
        "'))"), "load( 'All_inputs.RData' )", "Obj <- do.call(TMB::MakeADFun, All_inputs)")
    writeLines(txt, paste0(DLL, ".R"))
    Bdg_output <- gdbsource(paste0(DLL, ".R"))
    if (length(grep("#0", Bdg_output)) > 0) {
        message("Model has errors")
        print(Bdg_output)
        stop()
    }
    TMB::MakeADFun(..., DLL = DLL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
