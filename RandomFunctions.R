###### Loading libraries ########
# List of packages for session
.packages = c("plm", 
              "lmtest",
              "lfe",
              "multiwayvcov",
              "sandwich",
              "DataExplorer",
              "stringr",
              "devtools")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only = TRUE)


#### Options I Like ####
options(dplyr.width = Inf) # Makes dplyr show all columns in commands like head
options(stringsAsFactors = FALSE) # loads string variables as strings instead of factors


# Version of stargazer with booktabs
install_github("markwestcott34/stargazer-booktabs")
require(stargazer)

##### Felm corrected when fes are nested in clusters #######
felmNest <- function(form, data, ...)
{
  m <- felm(form, data, exactDOF = TRUE, keepCX = TRUE)
  
  mlm <- lm(m$cY ~ m$cX)
  
  clse <- cluster.vcov(mlm, m$clustervar, df_correction = FALSE) 
  
  ses <- coeftest(mlm , clse)
  m$cse <- ses[2:nrow(ses),2]
  m$ctval <- ses[2:nrow(ses),3]
  m$cpval <- ses[2:nrow(ses),4]
  m$cY <- NULL
  m$cX <- NULL
  return(m)
}

####### Rounding for table creation #######
meanPr <- function(data, digits = 3)
{
  m <- mean(data, na.rm = TRUE)
  return(formatC(m, digits = digits, format = "f"))
}

sdPr <- function(data, digits = 3)
{
  m <- sd(data, na.rm = TRUE)
  return(formatC(m, digits = digits, format = "f"))
}

tdif <- function(var1, var2, digits = 2)
{
  m <- t.test(var1, var2)
  return(formatC(m$statistic, digits = digits, format = "f"))
}

difPr <- function(var1, var2, digits = 2)
{
  m <- mean(var1, na.rm = TRUE) - mean(var2, na.rm = TRUE)
  return(formatC(m, digits = digits, format = "f"))
}

##### Get the mode of a data set #######
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

###### Better config for DataExplorer #######
# create_report(data, y = "", config = dconfig)
dconfig <- list(
  "introduce" = list(),
  "plot_str" = list(
    "type" = "diagonal",
    "fontSize" = 35,
    "width" = 1000,
    "margin" = list("left" = 350, "right" = 250)
  ),
  "plot_missing" = list(),
  "plot_histogram" = list(),
  "plot_bar" = list(),
  "plot_correlation" = list("use" = "pairwise.complete.obs"),
  "plot_boxplot" = list(),
  "plot_scatterplot" = list()
)

####### Tic/toc functionality ########
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic = 0
  tic <- proc.time()[type]
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

toc <- function()
{
  type <- get(".type", envir=baseenv())
  toc = 0
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  print(toc - tic)
  invisible(toc)
}

####### Add the following code to TexStudio options/build  ########
# user commands to allow for .Rnw files
# "Rscript.exe" -e  "knitr::knit2pdf('%.Rnw')" | pdflatex -synctex=1 -interaction=nonstopmode %.tex | txs:///view-pdf


####### Add scalebox and fix notes in stargazer ########
fixStargazer <- function(tab, nmodels, scalesize = .8)
{
  # tab should be a character vector of captured output from stargazer
  # this version requires you to use booktabs in your preamble of your tex doc
  # adding a scalebox
  scaleChar <- paste("\\scalebox{", scalesize, "}{", sep = "")
  a <- sub("\\begin{tabular}", paste(scaleChar, "\\begin{tabular}", sep = ""), tab, fixed = TRUE)
  b <- sub("\\end{tabular}", "\\end{tabular}}", a, fixed = TRUE)
  
  # making the notes go across the table
  noteChar <- paste("\\multicolumn{", nmodels + 1, "}{l}{\\parbox", sep = "")
  patChar <- paste(" & \\multicolumn{", nmodels, "}{l}{\\parbox", sep = "")
  c <- sub(patChar,  noteChar, b, fixed = TRUE)

  
  # dropping the comment at the top
  c[2:4] <- ""
  # printing output
  cat(c)
}

  
#### Fix the weird unlabeled dataset problem ####
# This is also fixed by loading the packages in the order from the master file
# The problem occurs when you load HMisc and tidyverse in the wrong order
labelDataset <- function(data) {
  correctLabel <- function(x) {
    
    if(!is.null(attributes(x)$labels)) {
      class(attributes(x)$labels) <- typeof(x)
    }
    return(x)
  }
  for(i in colnames(data)) {
    data[, i] <- correctLabel(data[, i])
  }
  return(data)
}


# twoSampleSumStats <- function(data, treatvar, varnames)
# {
#   tab <- c()
#   tab[1] <- "

















