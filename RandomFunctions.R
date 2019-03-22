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
# Mean
meanPr <- function(data, digits = 3)
{
  m <- mean(data, na.rm = TRUE)
  return(formatC(m, digits = digits, format = "f"))
}

# Standard Deviation
sdPr <- function(data, digits = 3)
{
  m <- sd(data, na.rm = TRUE)
  return(formatC(m, digits = digits, format = "f"))
}

# T-test for difference in means
tdif <- function(var1, var2, digits = 3)
{
  # Come back and add stars
  m <- t.test(var1, var2)
  return(formatC(m$statistic, digits = digits, format = "f"))
}

# Difference in means
difPr <- function(var1, var2, digits = 3)
{
  m <- mean(var1, na.rm = TRUE) - mean(var2, na.rm = TRUE)
  return(formatC(m, digits = digits, format = "f"))
}

# Within R-squared for FELM
WR2 <- function(model, digits = 3){
return(round(summary(model)$P.r.squared, digits = digits))
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

####### Add the following code to TexStudio options/build for .rnw file compiling  ########
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
  d <- sub(patChar,  noteChar, b, fixed = TRUE)

  
  # dropping the comment at the top
  d[2:4] <- ""
  
  # getting line endings
  a <- gsub("\\\\", "\\\\ \n", d, fixed = TRUE)
  # printing output (if trunctated, go to RStudiio Options, Code, Display, and set the number at the bottom higher.
  cat(a)
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


#### Produce a well-labeled table for sum stats of a treated and non-treated sample ####
twoSampleSumStats <- function(data, 
                              treatvar, 
                              varnames, 
                              fancyvarnames, 
                              caption = "", 
                              label = "", 
                              note = "", 
                              treatcolname = "Treated",
                              nontreatcolname = "Not Treated", 
                              scale = 1,
                              parbox = "21cm")
{
  # need to add a way to order the variables like varname orders them.
  # need to add a version that works for beamer (dropping cmidrule).
  # need to add a way to save it to a file instead of printing.
  # should add an option for digits
  
  # pull variables from data
  l <- which(names(data) %in% varnames)
  data <- data[,l]
  
  # separete treat vs nontreat samples
  treatdata    <- data[treatvar == 1, ]
  nontreatdata <- data[treatvar == 0, ]
  vars <- names(data)
  
  tab <- vector()
  # top of the table
  tab[1] <- paste0("\\begin{table}[h!] \\centering \n")
  tab[2] <- paste0("\\caption{",caption,"} \n")
  tab[3] <- paste0("\\label{",label,"} \n")
  tab[4] <- paste0("\\scalebox{", scale, "}{\\begin{tabular}{@{\\extracolsep{2pt}}lD{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3}} \n")
  tab[5] <- paste0("\\hline \n")
  tab[6] <- paste0("\\\\[0.5mm] \n")
  tab[7] <- paste0("Variable & \\multicolumn{2}{c}{Full Sample:} & \\multicolumn{2}{c}{",treatcolname,":} & \\multicolumn{2}{c}{",nontreatcolname,":}  & \\multicolumn{2}{c}{Diff. Means:}\\\\ \n")
  tab[8] <- paste0("\\cmidrule{2-3} \\cmidrule{4-5} \\cmidrule{6-7} \\cmidrule{8-9} \n")
  tab[9] <- paste0("\\\\[-1.8ex] \n")
  tab[10]<- paste0("& \\multicolumn{1}{c}{Mean} & \\multicolumn{1}{c}{St. Dev.} & \\multicolumn{1}{c}{Mean} & \\multicolumn{1}{c}{St. Dev.} & \\multicolumn{1}{c}{Mean} & \\multicolumn{1}{c}{St. Dev.} & \\multicolumn{1}{c}{Diff.} &\\multicolumn{1}{c}{t}\\\\ \n")
  tab[11]<- paste0("\\hline \\\\[-1.8ex] \n")
  
  row = length(tab)
  # gen sum stats
  for(i in 1:length(vars)){
    fm  <- meanPr(data[[i]])
    fs  <- sdPr(data[[i]])
    tm  <- meanPr(treatdata[[i]])
    ts  <- sdPr(treatdata[[i]])
    ntm <- meanPr(nontreatdata[[i]])
    nts <- sdPr(nontreatdata[[i]])
    d   <- difPr(treatdata[[i]], nontreatdata[[i]])
    dm  <- tdif(treatdata[[i]], nontreatdata[[i]])
    
    row <- row + 1
    tab[row] <- paste(fancyvarnames[i],"&",fm,"&",fs,"&",tm,"&",ts,"&",ntm,"&",nts,"&",d,"&", dm,"\\\\ \n")
  }
  
  # bottom
  tab[row+1] <- paste0("N & \\multicolumn{2}{c}{",nrow(data),"} & \\multicolumn{2}{c}{",nrow(treatdata),"} & \\multicolumn{2}{c}{",nrow(nontreatdata),"} & \\multicolumn{2}{c}{-} \\\\ \n")
  tab[row+2] <- paste0("\\hline \\\\[-1.8ex] \n")
  if (note != "") {
    tab[row+3] <- paste0("\\multicolumn{9}{l}{\\parbox[t]{",parbox,"}{\\textbf{Note:}",note,"}} \\\\ \n")
  }  else {
    tab[row+3] <- paste0(" \n")
  }
  tab[row+4] <- paste0("\\end{tabular}} \n")
  tab[row+5] <- paste0("\\end{table} \n")
  
  # print table
  cat(tab)
}


#### Function to make FELM formula from strings ####
felmForm <- function(y, treat, x, fe = 0, inst = 0, clus = 0)
{
  form <- paste(y, "~", paste(treat, collapse = "+"), "+", paste(x, collapse = " + "), "|", paste(fe, collapse = " + "), "|", inst, "|", paste(clus, collapse = " + "))
  form <- formula(form, env = <environment: R_GlobalEnv>)
  return(form)
}


















