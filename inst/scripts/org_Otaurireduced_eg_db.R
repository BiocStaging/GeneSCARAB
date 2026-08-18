#' org_Otaurireduced_eg_db
#'
#' SQlite containing GO and KO annnotation info for Ostreococcus tauri in 
#' AnnotationForge's output format. Original genome and annotation package 
#' were download from https://github.com/fran-romero-campero/AlgaeFUN and KO 
#' information was generated via KAAS 
#' (https://www.genome.jp/kaas-bin/kaas_main).
#'
#' @format OrgDb object
"org_Otaurireduced_eg_db"


# Once the annotation package is created, we will save its database 
library(AnnotationDbi)
library(org.Otaurireduced.eg.db)

saveDb(
  org.Otaurireduced.eg.db::org.Otaurireduced.eg.db,
  file = "../extdata/org_Otaurireduced_eg_db.sqlite"
)

# Clean SQLite database
library(DBI)
library(RSQLite)

con <- dbConnect(
  SQLite(),
  dbname = "../extdata/org_Otaurireduced_eg_db.sqlite"
)

dbExecute(con, "VACUUM")

dbDisconnect(con)

# Compress database for storing
gzfile <- gzfile("../extdata/org_Otaurireduced_eg_db.sqlite.gz", "wb")
infile <- file("../extdata/org_Otaurireduced_eg_db.sqlite", "rb")

while (length(x <- readBin(infile, "raw", 1024^2)) > 0) {
  writeBin(x, gzfile)
}

close(infile)
close(gzfile)

