# Script for generation of inst/extdata/org_Otaurireduced_eg_db.sqlite.gz
# dataset

# This dataset is an SQlite containing GO and KO annnotation info for 
# Ostreococcus tauri in AnnotationForge's output format

# Original annotation package org.Otauri.eg.db was download from
# https://github.com/fran-romero-campero/AlgaeFUN/tree/master/packages/annotation_packages
# and installed as source

# To reduce its size, first only GO info was retained
# Extract the original information and create table
library(org.Otauri.eg.db)
library(GO.db)
library(seqinr)
library(dplyr)
library(tidyr)
columns(org.Otauri.eg.db)


gos <- AnnotationDbi::select(org.Otauri.eg.db,
                             columns = c("GO"),
                             keys=keys(org.Otauri.eg.db,keytype = "GID"))

df.go <- data.frame(GID = gos$GID,GO = gos$GO, EVIDENCE = "ISS")
df.go <- df.go[df.go$GO != "",]
df.go <- df.go[!duplicated(df.go),]
df.go <- df.go[df.go$GID != "",]
df.go <- df.go[!is.na(df.go$GID),]

# After that, KO information for the corresponding genes was generated via KAAS
# (https://www.genome.jp/kaas-bin/kaas_main) on the complete genome available at
# https://github.com/fran-romero-campero/AlgaeFUN/tree/master/genomes and
# saved as a tsv named ko_table_base.tsv

# Next we extract these relationships

ko1 <- read.csv("ko_table_base.tsv", sep="\t", header=T)
head(ko1)

ko2 <- ko1[grep("ostta", ko1$gene),c(2,3)]

df.ko <- data.frame(GID = ko2$gene, KO = ko2$ko)
df.ko <- df.ko[!duplicated(df.ko),]
df.ko <- df.ko[df.ko$KO != "",]
df.ko <- df.ko[df.ko$GID != "",]
df.ko <- df.ko[!is.null(df.ko$KO),]
dim(df.ko)

# Once we have the GO and KO info, we use AnnotationForge to build an annotation
# package
AnnotationForge::makeOrgPackage(GO = df.go,
                                KO = df.ko,
                                version="1.0",
                                maintainer="Marcos Ramos-Gonzalez <mramos5@us.es>",
                                author="Marcos Ramos-Gonzalez",
                                outputDir = ".",
                                tax_id="70448",
                                genus="Ostreococcus",
                                species="taurireduced",
                                goTable="GO",
                                verbose = TRUE)

# And we install it
install.packages("./org.Otaurireduced.eg.db/", repos=NULL, type = "source")

# Once the annotation package is created, we clean its SQLite database 
library(AnnotationDbi)
library(org.Otaurireduced.eg.db)

saveDb(
  org.Otaurireduced.eg.db::org.Otaurireduced.eg.db,
  file = "../extdata/org_Otaurireduced_eg_db.sqlite"
)

library(DBI)
library(RSQLite)

con <- dbConnect(
  SQLite(),
  dbname = "../extdata/org_Otaurireduced_eg_db.sqlite"
)

dbExecute(con, "VACUUM")

dbDisconnect(con)

# Once cleaned, we compress the database for storing
gzfile <- gzfile("../extdata/org_Otaurireduced_eg_db.sqlite.gz", "wb")
infile <- file("../extdata/org_Otaurireduced_eg_db.sqlite", "rb")

while (length(x <- readBin(infile, "raw", 1024^2)) > 0) {
  writeBin(x, gzfile)
}

close(infile)
close(gzfile)

