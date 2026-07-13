library(testthat)
library(GeneSCARAB)

test_that("create_gene_list_go returns a list", {
    library(org.Otauri.eg.db)
    go.list.test <- create_gene_list_go(
        go_vector = "GO:0005515",
        org.package = "org.Otauri.eg.db", go_column = "GO",
        id_column = "GID"
    )
    expect_type(go.list.test, "list")
})

test_that("throws an error if GO column do not coincide with
          that of the package", {
    expect_error(
        go.list.test <- create_gene_list_go(
            go_vector = "GO:0005515",
            org.package = "org.Otauri.eg.db", go_column = "GOS",
            id_column = "GID"
        ),
        "ID or GO column names do not match"
    )
})

test_that("throws an error if ID column do not coincide with that
          of the package", {
    expect_error(
        go.list.test <- create_gene_list_go(
            go_vector = "GO:0005515",
            org.package = "org.Otauri.eg.db", go_column = "GO",
            id_column = "GIDS"
        ),
        "ID or GO column names do not match"
    )
})

test_that("throws an error if the selected annotation package
          is not available", {
    expect_error(
        go.list.test <- create_gene_list_go(
            go_vector = "GO:0005515",
            org.package = "org.Ag.eg.db", go_column = "GO",
            id_column = "GID"
        ),
        "Annotation package name do not correspond"
    )
})
