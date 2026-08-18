library(testthat)
library(GeneSCARAB)

test_that("create_gene_list_go returns a list", {
    org_Otaurireduced_eg_db <- load_example_annot()
    go.list.test <- create_gene_list_go(
        go_vector = "GO:0005515",
        org.package = org_Otaurireduced_eg_db, go_column = "GO",
        id_column = "GID"
    )
    expect_type(go.list.test, "list")
})

test_that("throws an error if GO column do not coincide with
          that of the package", {
    org_Otaurireduced_eg_db <- load_example_annot()
    expect_error(
        go.list.test <- create_gene_list_go(
            go_vector = "GO:0005515",
            org.package = org_Otaurireduced_eg_db, go_column = "GOS",
            id_column = "GID"
        ),
        "ID or GO column names do not match"
    )
})

test_that("throws an error if ID column do not coincide with that
          of the package", {
    org_Otaurireduced_eg_db <- load_example_annot()
    expect_error(
        go.list.test <- create_gene_list_go(
            go_vector = "GO:0005515",
            org.package = org_Otaurireduced_eg_db, go_column = "GO",
            id_column = "GIDS"
        ),
        "ID or GO column names do not match"
    )
})


