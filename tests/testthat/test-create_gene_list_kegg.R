library(testthat)
library(GeneSCARAB)

test_that("create_gene_list_kegg returns a list", {
    org_Otaurireduced_eg_db <- load_example_annot()
    ko.list.test <- create_gene_list_kegg(
        ko_vector = "K10666", org.package = org_Otaurireduced_eg_db,
        ko_column = "KO", id_column = "GID", ko_prefix = "map",
        species = "plants"
    )
    expect_type(ko.list.test, "list")
})

test_that("throws an error if KO column do not coincide with
          that of the package", {
    org_Otaurireduced_eg_db <- load_example_annot()
    expect_error(
        ko.list.test <- create_gene_list_kegg(
            ko_vector = "K10666", org.package = org_Otaurireduced_eg_db,
            ko_column = "KOT", id_column = "GID", ko_prefix = "map",
            species = "plants"
        ),
        "ID or KO column names do not match"
    )
})

test_that("throws an error if ID column do not coincide with
          that of the package", {
    org_Otaurireduced_eg_db <- load_example_annot()
    expect_error(
        ko.list.test <- create_gene_list_kegg(
            ko_vector = "K10666", org.package = org_Otaurireduced_eg_db,
            ko_column = "KO", id_column = "GIDS", ko_prefix = "map",
            species = "plants"
        ),
        "ID or KO column names do not match"
    )
})

test_that("throws an error if prefix is not one of the supported ones", {
    org_Otaurireduced_eg_db <- load_example_annot()
    expect_error(
        ko.list.test <- create_gene_list_kegg(
            ko_vector = "K10666", org.package = org_Otaurireduced_eg_db,
            ko_column = "KO", id_column = "GID", ko_prefix = "ole",
            species = "plants"
        ),
        "prefix must be either map or ko"
    )
})

test_that("throws an error if species is not one of the supported ones", {
    org_Otaurireduced_eg_db <- load_example_annot()
    expect_error(
        ko.list.test <- create_gene_list_kegg(
            ko_vector = "K10666", org.package = org_Otaurireduced_eg_db,
            ko_column = "KO", id_column = "GID", ko_prefix = "map",
            species = "horses"
        ),
        "Not supported species group"
    )
})
