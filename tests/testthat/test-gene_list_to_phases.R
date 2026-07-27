library(testthat)
library(GeneSCARAB)

test_that("gene_list_to_phases returns a list", {
    library(org.Otauri.eg.db)
    data("circa_table_genescarab")

    total_phases_table_sd <- data.frame(
        names = rownames(circa_table_genescarab),
        phase = as.numeric(
            circa_table_genescarab[["sd.peak.time.hours"]]
        )
    )
    go.list.test <- create_gene_list_go(
        go_vector = "GO:0005515",
        org.package = "org.Otauri.eg.db", go_column = "GO",
        id_column = "GID"
    )
    list_result <- gene_list_to_phases(go.list.test, total_phases_table_sd)
    expect_type(list_result, "list")
})

test_that("throws an error if phase table is not a two-column
          data.frame with columns named names and phase", {
    library(org.Otauri.eg.db)
    data("circa_table_genescarab")

    total_phases_table_sd <- data.frame(
        names = rownames(circa_table_genescarab),
        phase = as.numeric(
            circa_table_genescarab[["sd.peak.time.hours"]]
        )
    )

    go.list.test <- create_gene_list_go(
        go_vector = "GO:0005515",
        org.package = "org.Otauri.eg.db", go_column = "GO",
        id_column = "GID"
    )

    colnames(total_phases_table_sd) <- c("ids", "acrophase")
    expect_error(
        list_result <- gene_list_to_phases(
            go.list.test,
            total_phases_table_sd
        ),
        "names and phase"
    )
})
