library(testthat)
library(GeneSCARAB)

test_that("circular_dotplot runs without errors", {
    data("circa_table_genescarab")
    org_Otaurireduced_eg_db <- load_example_annot()
    total_phases_table_sd <- data.frame(
        names = rownames(circa_table_genescarab),
        phase = as.numeric(
            circa_table_genescarab[["sd.peak.time.hours"]]
        )
    )

    functional_data <- AnnotationDbi::select(org_Otaurireduced_eg_db,
        keys = AnnotationDbi::keys(org_Otaurireduced_eg_db, keytype = "GID"),
        columns = c("GID", "GO")
    )
    complete_gos <- unique(functional_data$GO)
    complete_gos <- complete_gos[!is.na(complete_gos)]

    set.seed(2345)
    subset_gos <- complete_gos[sample(
        seq_len(length(complete_gos)), 50,
        replace = FALSE
    )]
    go.list.test <- create_gene_list_go(
        go_vector = subset_gos,
        org.package = org_Otaurireduced_eg_db,
        go_column = "GO", id_column = "GID"
    )
    phases.list.sd <- gene_list_to_phases(
        go.list.test, total_phases_table_sd
    )
    phases.list.sd.clean <- phases.list.sd[which(
        sapply(phases.list.sd, nrow) != 0
    )][seq_len(3)]

    expect_no_error(
        circular_dotplot(phases.list.sd.clean, color.palette = "Tam")
    )
})
