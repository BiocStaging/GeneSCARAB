library(testthat)
library(GeneSCARAB)

test_that("circular_boxplot runs without errors", {
    library(org.Otauri.eg.db)
    data("circa_table_genescarab")
    total_phases_table_sd <- data.frame(
        names = rownames(circa_table_genescarab),
        phase = as.numeric(
            circa_table_genescarab[["sd.peak.time.hours"]]
        )
    )

    functional_data <- select(org.Otauri.eg.db,
        keys = keys(org.Otauri.eg.db, keytype = "GID"),
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
        org.package = "org.Otauri.eg.db",
        go_column = "GO", id_column = "GID"
    )
    phases.list.sd <- gene_list_to_phases(
        go.list.test,
        total_phases_table_sd
    )
    phases.list.sd.clean <- phases.list.sd[which(
        sapply(phases.list.sd, nrow) != 0
    )][seq_len(3)]
    go.circa.table.sd <- complete_circular_table(
        phase.list = phases.list.sd.clean,
        hr.on.large.sets.th = 200,
        force.hr.th = 0.04,
        rao.on.large.sets.th = 200,
        force.rao.th = 0.04,
        iter.rao = 999, iter.hr = 999
    )

    expect_no_error(
        circular_boxplot(go.circa.table.sd[seq_len(2), ],
            color.palette = "Tam"
        )
    )
})
