library(testthat)
library(GeneSCARAB)

test_that("multimodal_analysis returns a list", {
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
    subset_gos <- complete_gos[sample(seq_len(
        length(complete_gos)
    ), 50, replace = FALSE)]
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
    )]
    go.circa.table.sd <- complete_circular_table(
        phase.list = phases.list.sd.clean,
        hr.on.large.sets.th = 200,
        force.hr.th = 0.04,
        rao.on.large.sets.th = 200,
        force.rao.th = 0.04,
        iter.rao = 999, iter.hr = 999
    )
    multi_table <- go.circa.table.sd[which(
        go.circa.table.sd$hr_p_value < 0.05 &
            go.circa.table.sd$n > 2
    ), ]
    multi_go <- rownames(multi_table)[1]
    phase_table_multi_go <- phases.list.sd.clean[[multi_go]]
    multimodal_list <- multimodal_analysis(phase_table_multi_go)
    expect_type(multimodal_list, "list")
})

test_that("throws an error if phase_table does not contain a column
          named phase", {
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
    subset_gos <- complete_gos[sample(seq_len(
        length(complete_gos)
    ), 50, replace = FALSE)]
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
    )]
    go.circa.table.sd <- complete_circular_table(
        phase.list = phases.list.sd.clean,
        hr.on.large.sets.th = 200,
        force.hr.th = 0.04,
        rao.on.large.sets.th = 200,
        force.rao.th = 0.04,
        iter.rao = 999, iter.hr = 999
    )
    multi_table <- go.circa.table.sd[which(
        go.circa.table.sd$hr_p_value < 0.05 &
            go.circa.table.sd$n > 2
    ), ]
    multi_go <- rownames(multi_table)[1]
    phase_table_multi_go <- phases.list.sd.clean[[multi_go]]
    colnames(phase_table_multi_go) <- c("names", "phases")

    expect_error(
        multimodal_list <- multimodal_analysis(phase_table_multi_go),
        "No column named phase detected in phase_table"
    )
})
