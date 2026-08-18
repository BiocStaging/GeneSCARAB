library(testthat)
library(GeneSCARAB)

test_that("test_against_gen_dist returns a data.frame", {
    org_Otaurireduced_eg_db <- load_example_annot()
    data("circa_table_genescarab")
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
    not_experimentally_distributed_gos <- test_against_gen_dist(
        phases.list.sd.clean,
        total_phases_table_sd
    )

    expect_type(not_experimentally_distributed_gos, "list")
})

test_that("throws an error if phase.list is not a list of data
          frames with columns named names and phase", {
    org_Otaurireduced_eg_db <- load_example_annot()
    data("circa_table_genescarab")
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
    colnames(phases.list.sd.clean[[1]]) <- c("names", "phases")

    expect_error(
        not_experimentally_distributed_gos <- test_against_gen_dist(
            phases.list.sd.clean,
            total_phases_table_sd
        ),
        "Check that all tables in the list use names and phase"
    )
})

test_that("throws an error if total.phase.table does not contain
          a column named phase", {
    org_Otaurireduced_eg_db <- load_example_annot()
    data("circa_table_genescarab")
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
    colnames(total_phases_table_sd) <- c("names", "phases")

    expect_error(
        not_experimentally_distributed_gos <- test_against_gen_dist(
            phases.list.sd.clean,
            total_phases_table_sd
        ),
        "must have a column named phase"
    )
})
