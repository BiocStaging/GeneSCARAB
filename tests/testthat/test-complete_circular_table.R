library(testthat)
library(GeneSCARAB)

test_that("complete_circular_table returns a data.frame", {
    org_Otaurireduced_eg_db <- load_example_annot()
    data("circa_table_genescarab")
    total_phases_table_sd <- data.frame(
        names = rownames(circa_table_genescarab),
        phase = as.numeric(circa_table_genescarab[["sd.peak.time.hours"]])
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
    go.circa.table.sd <- complete_circular_table(
        phase.list = phases.list.sd.clean,
        hr.on.large.sets.th = 200,
        force.hr.th = 0.04,
        rao.on.large.sets.th = 200,
        force.rao.th = 0.04,
        iter.rao = 999, iter.hr = 999
    )
    expect_type(go.circa.table.sd, "list")
})

test_that("throws an error if phase.list is not a list of data
          frames with columns named names and phase", {
    org_Otaurireduced_eg_db <- load_example_annot()
    data("circa_table_genescarab")
    total_phases_table_sd <- data.frame(
        names = rownames(circa_table_genescarab),
        phase = as.numeric(circa_table_genescarab[["sd.peak.time.hours"]])
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
        go.circa.table.sd <- complete_circular_table(
            phase.list = phases.list.sd.clean,
            hr.on.large.sets.th = 200,
            force.hr.th = 0.04,
            rao.on.large.sets.th = 200,
            force.rao.th = 0.04,
            iter.rao = 999, iter.hr = 999
        ),
        "Check that all tables in the list use names and phase"
    )
})

test_that("throws an error if the column phase of any of the
          tables is not numeric", {
    org_Otaurireduced_eg_db <- load_example_annot()
    data("circa_table_genescarab")
    total_phases_table_sd <- data.frame(
        names = rownames(circa_table_genescarab),
        phase = as.numeric(circa_table_genescarab[["sd.peak.time.hours"]])
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

    phases.list.sd.clean[[1]]$phase <- "a"

    expect_error(
        go.circa.table.sd <- complete_circular_table(
            phase.list = phases.list.sd.clean,
            hr.on.large.sets.th = 200,
            force.hr.th = 0.04,
            rao.on.large.sets.th = 200,
            force.rao.th = 0.04,
            iter.rao = 999, iter.hr = 999
        ),
        "must be a numeric vector of phases"
    )
})

test_that("throws an error if phases_list is not a named
          list", {
    org_Otaurireduced_eg_db <- load_example_annot()
    data("circa_table_genescarab")
    total_phases_table_sd <- data.frame(
        names = rownames(circa_table_genescarab),
        phase = as.numeric(circa_table_genescarab[["sd.peak.time.hours"]])
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

    names(phases.list.sd.clean) <- NULL

    expect_error(
        go.circa.table.sd <- complete_circular_table(
            phase.list = phases.list.sd.clean,
            hr.on.large.sets.th = 200,
            force.hr.th = 0.04,
            rao.on.large.sets.th = 200,
            force.rao.th = 0.04,
            iter.rao = 999, iter.hr = 999
        ),
        "must be a named list"
    )
})
