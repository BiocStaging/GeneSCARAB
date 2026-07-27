library(testthat)
library(GeneSCARAB)

test_that("test_two_dist returns a data.frame with two columns", {
    library(org.Otauri.eg.db)
    library(GeneSCARAB)

    data("circa_table_genescarab")
    total_phases_table_sd <- data.frame(
        names = rownames(circa_table_genescarab),
        phase = as.numeric(
            circa_table_genescarab[["sd.peak.time.hours"]]
        )
    )
    total_phases_table_ld <- data.frame(
        names = rownames(circa_table_genescarab),
        phase = as.numeric(
            circa_table_genescarab[["ld.peak.time.hours"]]
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
        go.list.test, total_phases_table_sd
    )
    phases.list.sd.clean <- phases.list.sd[which(
        sapply(phases.list.sd, nrow) != 0
    )][seq_len(3)]
    phases.list.ld <- gene_list_to_phases(
        go.list.test, total_phases_table_sd
    )
    phases.list.ld.clean <- phases.list.sd[which(
        sapply(phases.list.sd, nrow) != 0
    )][names(phases.list.sd.clean)]
    diff_distributed_gos <- test_two_dist(
        phases.list.ld.clean, phases.list.sd.clean
    )

    expect_type(diff_distributed_gos, "list")
})

test_that("different lengths of input lists produce an error", {
    expect_error(
        test_two_dist(
            list(
                A = data.frame(names = paste0(
                    "gene", c("A", "B", "C")
                ), phase = rep(0, 3)),
                B = data.frame(names = paste0(
                    "gene", c("A", "B", "C")
                ), phase = rep(0, 3))
            ),
            list(A = data.frame(names = paste0(
                "gene", c("A", "B", "C")
            ), phase = rep(0, 3)))
        ),
        "Length of lists differ"
    )
})

test_that("column names of data frames are names and phase", {
    expect_error(
        test_two_dist(
            list(
                A = data.frame(names = paste0(
                    "gene", c("A", "B", "C")
                ), phase = rep(0, 3)),
                B = data.frame(ids = paste0(
                    "gene", c("A", "B", "C")
                ), phase = rep(0, 3))
            ),
            list(A = data.frame(names = paste0(
                "gene", c("A", "B", "C")
            ), phase = rep(0, 3)))
        ),
        "Length of lists differ"
    )
})

test_that("different names or order of input lists produce a warning", {
    expect_warning(
        test_two_dist(
            list(
                A = data.frame(names = paste0(
                    "gene", c("A", "B", "C")
                ), phase = rep(0, 3)),
                B = data.frame(names = paste0(
                    "gene", c("A", "B", "C")
                ), phase = rep(0, 3))
            ),
            list(
                C = data.frame(names = paste0(
                    "gene", c("A", "B", "C")
                ), phase = rep(0, 3)),
                D = data.frame(names = paste0(
                    "gene", c("A", "B", "C")
                ), phase = rep(0, 3))
            )
        ),
        "Set names do not match"
    )
})
