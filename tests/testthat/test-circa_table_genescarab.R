library(testthat)
library(GeneSCARAB)

test_that("example dataset has expected structure", {
    data("circa_table_genescarab")
    expect_s3_class(circa_table_genescarab, "data.frame")
    expect_gt(nrow(circa_table_genescarab), 0)
    expect_equal(ncol(circa_table_genescarab), 15)
    expect_true("ld.peak.time.hours" %in% colnames(circa_table_genescarab))
    expect_true("sd.peak.time.hours" %in% colnames(circa_table_genescarab))
    expect_type(circa_table_genescarab$sd.peak.time.hours, "double")
    expect_type(circa_table_genescarab$ld.peak.time.hours, "double")
})