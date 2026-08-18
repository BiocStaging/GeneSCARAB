library(testthat)
library(GeneSCARAB)

test_that("example dataset has expected structure", {
  org_Otaurireduced_eg_db <- load_example_annot()
  functional_data <- AnnotationDbi::select(org_Otaurireduced_eg_db,
                                           keys = AnnotationDbi::keys(
                                             org_Otaurireduced_eg_db, keytype = "GID"),
                                           columns = c("GID", "GO", "KO"))
  expect_s4_class(org_Otaurireduced_eg_db, "OrgDb")
  expect_gt(nrow(functional_data), 0)
  expect_true("GID" %in% AnnotationDbi::columns(org_Otaurireduced_eg_db))
  expect_true("GO" %in% AnnotationDbi::columns(org_Otaurireduced_eg_db))
  expect_true("KO" %in% AnnotationDbi::columns(org_Otaurireduced_eg_db))
  expect_type(functional_data$GID, "character")
  expect_type(functional_data$GO, "character")
  expect_type(functional_data$KO, "character")
})
