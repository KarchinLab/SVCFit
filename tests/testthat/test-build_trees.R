test_that("build_trees errors when run_tree=TRUE but run_clustering=FALSE", {
  expect_error(
    build_trees(run_tree = TRUE, run_clustering = FALSE),
    "run_tree = TRUE requires run_clustering = TRUE"
  )
})

test_that("build_trees errors when run_clustering=TRUE but pair_path is NULL", {
  expect_error(
    build_trees(run_clustering = TRUE, pair_path = NULL),
    "pair_path must be provided"
  )
})

test_that("build_trees errors when run_clustering=TRUE but pur_path is NULL", {
  expect_error(
    build_trees(run_clustering = TRUE, pair_path = "p.txt", pur_path = NULL),
    "pur_path must be provided"
  )
})

test_that("build_trees errors when run_clustering=TRUE but data_dir is NULL", {
  expect_error(
    build_trees(run_clustering = TRUE, pair_path = "p.txt",
                pur_path = "pur.txt", data_dir = NULL),
    "data_dir must be provided"
  )
})

test_that("build_trees returns two-element list with NULLs when no stages run", {
  result <- build_trees(run_clustering = FALSE, run_tree = FALSE)

  expect_type(result, "list")
  expect_length(result, 2)
  expect_null(result[[1]])
  expect_null(result[[2]])
})
