# Simple 2-clone CCF matrix used across tests:
#   clone 1: high CCF (0.8, 0.6) — candidate parent
#   clone 2: low CCF  (0.4, 0.2) — candidate child
mcf2 <- matrix(c(0.8, 0.4, 0.6, 0.2), nrow = 2, ncol = 2,
               dimnames = list(c("1", "2"), c("t1", "t2")))

test_that("prepareGraph includes a root edge for every clone", {
  g <- prepareGraph(mcf2, thresh = 0.1)

  expect_s3_class(g, "data.frame")
  expect_true(all(c("edge", "parent", "child") %in% names(g)))
  expect_setequal(g$child[g$parent == "root"], c("1", "2"))
})

test_that("prepareGraph edge format is 'parent->child'", {
  g <- prepareGraph(mcf2, thresh = 0.1)

  expect_true(all(grepl("->", g$edge)))
})

test_that("filterEdgesBasedOnCCFs removes edges where child CCF exceeds parent", {
  g  <- prepareGraph(mcf2, thresh = 0.1)
  fg <- filterEdgesBasedOnCCFs(g, mcf2, thresh = 0.1)

  # clone 2 (CCF 0.4, 0.2) cannot be parent of clone 1 (CCF 0.8, 0.6)
  expect_false(any(fg$parent == "2" & fg$child == "1"))
  # clone 1 can be parent of clone 2
  expect_true(any(fg$parent == "1" & fg$child == "2"))
})

test_that("verticesInGraph returns all unique node labels", {
  edges <- data.frame(
    edge   = c("root->1", "1->2"),
    parent = c("root", "1"),
    child  = c("1", "2"),
    stringsAsFactors = FALSE
  )
  v <- verticesInGraph(edges)
  expect_setequal(v, c("root", "1", "2"))
})

test_that("satisfiesSumCondition passes when children CCF <= parent CCF", {
  edges <- dplyr::tibble(parent = c("root", "1"), child = c("1", "2"),
                         edge   = c("root->1", "1->2"))
  expect_true(satisfiesSumCondition(edges, mcf2, thresh = 0.2))
})

test_that("satisfiesSumCondition fails when children CCF sum exceeds parent", {
  # Two children of root, together > 1
  mcf3 <- matrix(c(0.7, 0.7, 0.5, 0.5), nrow = 2, ncol = 2,
                 dimnames = list(c("1", "2"), c("t1", "t2")))
  edges <- dplyr::tibble(parent = c("root", "root"),
                         child  = c("1", "2"),
                         edge   = c("root->1", "root->2"))
  expect_false(satisfiesSumCondition(edges, mcf3, thresh = 0.2))
})
