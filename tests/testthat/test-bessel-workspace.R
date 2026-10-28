test_that("Bessel working arrays do not allocate an R vector per observation", {
  skip_if_not(capabilities("profmem"))
  d <- numerical_fixture("ST",N=120L)
  evaluate <- native("pmbfa_evaluate")
  allocated <- function(nu) {
    params <- d$params
    for(g in seq_along(params$args))params$args[[g]]$nu <- nu
    invisible(evaluate(d$x,params))
    gc()
    path <- tempfile()
    on.exit(unlink(path),add=TRUE)
    utils::Rprofmem(path)
    on.exit(utils::Rprofmem(NULL),add=TRUE)
    result <- evaluate(d$x,params)
    utils::Rprofmem(NULL)
    lines <- readLines(path,warn=FALSE)
    bytes <- sum(as.numeric(sub(" .*","",lines[grepl("^[0-9]+ :",lines)])))
    expect_true(all(is.finite(result$loglik)))
    bytes
  }
  # Both orders use the finite Rmath kernel, with quite different work lengths.
  # The output shapes are identical; only the reusable C++ buffer should grow.
  small <- allocated(6)
  large <- allocated(200)
  expect_lte(abs(large-small),8192)
})
