verbosity_fit <- function(entry, control, verbose=FALSE, criterion="BIC") {
  skewed <- grepl("pmsbfa",entry,fixed=TRUE)
  d <- numerical_fixture(if(skewed) "ST" else "Normal",N=60L)
  arguments <- list(x=d$x,G=2L,inits=d$params,control=control,
    modelSubset=d$params$moddat$name,verbose=verbose,criterion=criterion)
  if(startsWith(entry,"search_")) {
    arguments$max_q <- 1L
    arguments$max_r <- 1L
  } else {
    arguments$q <- 1L
    arguments$r <- 1L
  }
  do.call(getExportedValue("pmbfa",entry),arguments)
}

test_that("all four entry points separate model summaries from likelihood tracking", {
  controls <- quick_control(max.iter=3L,min.iter=1L,tol=1e6,min.component=0)
  for(entry in c("pmbfa","pmsbfa","search_pmbfa","search_pmsbfa")) {
    quiet <- capture.output(baseline <- verbosity_fit(entry,controls))
    expect_length(quiet,0L)
    expect_true(baseline$converged,info=entry)
    expect_false(baseline$control$track_likelihood,info=entry)

    summary <- capture.output(fit <- verbosity_fit(entry,controls,verbose=TRUE))
    m <- fit$params$moddat
    expect_identical(summary,sprintf("G=%d q=%d r=%d %s BIC=%.2f (converged)",
      m$G,m$q,m$r,m$name,fit$bic),info=entry)
    expect_identical(fit$loglik,baseline$loglik,info=entry)
    expect_identical(fit$params,baseline$params,info=entry)

    tracked_controls <- modifyList(controls,list(track_likelihood=TRUE))
    iterations <- capture.output(tracked <- verbosity_fit(entry,tracked_controls))
    expect_length(iterations,tracked$iter)
    expect_true(all(grepl("^Iteration [0-9]+: log likelihood = -?[0-9]+\\.[0-9]{2}, change = -?[0-9]+\\.[0-9]{2}$",iterations)),info=entry)
    expect_identical(tracked$loglik,baseline$loglik,info=entry)
    expect_identical(tracked$params,baseline$params,info=entry)

    combined <- capture.output(both <- verbosity_fit(entry,tracked_controls,verbose=TRUE))
    expect_identical(combined,c(iterations,summary),info=entry)
    expect_identical(both$posterior,baseline$posterior,info=entry)
  }
})

test_that("unconverged fits are quiet and summaries respect the requested criterion", {
  controls <- quick_control(max.iter=0L)
  for(entry in c("pmbfa","pmsbfa","search_pmbfa","search_pmsbfa")) {
    output <- capture.output(fit <- verbosity_fit(entry,controls,verbose=TRUE))
    expect_false(fit$converged,info=entry)
    expect_length(output,0L)
  }
  controls <- quick_control(max.iter=3L,min.iter=1L,tol=1e6,min.component=0)
  output <- capture.output(fit <- verbosity_fit("pmsbfa",controls,
    verbose=TRUE,criterion="ICL"))
  expect_true(fit$converged)
  expect_identical(output,sprintf("G=2 q=1 r=1 UUU-UUU-U-U ICL=%.2f (converged)",fit$icl))
})

test_that("adaptive caches and multiple starts produce one summary per successful model", {
  d <- numerical_fixture("Normal",N=60L)
  controls <- quick_control(max.iter=3L,min.iter=1L,tol=1e6,min.component=0,nrep=2L)
  output <- capture.output(fit <- pmbfa(d$x,G=2L,q=1L,r=1L,
    inits=list(d$params,d$params),control=controls,verbose=TRUE))
  successful <- fit$model.criteria$model[fit$model.criteria$converged]
  expect_gt(length(successful),0L)
  expect_length(output,length(successful))
  expect_true(all(grepl("^G=2 q=1 r=1 [CU]{3}-[CU]{3} BIC=-?[0-9]+\\.[0-9]{2} \\(converged\\)$",output)))
  models <- vapply(strsplit(output," ",fixed=TRUE),function(parts)parts[4L],character(1))
  expect_setequal(models,successful)
  expect_equal(anyDuplicated(models),0L)
  expect_equal(fit$initialization$nstarts,2L)
  expect_identical(names(fit$selection$paths),c("forward","backward"))
})

test_that("initialization and short refinements stay quiet when tracking is requested", {
  controls <- quick_control(track_likelihood=TRUE,verbose=TRUE)
  for(family in c("Normal","ST")) {
    d <- numerical_fixture(family,N=60L)
    output <- capture.output(bank <- native("pmbfa_initializations")(
      d$x,2L,1L,1L,family,controls))
    expect_length(output,0L)
    expect_length(bank,1L)
  }
})

test_that("likelihood controls cannot partially match the legacy track alias", {
  d <- numerical_fixture("Normal",N=60L)
  controls <- quick_control(max.iter=3L,min.iter=1L,tol=1e6,min.component=0)
  output <- capture.output(fit <- pmbfa(d$x,G=2L,q=1L,r=1L,
    modelSubset="UUU-UUU",inits=d$params,control=controls,track_likelihood=TRUE))
  expect_true(fit$converged)
  expect_true(fit$control$track_likelihood)
  expect_false(fit$control$verbose)
  expect_length(output,fit$iter)
  expect_true(all(startsWith(output,"Iteration ")))
  output <- capture.output(legacy <- pmbfa(d$x,G=2L,q=1L,r=1L,
    modelSubset="UUU-UUU",inits=d$params,control=controls,track=TRUE))
  expect_true(legacy$control$verbose)
  expect_false(legacy$control$track_likelihood)
  expect_length(output,1L)
  expect_true(startsWith(output,"G=2 q=1 r=1 UUU-UUU BIC="))
})
