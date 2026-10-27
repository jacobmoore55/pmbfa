verbosity_fit <- function(entry, control, verbose=FALSE, criterion="BIC",quiet=NULL) {
  skewed <- grepl("pmsbfa",entry,fixed=TRUE)
  d <- numerical_fixture(if(skewed) "ST" else "Normal",N=60L)
  arguments <- list(x=d$x,G=2L,inits=d$params,control=control,
    modelSubset=d$params$moddat$name,verbose=verbose,criterion=criterion)
  if(!is.null(quiet))arguments$quiet <- quiet
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

test_that("adaptive constraints and multiple starts produce only the tuple winner", {
  d <- numerical_fixture("Normal",N=60L)
  controls <- quick_control(max.iter=3L,min.iter=1L,tol=1e6,min.component=0,nrep=2L)
  output <- capture.output(fit <- pmbfa(d$x,G=2L,q=1L,r=1L,
    inits=list(d$params,d$params),control=controls,verbose=TRUE))
  expect_gt(sum(fit$model.criteria$converged),1L)
  expect_identical(output,sprintf("G=2 q=1 r=1 %s BIC=%.2f (converged)",
    fit$params$moddat$name,fit$bic))
  expect_equal(fit$bic,min(fit$model.criteria$score))
  expect_equal(fit$initialization$nstarts,2L)
  expect_identical(names(fit$selection$paths),c("forward","backward"))
})

verbosity_summary <- function(fit) {
  m <- fit$params$moddat
  sprintf("G=%d q=%d r=%d %s %s=%.2f (converged)",
    m$G,m$q,m$r,m$name,fit$criterion,fit$score)
}

test_that("direct fitting reports one selected model for every tested G,q,r tuple", {
  controls <- quick_control(max.iter=3L,min.iter=1L,tol=1e6,min.component=0)
  G <- c(2L,1L)
  tuples <- expand.grid(G=G,q=1:2,r=1:2)
  expected_keys <- with(tuples,sprintf("G=%d q=%d r=%d",G,q,r))
  for(entry in c("pmbfa","pmsbfa")) {
    skewed <- entry == "pmsbfa"
    d <- numerical_fixture(if(skewed) "ST" else "Normal",n=6L,p=6L,N=60L)
    subset <- if(skewed)c("CCC-CCC-C-C","UUU-UUU-U-U") else c("CCC-CCC","UUU-UUU")
    output <- capture.output(fits <- do.call(getExportedValue("pmbfa",entry),
      list(x=d$x,G=G,q=1:2,r=1:2,modelSubset=subset,control=controls,verbose=TRUE)))
    expect_true(all(vapply(fits,function(fit)fit$converged,logical(1))),info=entry)
    expect_length(output,nrow(tuples))
    expect_setequal(output,vapply(fits,verbosity_summary,character(1)))
    keys <- vapply(strsplit(output," ",fixed=TRUE),
      function(parts)paste(parts[1:3],collapse=" "),character(1))
    expect_identical(keys,expected_keys,info=entry)
    for(fit in fits) expect_equal(fit$score,min(fit$model.criteria$score))
  }
})

test_that("factor searches report only the best model for each tested G", {
  controls <- quick_control(max.iter=3L,min.iter=1L,tol=1e6,min.component=0)
  G <- c(2L,1L)
  for(entry in c("search_pmbfa","search_pmsbfa")) {
    skewed <- entry == "search_pmsbfa"
    d <- numerical_fixture(if(skewed) "ST" else "Normal",n=6L,p=6L,N=60L)
    subset <- if(skewed)c("CCC-CCC-C-C","UUU-UUU-U-U") else c("CCC-CCC","UUU-UUU")
    for(method in c("adaptive","grid")) {
      arguments <- list(x=d$x,G=G,max_q=2L,max_r=2L,modelSubset=subset,
        control=controls,verbose=TRUE)
      if(method == "grid")arguments$search.method <- method
      output <- capture.output(fits <- do.call(getExportedValue("pmbfa",entry),arguments))
      expect_true(all(vapply(fits,function(fit)fit$converged,logical(1))),info=entry)
      expect_length(output,length(G))
      by_G <- fits[match(G,vapply(fits,function(fit)as.integer(fit$params$moddat$G),integer(1)))]
      expect_identical(output,unname(vapply(by_G,verbosity_summary,character(1))),info=entry)
      for(fit in fits) {
        expect_gt(nrow(fit$factor.search),1L)
        expect_equal(fit$score,min(fit$factor.search$score))
        expect_identical(fit$factor.search.method,method)
        expect_true(fit$control$verbose)
      }
    }
  }
})

test_that("parallel fitting reports tuple winners in the parent process", {
  d <- numerical_fixture("Normal",N=60L)
  controls <- quick_control(max.iter=3L,min.iter=1L,tol=1e6,min.component=0)
  arguments <- list(x=d$x,G=c(2L,1L),q=1L,r=1L,
    modelSubset=c("CCC-CCC","UUU-UUU"),control=controls,verbose=TRUE)
  serial <- capture.output(expected <- do.call(pmbfa,arguments))
  expect_length(serial,2L)
  arguments$parallelize <- TRUE
  arguments$ncores <- 2L
  if(.Platform$OS.type == "unix") {
    output <- capture.output(forked <- do.call(pmbfa,arguments))
    expect_identical(output,serial)
    expect_identical(lapply(forked,function(fit)fit$params),
      lapply(expected,function(fit)fit$params))
  }
  engine <- clone_engine()
  engine$.Platform <- modifyList(.Platform,list(OS.type="windows"))
  output <- capture.output(socket <- do.call(engine$pmbfa,arguments))
  expect_identical(output,serial)
  expect_identical(lapply(socket,function(fit)fit$params),
    lapply(expected,function(fit)fit$params))
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

test_that("quiet is the inverse verbosity alias and overrides conflicting flags", {
  controls <- quick_control(max.iter=3L,min.iter=1L,tol=1e6,min.component=0)
  for(entry in c("pmbfa","pmsbfa","search_pmbfa","search_pmsbfa")) {
    reference <- capture.output(expected <- verbosity_fit(entry,controls,verbose=TRUE))
    output <- capture.output(fit <- verbosity_fit(entry,
      modifyList(controls,list(verbose=FALSE)),quiet=FALSE))
    expect_identical(output,reference,info=entry)
    expect_length(output,1L)
    expect_true(fit$control$verbose,info=entry)
    expect_identical(fit$params,expected$params,info=entry)
    expect_identical(fit$loglik,expected$loglik,info=entry)
    output <- capture.output(silent <- verbosity_fit(entry,
      modifyList(controls,list(verbose=TRUE)),verbose=TRUE,quiet=TRUE))
    expect_length(output,0L)
    expect_false(silent$control$verbose,info=entry)
    expect_identical(silent$params,expected$params,info=entry)
    for(invalid in list(NA,c(FALSE,TRUE),0))
      expect_error(verbosity_fit(entry,controls,quiet=invalid),"quiet must be TRUE or FALSE")
  }
})

test_that("quiet=FALSE produces exactly one search winner for each G in 1:4", {
  controls <- quick_control(max.iter=3L,min.iter=1L,tol=1e6,min.component=0)
  d <- numerical_fixture("Normal",N=120L,G=4L)
  for(entry in c("search_pmbfa","search_pmsbfa")) {
    subset <- if(entry == "search_pmbfa")c("CCC-CCC","UUU-UUU") else
      c("CCC-CCC-C-C","UUU-UUU-U-U")
    output <- capture.output(fits <- do.call(getExportedValue("pmbfa",entry),
      list(x=d$x,G=1:4,max_q=1L,max_r=1L,modelSubset=subset,
        control=controls,quiet=FALSE)))
    expect_true(all(vapply(fits,function(fit)fit$converged,logical(1))),info=entry)
    expect_length(output,4L)
    by_G <- fits[order(vapply(fits,function(fit)fit$params$moddat$G,numeric(1)))]
    expect_identical(output,unname(vapply(by_G,verbosity_summary,character(1))),info=entry)
  }
})
