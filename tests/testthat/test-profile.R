test_that("profiling preserves all fitted numeric values and native printing precision", {
pmbfa_initializations <- native("pmbfa_initializations")
pmsbfaEM <- native("pmsbfaEM")
pmbfa_model_npar <- native("pmbfa_model_npar")
set.seed(491)
x <- array(rnorm(4*4*60),c(4,4,60))
ctl <- list(ncandid=6,nstart=2,nrep=1,init.iter=2,min.component=0,
            max.iter=5,tol=1e-12,profile=FALSE)
p <- pmbfa_initializations(x,2,1,1,'ST',ctl)[[1]]
a <- pmsbfaEM(x,p,ctl)
ctl$profile <- TRUE
b <- pmsbfaEM(x,p,ctl)
stopifnot(identical(a$loglik,b$loglik),identical(a$params,b$params),
          identical(a$posterior,b$posterior),is.null(a$diagnostics$timing),
          length(b$diagnostics$timing)==7L,all(b$diagnostics$timing>=0),
          sum(b$diagnostics$timing)<=b$diagnostics$total_native_seconds,
          pmbfa_model_npar(20,20,2,3,2,'CUC-CUC-C-U')==1302,
          pmbfa_model_npar(20,20,2,3,2,'CCC-CCC-C-U')==1300,
          b$npar==pmbfa_model_npar(4,4,2,1,1,p$moddat$name,'ST'))
ctl$verbose <- TRUE
printed <- capture.output(pmsbfaEM(x,p,ctl))
lines <- printed[grepl('^Iteration ',printed)]
stopifnot(length(lines)>0L,all(grepl('log likelihood = -?[0-9]+\\.[0-9]{2}, change = -?[0-9]+\\.[0-9]{2}$',lines)))
expect_length(b$diagnostics$pcg_iterations_by_cycle,b$iter)
expect_equal(sum(b$diagnostics$pcg_iterations_by_cycle),b$diagnostics$pcg_iterations)
cat('Profiling invariance, timing, native count, and native printing checks passed.\n')

expect_true(TRUE)
})
