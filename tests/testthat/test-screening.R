test_that("nested screening records decisions and respects its assumptions", {
e <- asNamespace("pmbfa")
checks <- 0L
check <- function(x,label) {checks <<- checks+1L;expect_true(isTRUE(x),info=label)}
fit <- function(name,bic,k,N=200,G=2,family="ST",ll=NULL,converged=TRUE) {
  if(is.null(ll))ll <- (k*log(N)-bic)/2
  structure(list(params=list(moddat=list(n=20,p=20,N=N,G=G,q=3,r=2,
       family=family,name=name)),loglik=ll,iter=1L,bic=bic,icl=bic+5,
       score=if(converged)bic else Inf,npar=k,converged=converged,
       status=if(converged)"converged"else"max_iter",message=""),
       class=c("pmsbfa","pmbfa"))
}
child <- "CCC-CCC-C-U"
parents <- list(fit("CUC-CCC-C-U",333678.9,1301),
                 fit("CCC-CUC-C-U",333699.3,1301))
rule <- e$.pmbfa_nested_screen(child,parents,333137,"BIC",function(name)1300)
check(rule$screened,"User example should screen the joint model")
check(abs(rule$lower.bound-(333699.3-log(200)))<1e-8,"Incorrect nested BIC lower bound")
check(333858.6>=rule$lower.bound,"True child BIC should exceed the bound in example")
check(e$.pmbfa_nested_screen(child,parents,333142,"ICL",function(name)1300)$screened,
       "BIC lower bound may also screen ICL")
failed <- parents;failed[[1]]$converged <- FALSE
check(!e$.pmbfa_nested_screen(child,failed,333137,"BIC",function(name)1300)$screened,
       "Unconverged parents cannot screen")
mismatch <- parents;mismatch[[2]]$params$moddat$q <- 4
check(!e$.pmbfa_nested_screen(child,mismatch,333137,"BIC",function(name)1300)$screened,
       "Mismatched parent tuples cannot screen")
check(!e$.pmbfa_nested_screen("CCC-CCC-U-U",parents,333137,"BIC",function(name)1300)$screened,
       "Non-nested parents cannot screen")
check(!e$.pmbfa_nested_screen(child,parents,Inf,"BIC",function(name)1300)$screened,
       "No finite incumbent means no screening")
check(!e$.pmbfa_nested_screen(child,parents,rule$lower.bound,"BIC",function(name)1300)$screened,
       "A tie must not screen")
# A finite-start parent local maximum is not an upper bound: the unseen child
# can have a better attained likelihood than that parent. This deliberately
# demonstrates why nested screening must not be described as a certificate.
child_local_ll <- min(vapply(parents,e$.pmbfa_final_ll,numeric(1)))+300
child_local_bic <- -2*child_local_ll+1300*log(200)
check(child_local_bic<333137 && rule$screened,
       "Local-optimum counterexample must remain visible in the audit")
check(identical(eval(formals(e$pmsbfa)$screening),c("none","nested")),
       "Unscreened search must remain the default")

run_path <- function(screening,criterion="BIC") {
  cache <- new.env(parent=emptyenv());calls <- character()
  npar <- function(name)100+sum(strsplit(name,"",fixed=TRUE)[[1]] == "U")
  fitter <- function(name) {
    if(exists(name,cache,inherits=FALSE))return(get(name,cache))
    calls <<- c(calls,name)
    count <- sum(strsplit(name,"",fixed=TRUE)[[1]] == "C")
    ll <- 10000-1000*count;k<-npar(name);bic <- -2*ll+k*log(60)
    out <- fit(name,bic,k,N=60,family="Normal",ll=ll)
    out$score <- out[[tolower(criterion)]]
    assign(name,out,cache);out
  }
  result <- e$.pmbfa_adaptive(fitter,2,"Normal","backward",criterion,
    screening=screening,npar_model=npar,
    has_fit=function(name)exists(name,cache,inherits=FALSE))
  list(result=result,calls=calls)
}
none <- run_path("none");nested <- run_path("nested")
check(identical(none$result$params$moddat$name,nested$result$params$moddat$name),
       "Screened and unscreened exact-optimum mock selections differ")
check(length(none$calls)==10L && length(nested$calls)==7L,
       "Three joint fits should be saved in backward mock path")
check(nrow(nested$result$selection$screened.models)==3L,
       "Every skipped joint model needs an explicit record")
for(stage in nested$result$selection$paths$backward) {
  skipped <- stage[stage$screened,,drop=FALSE]
  check(nrow(skipped)==1L && all(is.na(skipped$score)) &&
        all(is.na(skipped$converged)) && skipped$status=="screened",
        "Skipped model must not receive fabricated fit quantities")
  check(!any(skipped$model %in% nested$calls),"A skipped model was nevertheless fitted")
}
icl <- run_path("nested","ICL")
check(nrow(icl$result$selection$screened.models)==3L,"ICL screening path failed")
# Existing cached fits are used, not converted into skipped placeholders.
cache_case <- e$.pmbfa_adaptive(function(name) {
  out<-fit(name,100,10,N=60,family="Normal");out$score<-100;out
},2,"Normal","backward","BIC",screening="nested",
 npar_model=function(name)1,has_fit=function(name)TRUE)
check(nrow(cache_case$selection$screened.models)==0L,"Cached models must not be screened")

printed <- capture.output(e$print.pmbfa(fit(child,333858.612345,1300)))
check(any(grepl("333858.61",printed,fixed=TRUE)),"Fit BIC must print two decimals")
check(!any(grepl("333858.612",printed,fixed=TRUE)),"Fit printed past two decimals")
# Methods use numeric formatting explicitly even when global print options vary.
printed_grid <- capture.output(e$print.pmbfa_fits(list(fit(child,333858.612345,1300))))
check(any(grepl("333858.61",printed_grid,fixed=TRUE)) &&
       !any(grepl("333858.612",printed_grid,fixed=TRUE)),"Grid criterion must print two decimals")
cat("Screening and print checks:",checks,"passed.\n")
cat("Finite-start counterexample: screened bound",rule$lower.bound,
    "but possible child local BIC",child_local_bic,"beats incumbent333137.\n")


})
