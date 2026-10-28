# Shared functions for the sharp RD analysis.
#
# Data: NHIS 1992-1996 from the Card, Dobkin and Maestas replication file
# (data/c9203.sas7bdat). Survey weights are the supplied `wt`, unchanged.
# Missing outcomes are kept at load time; each regression drops only the
# variables it needs, so the doctor and hospital samples are outcome-specific.

find_root <- function(start=getwd()) {
  root <- normalizePath(start,winslash='/',mustWork=TRUE)
  while (!file.exists(file.path(root,'data','c9203.sas7bdat'))) {
    parent <- dirname(root)
    if (identical(parent,root)) stop('Cannot locate project data.')
    root <- parent
  }
  root
}

load_rd_data <- function(root) {
  d <- as.data.frame(haven::read_sas(file.path(root,'data','c9203.sas7bdat')))
  names(d) <- tolower(names(d))
  stopifnot(all(c('year','age','age4','post65','sawdr','inhosp','emp','mcare',
                  'covered2','female','wnh','bnh','hispanic','onh','educ',
                  'region','group','wt') %in% names(d)))
  d$.source_row <- seq_len(nrow(d))
  d <- subset(d,year %in% 1992:1996 & is.finite(age4) & is.finite(wt) & wt>0)
  stopifnot(all(d$post65==as.integer(d$age4>=65)))
  # Below high school; missing education stays missing.
  d$lowed_clean <- ifelse(is.na(d$educ),NA_real_,as.numeric(d$educ<12))
  # Race x education groups rebuilt from observed education. The supplied
  # `group` is kept for reference only.
  d$group_supplied <- d$group
  d$group <- ifelse(is.na(d$educ) | is.na(d$wnh),NA_integer_,
                    as.integer(1+2*(1-d$wnh)+d$lowed_clean))
  for(r in 1:4) d[[paste0('region_',r)]] <- as.numeric(d$region==r)
  for(v in c('sawdr','inhosp','mcare','covered2','emp','educ'))
    d[[paste0('observed_',v)]] <- as.numeric(!is.na(d[[v]]))
  d
}

rd_labels <- c(sawdr='Saw a doctor (past year)',inhosp='Hospital stay (past year)',
  emp='Employed (past two weeks)',mcare='Medicare coverage',covered2='Any health insurance',
  female='Female',wnh='White non-Hispanic',bnh='Black non-Hispanic',hispanic='Hispanic',
  lowed_clean='Below high school',region_1='Northeast',region_2='Midwest',
  region_3='South',region_4='West',
  observed_sawdr='Doctor-visit answer observed',observed_inhosp='Hospital-stay answer observed',
  observed_mcare='Medicare answer observed',observed_covered2='Coverage answer observed',
  observed_emp='Employment answer observed',
  observed_educ='Education observed',
  sawdr_index='Predicted doctor visit (covariate index)',
  inhosp_index='Predicted hospital stay (covariate index)')

group_labels <- c('1'='White non-Hispanic, high education','2'='White non-Hispanic, low education',
                  '3'='Other race/ethnicity, high education','4'='Other race/ethnicity, low education')
group_short <- c('1'='White, high ed','2'='White, low ed','3'='Other, high ed','4'='Other, low ed')

wmean <- function(v,w) {ok <- is.finite(v) & is.finite(w) & w>0; weighted.mean(v[ok],w[ok])}
pct <- function(v,digits=2) sprintf(paste0('%.',digits,'f'),100*v)
pv <- function(p) format.pval(p,digits=3,eps=.001)
ci <- function(r) sprintf('[%s, %s]',pct(r$ci_lower),pct(r$ci_upper))
fmt_n <- function(n) format(n,big.mark=',',trim=TRUE)
# "N left / N right" for one or more fitted rows.
n_lr <- function(z) sprintf('%s / %s',fmt_n(z$n_left),fmt_n(z$n_right))

# Log of every model's sample, exported for auditing.
model_samples <- list()

# Local polynomial RD: separate slopes each side of the cutoff, triangular
# kernel times the supplied survey weights, survey-year fixed effects, HC1
# standard errors and normal critical values (as in the fuzzy RD section).
# h: bandwidth in years (keeps |age - cutoff| < h).
# donut: drops cells with |age - cutoff| < donut (0.25 drops only the 65.00 cell;
#        0.5 drops 64.75, 65.00 and 65.25).
# covars: optional numeric controls; predictions hold them at weighted means.
fit_rd <- function(d,y,h=5,cutoff=65,donut=0,poly=1,covars=NULL,
                   weighted=TRUE,year_fe=TRUE,label='') {
  stopifnot(cutoff-h>=55.25,cutoff+h<=74.75)
  vars <- unique(c(y,'age4','year','wt',covars))
  x <- d$age4-cutoff
  complete <- complete.cases(d[,vars]) & is.finite(d[[y]]) & is.finite(x) &
    is.finite(d$wt) & d$wt>0
  window <- is.finite(x) & abs(x)<h & abs(x)>=donut
  z <- d[complete & window,,drop=FALSE]; z$x <- x[complete & window]
  z$eligible <- as.integer(z$x>=0); z$y <- z[[y]]
  sw <- if(weighted) z$wt else rep(1,nrow(z))
  z$w <- (1-abs(z$x)/h)*sw/mean(sw)
  trend <- if(poly==1) 'x' else paste0('(',paste(c('x',paste0('I(x^',2:poly,')')),collapse='+'),')')
  rhs <- c(paste0('eligible*',trend),if(year_fe) 'factor(year)',covars)
  f <- lm(as.formula(paste('y ~',paste(rhs,collapse='+'))),data=z,weights=w)
  X <- model.matrix(f); n <- nrow(X); k <- ncol(X)
  if(f$rank<k) stop('Rank-deficient fit: ',y)
  bread <- solve(crossprod(X,X*z$w))
  score <- X*as.numeric(z$w*residuals(f))
  V <- n/(n-k)*bread%*%crossprod(score)%*%bread
  tau <- unname(coef(f)['eligible']); se <- sqrt(V['eligible','eligible'])
  crit <- qnorm(.975); yr <- aggregate(w~year,z,sum)
  cm <- if(length(covars)) vapply(z[covars],weighted.mean,numeric(1),w=z$w) else NULL
  pred <- function(x,side=as.integer(x>=0)) {
    side <- rep(side,length.out=length(x))
    vapply(seq_along(x),function(j) {
      nd <- data.frame(x=x[j],eligible=side[j],year=yr$year)
      for(v in names(cm)) nd[[v]] <- cm[[v]]
      weighted.mean(predict(f,newdata=nd),yr$w)
    },numeric(1))
  }
  below <- pred(0,0)
  row <- data.frame(outcome=y,bandwidth=h,cutoff=cutoff,donut=donut,poly=poly,
    estimate=tau,se=se,ci_lower=tau-crit*se,ci_upper=tau+crit*se,
    p_value=2*pnorm(-abs(tau/se)),n_left=sum(z$x<0),n_right=sum(z$x>=0),N=n,
    below_rate=below,relative=if(abs(below)>1e-12) tau/below else NA_real_)
  model_samples[[length(model_samples)+1L]] <<- data.frame(
    label=if(nzchar(label)) label else y,outcome=y,required_variables=paste(vars,collapse=', '),
    h=h,cutoff=cutoff,donut=donut,polynomial=poly,year_FE=year_fe,survey_weighted=weighted,
    N_input=nrow(d),N_used=n,n_left=sum(z$x<0),n_right=sum(z$x>=0),
    age_min=min(z$age4),age_max=max(z$age4))
  list(row=row,predict=pred,model=f,data=z)
}

rd_rows <- function(fits) do.call(rbind,lapply(fits,function(f) f$row))

table_rd <- function(z,name,relative=TRUE) {
  tab <- data.frame(Specification=name,
    `Jump (pp)`=pct(z$estimate),
    `95% CI (pp)`=sprintf('[%s, %s]',pct(z$ci_lower),pct(z$ci_upper)),
    `p-value`=pv(z$p_value),check.names=FALSE)
  if(relative) {
    tab$`Rate below cutoff` <- paste0(pct(z$below_rate,1),'%')
    tab$`Relative change` <- paste0(pct(z$relative,1),'%')
  }
  if('p_holm' %in% names(z)) tab$`Holm p-value` <- pv(z$p_holm)
  tab$N <- fmt_n(z$N)
  tab$`N left` <- fmt_n(z$n_left)
  tab$`N right` <- fmt_n(z$n_right)
  knitr::kable(tab,align=c('l',rep('r',ncol(tab)-1)),row.names=FALSE)
}

# Weighted quarter-age means of the model's own sample, with the local fit
# on each side of the cutoff.
plot_rd <- function(d,y,h=5,cutoff=65,main=rd_labels[y],fit=fit_rd(d,y,h,cutoff)) {
  z <- fit$data
  g <- split(z,z$age4)
  a <- vapply(g,function(v) v$age4[1],numeric(1))
  m <- vapply(g,function(v) 100*wmean(v[[y]],v$wt),numeric(1))
  xl <- seq(-h+.01,-.001,length.out=100); xr <- seq(0,h-.01,length.out=100)
  yl <- 100*fit$predict(xl); yr <- 100*fit$predict(xr)
  plot(a,m,type='n',xlab='Age',ylab='Survey-weighted percentage',
       main=main,ylim=range(c(m,yl,yr)),bty='l',
       col.axis='#637783',col.lab='#344f60',col.main='#17465c')
  grid(nx=NA,ny=NULL,col='#edf1f3',lty=1)
  abline(v=cutoff,lty=2,col='#82949f')
  points(a,m,pch=21,bg='#96adb7',col='white',cex=1.15)
  lines(cutoff+xl,yl,col='#176d85',lwd=2.5)
  lines(cutoff+xr,yr,col='#c17948',lwd=2.5)
}

# Matched subgroup sample: group members with the outcome, any coverage,
# Medicare and employment all observed. Every subgroup regression for that
# outcome uses exactly these respondents.
subgroup_sample <- function(d,y,g)
  d[!is.na(d$group) & d$group==g & complete.cases(d[c(y,'covered2','mcare','emp')]),]

# First-stage and employment jumps by group on the matched samples. This does
# not estimate the outcome jump; that is left to the estimation report.
subgroup_first_stage <- function(d,y,h=5) {
  do.call(rbind,lapply(1:4,function(g) {
    dg <- subgroup_sample(d,y,g)
    fc <- fit_rd(dg,'covered2',h,label=paste('Subgroup',g,y,'coverage'))
    fm <- fit_rd(dg,'mcare',h,label=paste('Subgroup',g,y,'Medicare'))
    fe <- fit_rd(dg,'emp',h,label=paste('Subgroup',g,y,'employment'))
    stopifnot(identical(fc$data$.source_row,fm$data$.source_row),
              identical(fc$data$.source_row,fe$data$.source_row))
    pre <- subset(fc$data,age4>=63 & age4<65)
    c1 <- fc$row; m <- fm$row; e <- fe$row
    data.frame(sample=y,group=g,label=group_labels[[as.character(g)]],N=c1$N,
      n_left=c1$n_left,n_right=c1$n_right,
      pre_uninsured=1-wmean(pre$covered2,pre$wt),pre_mcare=wmean(pre$mcare,pre$wt),
      cov_jump=c1$estimate,cov_se=c1$se,cov_lo=c1$ci_lower,cov_hi=c1$ci_upper,cov_p=c1$p_value,
      mcare_jump=m$estimate,mcare_se=m$se,mcare_lo=m$ci_lower,mcare_hi=m$ci_upper,mcare_p=m$p_value,
      emp_jump=e$estimate,emp_se=e$se,emp_lo=e$ci_lower,emp_hi=e$ci_upper,emp_p=e$p_value)
  }))
}

export_models <- function(out) write.csv(do.call(rbind,model_samples),
  file.path(out,'model_sample_sizes.csv'),row.names=FALSE)
