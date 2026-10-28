# Revised copy; the original rd_functions.R is unchanged.
# Shared functions for the revised sharp-RD reports. No generated outputs are inputs.
find_root <- function(start=getwd()) {
  root <- normalizePath(start,winslash='/',mustWork=TRUE)
  while (!file.exists(file.path(root,'data','c9203.sas7bdat'))) {
    parent <- dirname(root)
    if (identical(parent,root)) stop('Cannot locate project data.')
    root <- parent
  }
  root
}

# Restore the weight fields recommended in the 1993 and 1996 NHIS insurance
# codebooks. combo1.sas divided BOTH by two. Preserve the supplied wt explicitly.
# 1993 is a second-half-year population weight, NOT a full-year observation.
nhis_weights <- function(d,scheme='nhis') {
  if (!'wt_replication' %in% names(d)) d$wt_replication <- d$wt
  stopifnot(all(d$year %in% 1992:1996))
  d$weight_restore_factor <- ifelse(d$year %in% c(1993,1996),2,1)
  d$wt_nhis <- d$wt_replication*d$weight_restore_factor
  d$wt <- if(scheme=='nhis') d$wt_nhis else d$wt_replication
  d$wt_pool <- d$wt/5
  d$weight_scheme <- scheme
  d
}

# Retain missing responses here; each model drops only its required variables.
load_rd_data <- function(root) {
  d <- as.data.frame(haven::read_sas(file.path(root,'data','c9203.sas7bdat')))
  names(d) <- tolower(names(d)); d$.source_row <- seq_len(nrow(d))
  d <- subset(d,year %in% 1992:1996 & is.finite(age4) & is.finite(wt) & wt>0)
  d <- nhis_weights(d)
  d$lowed_clean <- ifelse(is.na(d$educ),NA_real_,as.numeric(d$educ<12))
  d$group_original <- d$group
  d$group <- ifelse(is.na(d$educ)|is.na(d$wnh),NA_integer_,
                    as.integer(1+2*(1-d$wnh)+d$lowed_clean))
  for(r in 1:4) d[[paste0('region_',r)]] <- as.numeric(d$region==r)
  for(v in c('sawdr','inhosp','mcare','educ')) d[[paste0('observed_',v)]] <- as.numeric(!is.na(d[[v]]))
  d$included_doctor_iv <- as.numeric(!is.na(d$sawdr)&!is.na(d$mcare))
  stopifnot(all(d$post65==as.integer(d$age4>=65)))
  d
}
rd_labels <- c(sawdr='Doctor contact (past year)',inhosp='Hospital stay (past year)',
  emp='Employed (past two weeks)',mcare='Medicare coverage',covered2='Any coverage',
  female='Female',wnh='White non-Hispanic',bnh='Black non-Hispanic',hispanic='Hispanic',
  lowed_clean='Below high school',region_1='Northeast',region_2='Midwest',region_3='South',
  observed_sawdr='Doctor response observed',observed_inhosp='Hospital response observed',
  observed_mcare='Medicare response observed',observed_educ='Education observed',
  included_doctor_iv='Doctor and Medicare both observed')
group_labels <- c('1'='White non-Hispanic, high education','2'='White non-Hispanic, low education',
                  '3'='Other race/ethnicity, high education','4'='Other race/ethnicity, low education')
model_samples <- list()
wmean <- function(v,w) {ok <- is.finite(v)&is.finite(w)&w>0; weighted.mean(v[ok],w[ok])}
pct <- function(v,digits=2) sprintf(paste0('%.',digits,'f'),100*v)
pv <- function(p) format.pval(p,digits=3,eps=.001)
ci <- function(r) sprintf('[%s, %s]',pct(r$ci_lower),pct(r$ci_upper))

# Separate local polynomial slopes; year FE; triangular x survey weights; HC1.
# A donut of .25 removes age 65.00; .5 removes 64.75,65.00,65.25.
fit_rd <- function(d,y,h=5,cutoff=65,donut=0,poly=1,covars=NULL,
                   weighted=TRUE,year_fe=TRUE,label='') {
  stopifnot(cutoff-h>=55.25,cutoff+h<=74.75)
  vars <- unique(c(y,'age4','year','wt_pool',covars))
  x <- d$age4-cutoff
  complete <- complete.cases(d[,vars]) & is.finite(d[[y]]) & is.finite(x) &
    is.finite(d$wt_pool)&d$wt_pool>0
  window <- is.finite(x)&abs(x)<h&abs(x)>=donut
  z <- d[complete&window,,drop=FALSE]; z$x <- x[complete&window]
  z$eligible <- as.integer(z$x>=0); z$y <- z[[y]]
  sw <- if(weighted) z$wt_pool else rep(1,nrow(z))
  z$w <- (1-abs(z$x)/h)*sw/mean(sw)
  trend <- if(poly==1) 'x' else '(x+I(x^2))'
  rhs <- c(paste0('eligible*',trend),if(year_fe)'factor(year)',covars)
  f <- lm(as.formula(paste('y ~',paste(rhs,collapse='+'))),data=z,weights=w)
  X <- model.matrix(f); n <- nrow(X); k <- ncol(X)
  if(f$rank<k) stop('Rank-deficient fit: ',y)
  stopifnot(n>k,all(z$year %in% 1992:1996),all(z$age4>=55.25&z$age4<=74.75))
  bread <- solve(crossprod(X,X*z$w))
  score <- X*as.numeric(z$w*residuals(f))
  V <- n/(n-k)*bread%*%crossprod(score)%*%bread
  tau <- unname(coef(f)['eligible']); se <- sqrt(V['eligible','eligible'])
  if(length(unique(z$y))==1) {tau<-0;se<-0}
  p <- if(se==0&&tau==0) 1 else 2*pnorm(-abs(tau/se))
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
    estimate=tau,se=se,ci_lower=tau-crit*se,ci_upper=tau+crit*se,p_value=p,
    n_left=sum(z$x<0),n_right=sum(z$x>=0),N=n,below_rate=below,
    relative=if(abs(below)>1e-12)tau/below else NA_real_,year_FE=year_fe,
    survey_weighted=weighted,weight_scheme=paste(unique(z$weight_scheme),collapse=';'))
  model_samples[[length(model_samples)+1L]] <<- data.frame(
    label=if(nzchar(label))label else y,outcome=y,required_variables=paste(vars,collapse=', '),
    h=h,cutoff=cutoff,donut=donut,polynomial=poly,year_FE=year_fe,
    weight_scheme=if(weighted)paste(unique(z$weight_scheme),collapse=';')else'kernel only',
    N_input=nrow(d),N_complete_before_window=sum(complete),
    N_window_before_variable_drops=sum(window),N_used=n,n_left=sum(z$x<0),n_right=sum(z$x>=0),
    age_min=min(z$age4),age_max=max(z$age4))
  list(row=row,predict=pred,model=f,data=z)
}
rd_rows <- function(fits) do.call(rbind,lapply(fits,function(f)f$row))
table_rd <- function(z,name,relative=TRUE) {
  tab <- data.frame(Specification=name,`Jump (pp)`=pct(z$estimate),
    `95% CI (pp)`=sprintf('[%s, %s]',pct(z$ci_lower),pct(z$ci_upper)),
    `p-value`=pv(z$p_value),check.names=FALSE)
  if(relative) {
    tab$`Below-cutoff rate` <- paste0(pct(z$below_rate,1),'%')
    tab$`Relative change` <- paste0(pct(z$relative,1),'%')
  }
  if('p_holm'%in%names(z)) tab$`Holm p` <- pv(z$p_holm)
  tab$N <- z$N; tab$`N left / right` <- paste(z$n_left,z$n_right,sep=' / ')
  knitr::kable(tab,row.names=FALSE)
}
plot_rd <- function(d,y,h=5,cutoff=65,main=rd_labels[y],fit=fit_rd(d,y,h,cutoff)) {
  z <- fit$data # Plot exactly the model's complete-case sample, including donuts.
  g <- split(z,z$age4)
  a <- vapply(g,function(v)v$age4[1],numeric(1))
  m <- vapply(g,function(v)100*wmean(v[[y]],v$wt),numeric(1))
  xl <- seq(-h+.01,-.001,length.out=100); xr <- seq(0,h-.01,length.out=100)
  yl <- 100*fit$predict(xl); yr <- 100*fit$predict(xr)
  plot(a,m,type='n',xlab='Quarter-based age',ylab='Survey-weighted percentage',
       main=main,ylim=range(c(m,yl,yr)),bty='l',col.main='#17465c')
  grid(nx=NA,ny=NULL,col='#edf1f3',lty=1); abline(v=cutoff,lty=2,col='#82949f')
  points(a,m,pch=21,bg='#96adb7',col='white',cex=1.1)
  lines(cutoff+xl,yl,col='#176d85',lwd=2.5);lines(cutoff+xr,yr,col='#c17948',lwd=2.5)
}

# Mechanism comparisons use a common sample within each outcome and group.
subgroup_analysis <- function(d,y,h=5) {
  do.call(rbind,lapply(1:4,function(g) {
    dg <- d[!is.na(d$group)&d$group==g&complete.cases(d[c(y,'mcare','emp')]),]
    fy <- fit_rd(dg,y,h,label=paste('Matched subgroup',g,y))
    fm <- fit_rd(dg,'mcare',h,label=paste('Matched subgroup',g,y,'Medicare'))
    fe <- fit_rd(dg,'emp',h,label=paste('Matched subgroup',g,y,'employment'))
    stopifnot(identical(fy$data$.source_row,fm$data$.source_row),
              identical(fy$data$.source_row,fe$data$.source_row))
    pre <- subset(fy$data,age4>=63 & age4<65)
    a<-fy$row;m<-fm$row;e<-fe$row
    data.frame(outcome=y,group=g,label=group_labels[as.character(g)],N=a$N,
      n_left=a$n_left,n_right=a$n_right,N_mcare=m$N,N_emp=e$N,
      estimate=a$estimate,se=a$se,lo=a$ci_lower,hi=a$ci_upper,p=a$p_value,
      mcare_jump=m$estimate,mcare_lo=m$ci_lower,mcare_hi=m$ci_upper,
      emp_jump=e$estimate,emp_lo=e$ci_lower,emp_hi=e$ci_upper,emp_p=e$p_value,
      pre_uninsured=1-wmean(pre$covered2,pre$wt),N_pre_insurance=sum(!is.na(pre$covered2)))
  }))
}
weight_audit <- function(d) do.call(rbind,lapply(split(d,d$year),function(v)
  data.frame(year=v$year[1],N=nrow(v),restore_factor=v$weight_restore_factor[1],
    supplied_weight_sum=sum(v$wt_replication),source_scale_sum=sum(v$wt_nhis),
    pooled_weight_sum=sum(v$wt_pool))))
export_models <- function(out) write.csv(do.call(rbind,model_samples),
                                         file.path(out,'model_sample_sizes.csv'),row.names=FALSE)
