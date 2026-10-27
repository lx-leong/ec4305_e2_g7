# Shared functions for the sharp RD analysis.

find_root <- function(start=getwd()) {
  root <- normalizePath(start,winslash='/',mustWork=TRUE)
  while (!file.exists(file.path(root,'data','nhis_clean_full.csv'))) {
    parent <- dirname(root)
    if (identical(parent,root)) stop('Cannot locate project data.')
    root <- parent
  }
  root
}

# Analysis sample: 1992-1996, doctor-visit outcome observed, positive weight.
load_rd_data <- function(root) {
  d <- read.csv(file.path(root,'data','nhis_clean_full.csv'))
  stopifnot(all(c('year','age4','post65','sawdr','inhosp','emp','mcare',
                  'covered2','female','wnh','bnh','hispanic','onh',
                  'lowed_clean','region','group','wt') %in% names(d)))
  d <- subset(d,year %in% 1992:1996 & !is.na(sawdr) & is.finite(wt) & wt>0)
  for(r in 1:4) d[[paste0('region_',r)]] <- ifelse(is.na(d$region),NA,as.integer(d$region==r))
  d
}

rd_labels <- c(sawdr='Saw a doctor (past year)',inhosp='Hospital stay (past year)',
               emp='Employed',mcare='Medicare coverage',covered2='Any coverage',
               female='Female',wnh='White non-Hispanic',bnh='Black non-Hispanic',
               hispanic='Hispanic',lowed_clean='Below high school',
               region_1='Northeast',region_2='Midwest',region_3='South',region_4='West',
               sawdr_index='Predicted doctor visit (covariate index)')

group_labels <- c('1'='White, high education','2'='White, low education',
                  '3'='Minority, high education','4'='Minority, low education')

# Local polynomial RD: separate slopes each side of the cutoff, triangular
# kernel times survey weights, survey-year fixed effects, HC1 standard errors.
# h: bandwidth in years (keeps |age - cutoff| < h).
# donut: drops cells with |age - cutoff| < donut (0.25 drops only the 65.00 cell).
# covars: optional numeric controls; predictions hold them at weighted means.
fit_rd <- function(d,y,h=5,cutoff=65,donut=0,poly=1,covars=NULL) {
  x <- d$age4-cutoff
  keep <- complete.cases(d[,c(y,'age4','year','wt',covars)]) &
    is.finite(d[[y]]) & is.finite(x) & is.finite(d$wt) & d$wt>0 &
    abs(x)<h & abs(x)>=donut
  z <- d[keep,]; z$x <- x[keep]; z$eligible <- as.integer(z$x>=0)
  z$y <- z[[y]]; z$w <- (1-abs(z$x)/h)*z$wt/mean(z$wt)
  trend <- if(poly==1) 'x' else paste0('(',paste(c('x',paste0('I(x^',2:poly,')')),collapse='+'),')')
  rhs <- c(paste0('eligible*',trend),'factor(year)',covars)
  f <- lm(as.formula(paste('y ~',paste(rhs,collapse='+'))),data=z,weights=w)
  X <- model.matrix(f); n <- nrow(X); k <- ncol(X)
  if(f$rank<k) stop('Rank-deficient fit: ',y)
  bread <- solve(crossprod(X,X*z$w))
  score <- X*as.numeric(z$w*residuals(f))
  V <- n/(n-k)*bread%*%crossprod(score)%*%bread
  tau <- unname(coef(f)['eligible']); se <- sqrt(V['eligible','eligible'])
  crit <- qt(.975,n-k); yr <- aggregate(w~year,z,sum)
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
    p_value=2*pt(-abs(tau/se),n-k),n_left=sum(z$x<0),n_right=sum(z$x>=0),
    below_rate=below,relative=tau/below)
  list(row=row,predict=pred,model=f)
}

rd_rows <- function(fits) do.call(rbind,lapply(fits,function(f) f$row))

pct <- function(v,digits=2) sprintf(paste0('%.',digits,'f'),100*v)

table_rd <- function(z,name,relative=TRUE) {
  tab <- data.frame(Specification=name,
    `Jump (pp)`=pct(z$estimate),
    `95% CI (pp)`=sprintf('[%s, %s]',pct(z$ci_lower),pct(z$ci_upper)),
    `p-value`=format.pval(z$p_value,digits=3,eps=.001),check.names=FALSE)
  if(relative) {
    tab$`Rate below cutoff` <- paste0(pct(z$below_rate,1),'%')
    tab$`Relative change` <- paste0(pct(z$relative,1),'%')
  }
  if('p_holm' %in% names(z)) tab$`Holm p-value` <- format.pval(z$p_holm,digits=3,eps=.001)
  tab$N <- format(z$n_left+z$n_right,big.mark=',',trim=TRUE)
  knitr::kable(tab,align=c('l',rep('r',ncol(tab)-1)),row.names=FALSE)
}

# Weighted quarter-age means with the local fit on each side of the cutoff.
plot_rd <- function(d,y,h=5,cutoff=65,main=rd_labels[y],fit=fit_rd(d,y,h,cutoff)) {
  z <- subset(d,is.finite(age4)&is.finite(wt)&wt>0&abs(age4-cutoff)<=h)
  z <- z[!is.na(z[[y]]),]
  g <- split(z,z$age4)
  a <- vapply(g,function(v) v$age4[1],numeric(1))
  m <- vapply(g,function(v) 100*weighted.mean(v[[y]],v$wt),numeric(1))
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
