# Negative binomial distribution and zero-inflated negative
# binomial distribution function,
# original scale and log-transformed.

# Use both parametrization with
# 1) r=number of failures, p=success rate
# 2) m=expected value, kappa=overdispersion
# (0=no overdispersion, increasing kappa means increasing overdispersion)

# Relationship between the two:
# m=pr/(1-p), kappa=1/r
# p=m/(m+1/kappa), r=1/kappa

# Zero-inflation only used in the m, kappa definition

# Note that R has a different definition than Wikipedia
# Wikipedia has the outcome to be number of successes
# before a given number of failures happen, while
# R has it as the number of failures before a
# given number of successes happen. Thus p.R=1-p.wikipedia
# and vice versa.
# Since my calculations are based on the Wikipedia
# definition, I will use the Wikipedia definition here also.


# First parametrization:
dnegbinom1=function(k,r,p,uselog=FALSE)
  dnbinom(k,size=r,prob=1-p,log=uselog)

# Second parametrization:
dnegbinom2=function(k,m,kappa,uselog=FALSE)
  dnegbinom1(k, 1/kappa, m/(m+1/kappa), uselog=uselog)
  


# Zero-inflated negative binomial (only second parametrization):

dnegbinom.zero=function(k,m,kappa,p.nonzero)
  p.nonzero*dnegbinom2(k,m,kappa)+(1-p.nonzero)*(k==0)

dlnegbinom.zero=function(k,m,kappa,p.nonzero)
  log(p.nonzero*dnegbinom2(k,m,kappa)+(1-p.nonzero)*(k==0))

