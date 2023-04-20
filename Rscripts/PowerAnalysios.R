# Power Analysis

# https://www.statmethods.net/stats/power.html
library(pwr)

#function	power calculations for
#pwr.2p.test	two proportions (equal n)
#pwr.2p2n.test	two proportions (unequal n)
#pwr.anova.test	balanced one way ANOVA
#pwr.chisq.test	chi-square test
#pwr.f2.test	general linear model
#pwr.p.test	proportion (one sample)
#pwr.r.test	correlation
#pwr.t.test	t-tests (one sample, 2 sample, paired)
#pwr.t2n.test	t-test (two samples with unequal n)

# Enter three of the four quantities (effect size, sample size, significance level, power) and the fourth is calculated

#Tests of Proportions

pwr.2p.test(h = , n = , sig.level =, power = )
# h=effect size, n=common sample size in each group
# h= 2arcsin(sqrt(p1))-2arcsin(sqrt(p2))
# Cohen suggests that h values of 0.2, 0.5, and 0.8 represent small, medium, and large effect sizes respectively.

#For unequal n's use
pwr.2p2n.test(h = , n1 = , n2 = , sig.level = , power = )

#you can specify alternative="two.sided", "less", or "greater" to indicate a two-tailed, or one-tailed test.(defalut=two-tailed)


# Chi-square Tests

pwr.chisq.test(w =, N = , df = , sig.level =, power = )

#where w is the effect size, N is the total sample size, and df is the degrees of freedom. The effect size w is defined as



# Test the power of group prportion change of genetic clusters

pwr.2p2n.test(h =0.2 , n1 = 58 , n2 =72 , sig.level = 0.05)
