rm(list=ls())

library(sharp)

for (filename in list.files("Scripts/Functions/")){
  source(paste0("Scripts/Functions/", filename))
}

M=CoMembership(groups=sort(rep(1:3,5)))
sum(M)
plot(StabilityScore(selprop=M[upper.tri(M)], K=100, n_cat = 2, pi_list=seq(0.01,0.99,by=0.01)))
M=CoMembership(groups=sort(rep(1:5,3)))
sum(M)
lines(StabilityScore(selprop=M[upper.tri(M)], K=100, n_cat=2, pi_list=seq(0.01,0.99,by=0.01)), col="red")


M=CoMembership(groups=sort(rep(1:3,5)))
sum(M)
s=Score(M=M, L=3)
print(s)
plot(s)
M=CoMembership(groups=sort(rep(1:5,3)))
sum(M)
s=Score(M=M, L=5)
print(s)
lines(s, col="red")

pbinom(8, size=10, prob=0.2, lower.tail = FALSE)
pbinom(1, size=10, prob=0.8, lower.tail = TRUE)



pi=0.8
K=100
q=sum(M)/2
N=sum(upper.tri(M))
thr=round(pi*K)
p1=pbinom(thr-1, size=K, prob=q/N, lower.tail = FALSE, log.p = TRUE)
p0=pbinom(thr-1, size=K, prob=q/N, lower.tail = TRUE, log.p = TRUE)
(exp(p0)+exp(p1)-1)



