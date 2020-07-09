### R code associated with Vannatta and Minchella. THE CRYPTIC INFLUENCE OF PARASITISM ON ECOSYSTEM PROCESSES

###set directory

setwd("F:/A_in_prep/Microcosms/Data")


###import data
mydata <- read.csv("mesocosm.master.csv", header = TRUE)


### mean infection intensity in species
phymet <- lm(log(mydata$physa.mean.metas + 1) ~ mydata$infected.days)
summary(phymet)
plot(phymet)
shapiro.test(phymet$residuals)

promet <- lm(log(mydata$pro.mean.metas + 1) ~ mydata$infected.days)
summary(promet)
plot(promet)
shapiro.test(promet$residuals)



###prevalence functions
phyprev <- lm(mydata$physa.prev ~ mydata$infected.days)
summary(phyprev)

proprev <- lm(mydata$pro.prev ~ mydata$infected.days)
summary(proprev)

# cannot responsibly use linear regression since predicts beyond the bounds of 100%
#create nonlinear function
fx <- function(x, k) 1*(1-(exp(-k*x)))
plot(mydata$physa.prev ~ mydata$infected.days)
curve(fx(x, k=0.1), add=T)

nlsfit <- nls(physa.prev ~ fx(infected.days, k),
              data=mydata,
              start=list(k=0.1))
p <- coef(nlsfit)

plot(mydata$physa.prev ~ mydata$infected.days, ylab = 'Infection prevalence in Physa', xlab = 'Infected-snail days')
curve(fx(x, k=p["k"]), 
      add=T, col="red")

#does linear or asmyptotic work best
library(pgirmess)
out <- selMod(list(nlsfit, phyprev), Order="AICc")
out

#now do for promenetus
plot(mydata$pro.prev ~ mydata$infected.days)
curve(fx(x, k=0.01), add=T)

nls.pro <- nls(pro.prev ~ fx(infected.days, k),
               data=mydata,
               start=list(k=0.01))
p <- coef(nls.pro)

plot(mydata$pro.prev ~ mydata$infected.days, ylab = 'Infection prevalence in Promenetus', xlab = 'Infected-snail days')
curve(fx(x, k=p["k"]), 
      add=T, col="red")
abline(proprev)

out <- selMod(list(nls.pro, proprev), Order="AICc")
out



### look at snail abundance across the treatment sets
physa <- lm(log(mydata$physa.abun) ~ mydata$infected.days)
summary(physa)
shapiro.test(physa$residuals)

promen <- lm(mydata$pro.abun ~ mydata$infected.days)
summary(promen)
shapiro.test(promen$residuals)

total.ab <- lm(log(mydata$total.abun) ~ mydata$infected.days)
summary(total.ab)
shapiro.test(total.ab$residuals)

phys.pro <- lm(log(mydata$phys.pro) ~ mydata$infected.days)
summary(phys.pro)
shapiro.test(phys.pro$residuals)

phys.all <- lm((mydata$physa.abun/mydata$total.abun) ~ mydata$infected.days)
summary(phys.all)
shapiro.test(phys.all$residuals)


### periphyton characteristics
#create quadratic value for infected-snail days
infected2 <- mydata$infected.days * mydata$infected.days

#periphyton percent ash-free dry mass
AFDMq <-lm(mydata$X.AFDM ~ mydata$infected.days + infected2)
summary(AFDMq)
plot(AFDMq)
shapiro.test(AFDMq$residuals)

### chlorophyll a
chloro <- lm(mydata$chla.per.cm2 ~ mydata$infected.days)
summary(chloro)
plot(chloro)
shapiro.test(chloro$residuals)

#periphyton dry mass per square centimeter
periq <-lm(log(mydata$dm.per.cm2) ~ mydata$infected.days + infected2)
summary(periq)
shapiro.test(periq$residuals)
plot(periq)



###Wolfiella biomass
wolf.mod <- lm(mydata$total.wolf ~ mydata$infected.days)
summary(wolf.mod)
plot(wolf.mod)
shapiro.test(wolf.mod$residuals)


###in situ diel primary production
#week 6
DO <- lm(log(mydata$do.delta.s) ~ mydata$infected.days)
summary(DO)
plot(mydata$do.delta.s ~ mydata$infected.days)
shapiro.test(DO$residuals)



### nutrients in the mesocosms
TN<- read.csv("nutrients.csv", header = TRUE)

#nutrients for week 6 and later (once parasitism could have impact)
nut.late <- TN[25:97,]

### DOC mixed model
library(lattice)
library(lme4)
library(blmeco)
library(car)
library(effects)
M2 <- lmer(DOC ~ infected.days + (1 | tank), data = nut.late)
summary(M2)
# for reasonable sample sizes, can use Wald inference:
Betas  <- fixef(M2)                  #Get the betas
SE     <-  sqrt(diag(vcov(M2)))      #Get the SEs
pval   <- 2*pnorm(-abs(Betas  / SE)) #Z distribution
Output <- cbind(Betas,SE, pval)
print(Output, digits = 3) # P values for fixed effect
plot(M2)
qqmath(M2)
Anova(M2, type="III")

#total dissolved nitrogen dynamics
N.mod <- lm(log(TN$TN) ~ TN$infected.days * TN$week)
summary(N.mod)
plot(N.mod)
shapiro.test(N.mod$residuals)


#### Additional analyses within supplementary material ####

### invertebrate community NMDS
library(vegan)
data <- read.csv("community.csv")

#rare species occuring in only one tank or ubiquitous species are removed as these provide little information
#or spurious correlations. Cladocerans are analysed at the genus level.
zoos <- metaMDS(data[,c(6,9:12,14:15,17:23)], distance = "jaccard" , k=2, trymax=100)		## Ordinate in 2 dim, 20 tries for best solution
plot(zoos, type = 't', xlim = c(-1.75,1))
plot(zoos$points)
points(zoos, pch = 16, cex = 2, col=as.numeric(data$col + 2))
stressplot(zoos)	
#prep data for adonis
treat <- as.factor(data$col)
species <- data[,c(6,9:12,14:15,17:23)]
# Do adonis to check for differences
adon <- adonis(species ~ treat, method="jaccard", data=data, control=permControl(strata=treat), permuations=9999)
adon

### TDP model
P.mod <- lm(log(TN$TP) ~ TN$infected.days * TN$week)
summary(P.mod)
plot(P.mod)
shapiro.test(P.mod$residuals)


### look at stoichmetry
CP.mod <- lm(log(TN$CP) ~ TN$infected.days * TN$week)
summary(CP.mod)
plot(CP.mod)
shapiro.test(CP.mod$residuals)

CN.mod <- lm(TN$CN ~ TN$infected.days * TN$week)
summary(CN.mod)
plot(CN.mod)
shapiro.test(CN.mod$residuals)


NP.mod <- lm(TN$NP ~ TN$infected.days * TN$week)
summary(NP.mod)
plot(NP.mod)
shapiro.test(NP.mod$residuals)

#N:P stoich residuals are not normal. Attempt a box-cox transformation
my.bc <- boxcox(NP.mod, lambda=seq(-0.5, 0.5, by=0.01))
max(my.bc$y)
my.bc$x[which(my.bc$y==(max(my.bc$y)))]

#suggests raising NP data to -0.04 will be closest to normal.
NP.mod <- lm((TN$NP^ -0.04) ~ TN$infected.days * TN$week)
summary(NP.mod)
plot(NP.mod)
shapiro.test(NP.mod$residuals)

#cannot get NP to normal so use nonparametric ancova from sm pacakge; but have to treat treatment as a factor
library(sm)
np.NP.mod <- sm.ancova( TN$week, TN$NP, TN$factor, model='equal')









