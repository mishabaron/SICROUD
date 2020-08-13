#################################################################
##### SICROUD model #############################################
# This is actually SICROUD, focusing on Ttheta ##################
# No plot version, just for return()           ##################
#################################################################

SICROUD = function(Country){  

# theta = daily infection rate per person, we drop the Bern or Poi meet part of it
# h = disease duration. Nobody should stay more than h days in state C.
#     After h days, confirmed live cases will be considered U = unobserved recoveries

h = 14; 

# Choose the moving average window:
MAh = 7; 

# Choose the minimum number of cases to perform statistical inference:
MinCases = 30;
SmallCases = 1;

# Indices for the same country are sometimes different in some countries (ex.: France)
CountryIndexC = which(ConfirmedGlobal$Country.Region == Country);
CountryIndexR = which(RecoveredGlobal$Country.Region == Country);
CountryIndexD = which(DiedGlobal$Country.Region == Country);

L = dim(ConfirmedGlobal)[2];   # Number of days. Data start from column 5 (01-22-2020)
ConfirmedMatrix = ConfirmedGlobal[CountryIndexC,]
DiedMatrix      = DiedGlobal[CountryIndexD,]
RecoveredMatrix = RecoveredGlobal[CountryIndexR,]

# Organize observed processes for the given country
# C, O, D = CUMULATIVE confirmed, recovered-observed, died total @ day t
# c, o, d = NEW confirmed, recovered-observed, died @ day t ( = dC, dR, dD)

if (nrow(ConfirmedMatrix)==1){ 
 C = as.numeric(ConfirmedMatrix[,5:L]); 
 O = as.numeric(RecoveredMatrix[,5:L]); 
 D = as.numeric(DiedMatrix[,5:L]);
 } else { if (sum(is.na(ConfirmedMatrix[,1]))==0){
               C = as.numeric(colSums(ConfirmedMatrix[,5:L]));
               O = as.numeric(colSums(RecoveredMatrix[,5:L]));
               D = as.numeric(colSums(DiedMatrix[,5:L]));
          } else {k = which(is.na(ConfirmedMatrix[,1]));
               C = as.numeric(ConfirmedMatrix[min(k),5:L])
               O = as.numeric(RecoveredMatrix[min(k),5:L])
               D = as.numeric(DiedMatrix[min(k),5:L])}}
T = L-4;
Days = 1:T;

# Cumulative counts cannot decrease! #
for (t in 1:T){
       C[t] = max(C[1:t]);
       O[t] = max(O[1:t]);
       D[t] = max(D[1:t]);   }

c = c(0, C[2:(T)] - C[1:(T-1)]);   # Start from 0, let it have the same length T
d = c(0, D[2:(T)] - D[1:(T-1)]);
o = c(0, O[2:(T)] - O[1:(T-1)]);

# Skip low-count days in the beginning and in the end!!!

LowCountEarly = 1*(C < MinCases);
if (C[T] < MinCases){ # print(paste(Country,"reported only",C[T],"cases"));
   LowCountEarly = 1*(C < SmallCases); }
T0 = max(1,sum(LowCountEarly));
T1 = T

# MAh-day moving-average smoothing  (MAh can be changed in the beginning)
MAd = (MAh-1)/2;
U = rep(0,T);  Obar = U;   Dbar = U;   Cbar = U;
u = U;         obar = U;   dbar = U;   cbar = U;

for (t in 1:T){       # Moving averages
  M1 = max(1,t-MAd); M2 = min(t+MAd,T);
  Cbar[t] = round(mean(C[M1:M2]));
  Dbar[t] = round(mean(D[M1:M2]));
  Obar[t] = round(mean(O[M1:M2]));
  obar[t] = round(mean(o[M1:M2]));
  dbar[t] = round(mean(d[M1:M2]));
  cbar[t] = round(mean(c[M1:M2]));
                 }
C=Cbar; D=Dbar; O=Obar; c=cbar; o=obar; d=dbar; CC = C - O - D;

for (t in (h+1):T){  
  Cextra  = max( 0, min(  C[t:T] - O[t:T] - D[t:T] )) - U[t-1] ; 
  CCextra = max( 0, min( CC[t:T] - o[t:T] - d[t:T] )) - U[t-1] ; 
  Ucandidate = max(0, C[t-h] - O[t] - D[t]- U[t-1]);
  u[t]    = min( Ucandidate, Cextra, CCextra ); 

# u = new daily unobserved confirmed recovered
# u = number of people who were confirmed h days ago and haven't been reported as recovered or perished

  U[t] = U[t-1] + u[t];  # U = cumulative unobserved confirmed recovered
} 

R = O + U;   r = o + u;  # All confirmed recovered (R=cumulative, r=new daily counts)
CC = C - R - D;

###############################
### Estimate SIR parameters ###
###############################

#counts(Country);
CountryIndex = which(Population$Location == Country);
Npop = Population$PopTotal[CountryIndex]*1000; 
Continent = as.factor(Population$Continent[CountryIndex]);

# The new steps start here. All before this was as in SIRDmcmc.

# We no longer assume that all I get tested and that I = C.
# However, we still assume that the test is 100% accurate, that
# the recovery rates and the death rates are the same among tested
# and untested infected people, and that all confirmed deaths and recoveries
# are reported.

# SICROUD model with parameters:
# gamma = daily recovery rate
# delta = daily mortality rate
# theta = daily infection rate (per 1 infected contact)
# ttheta (theta-tilde) = overall daily infection rate = probability for any susceptible to become infected on a given day
# zeta  = testing rate, P(an infected person gets tested)
# eta   = recovery reporting probability, P(reported | recovered confirmed) = new parameter!
# rho   = new parameter = the average number of people in state I

# Empirical prior Beta(a,b) on each parameter, estim. by MOM
# Gamma distribution of rho, the Poisson parameter 

S = Npop-C;                    # Just for now, to guess the prior
C[C==0] = 1;  S[S==0] = 1;     # To avoid div/0
D[D>=Npop]=Npop-1;             # Same... Hopefully, pop. does not extinct

eps = 0.0000000001;

CC1 = CC + (CC==0);            # To avoid CC=0 in the denominator
r1  = r + (r==0);              # To avoid r=0 in the denominator
Delta = d[(T0+1):T1]/CC1[(T0+1):T1]; # Observed daily mortality rate
Ttheta = c[(T0+1):T1] / S[(T0+1):T1];   
Gamma = r[(T0+1):T1]/CC1[(T0+1):T1]; # Observed daily recovery rate
Eta   = o[(T0+1):T1]/r1[(T0+1):T1];  # Reporting rate = detectability of recoveries
Rho   = c[(T0+1):T1];                # rho = E(I)

Delta  = Delta + eps*(Delta==0);
Ttheta = Ttheta+ eps*(Ttheta==0);
Gamma  = Gamma + eps*(Gamma==0);
Eta    = Eta   + eps*(Eta  ==0);
Rho    = Rho   + eps*(Rho  ==0);

# Here we assume that report occurs on the same day, if a recovery is reported  
# Theta = Observed daily infection rate, per contact; 
# Here we "assume" infected i = confirmed c; then CC = currently @ Confirmed state
# Assuming i = c means that 1/2 of infected gets tested, 1/2 not tested

# For parameter estimation, we use the interval [T0, T1] 

Mgamma = max(eps, mean(Gamma)); Vgamma = max(eps^2, var(Gamma)); 
Agamma = max(eps, Mgamma*(Mgamma*(1-Mgamma)/Vgamma-1)); 
Bgamma = max(eps, Agamma*(1-Mgamma)/Mgamma);            

Mdelta = max(eps, mean(Delta)); Vdelta = max(eps^2, var(Delta)); 
Adelta = max(eps, Mdelta*(Mdelta*(1-Mdelta)/Vdelta-1));
Bdelta = max(eps, Adelta*(1-Mdelta)/Mdelta);

Mttheta = max(eps, mean(Ttheta)); Vttheta = max(eps^2, var(Ttheta)); 
Attheta = max(eps, Mttheta*(Mttheta*(1-Mttheta)/Vttheta-1));
Bttheta = max(eps, Attheta*(1-Mttheta)/Mttheta);

Meta = max(eps, mean(Eta)); Veta = max(eps^2, var(Eta)); 
Aeta = max(eps, Meta*(Meta*(1-Meta)/Veta-1));
Beta = max(eps, Aeta*(1-Meta)/Meta);

Mzeta = 0.1; Vzeta = 0.01;   # Pure guess. About 10% are tested in USA. 
# This is only a prior. With large Npop, it does not affect the result
Azeta = max(eps,Mzeta*(Mzeta*(1-Mzeta)/Vzeta-1));
Bzeta = max(eps,Azeta*(1-Mzeta)/Mzeta);

Mrho = max(eps, mean(Rho)); Vrho = max(eps, var(Rho));    # This prior "assumes" i=c
Arho = Mrho^2/Vrho; Brho = Mrho/Vrho; # B = frequency parameter of the Gamma dist of rho 

### M.C.M.C. #########
### Initialization ###

 TT = T1-T0+1; 
 I = rep(0,TT); R2 = I; D2 = I; i = I; r2 = I; d2 = I;  
 R0 = I;  # Basic reproduction number
 tthetagen = I; zetagen = I; etagen = I; 
 gammagen = I; deltagen = I; rhogen = I; 
 
 ### We'll also generate these counts, and then compare with observed. ###
 cgen = I; ogen = I; dgen = I;

S = Npop-C-R-D; I[1] = C[T0]; i[1] = I[1];   
tthetagen[1] = Mttheta; zetagen[1] = Mzeta; etagen[1] = Meta;
gammagen[1] = Mgamma; deltagen[1] = Mdelta; rhogen[1] = max(eps,CC[T0]);
cgen[1] = c[T0]; ogen[1] = o[T0]; dgen[1] = d[T0]; 

############################################################################
### MCMC starts here #######################################################
############################################################################

for (t in 2:(T1-T0+1)){ 
   tt = t+T0-1;      # t = day since T0, tt = day since Jan 22

   rho   = min(Npop,rgamma(1,Arho,Brho)); rhogen[t]  = rho;
   ttheta = rbeta(1,Attheta,Bttheta); tthetagen[t] = ttheta;  
   zeta  = rbeta(1,Azeta,Bzeta);      zetagen[t]  = zeta;
   gamma = rbeta(1,Agamma,Bgamma);    gammagen[t] = gamma;
   eta   = rbeta(1,Aeta,Beta);        etagen[t]   = eta;
   delta = rbeta(1,Adelta,Bdelta);    deltagen[t] = delta;

   R0[t] = ttheta / (gamma+delta);   

   I[t]  = c[tt] + rpois(1,rho);     # Based on the derived posterior    
   r2[t] = rbinom(1,I[t-1],gamma); R2[t] = R2[t-1] + r2[t];  # Does not need to be multiplied by (1-zeta)
   d2[t] = rbinom(1,I[t-1],delta); D2[t] = D2[t-1] + d2[t];
   S[tt]  = Npop - C[tt] - R2[t] - D2[t] - I[t];
   i[t]  = S[tt-1] - S[tt];

# Generate counts for cross-validation   
   cgen[t] = rbinom(1,I[t-1],zeta);
   rgen    = rbinom(1,CC[tt-1],gamma);
   ogen[t] = rbinom(1,rgen,eta);
   dgen[t] = rbinom(1,CC[tt-1],delta);

   # Now, the hyperparameter update, prior -> posterior
   Agamma  = Agamma + r[tt]; Bgamma = max(Bgamma + CC[tt-1] - r[tt],eps);
   Adelta  = Adelta + d[tt]; Bdelta = max(Bdelta + CC[tt-1] - d[tt],eps);
   Aeta    = Aeta   + o[tt]; Beta   = max(Beta   + r[tt]    - o[tt],eps);
   Azeta   = Azeta  + c[tt]; Bzeta  = max(Bzeta  + I[t-1]   - c[tt],eps);
   Attheta = Attheta+ i[t];  Bttheta= max(Bttheta+ S[tt-1]  - i[t], eps);
   Arho    = Arho   + I[t];  Brho   = Brho  + 1;

  #       print(paste("t=",t,"tt=",tt))
  #       print("parameters ttheta,rho,zeta,gamma,delta,R0")
  #       print(c(ttheta,rho,zeta,gamma,delta,R0[t]))
  #       print("Generated counts S[t-1],I[t-1],i[t],r2[t],S[t],I[t]")
  #       print(c(S[tt-1],I[t-1],i[t],r2[t],S[tt],I[t]))
  #       print("Hyperparameters Attheta,Bttheta,Alambda,Blambda,Azeta,Bzeta,Agamma,Bgamma,Adelta,Bdelta,Aeta,Beta")
  #       print(c(Attheta,Bttheta,Alambda,Blambda,Azeta,Bzeta,Agamma,Bgamma,Adelta,Bdelta,Aeta,Beta))
    
# We estimate ttheta from i[t] ~ Binomial(S,ttheta)
# For gamma and delta posterior updates, we used observed counts r, d
# and omitted unobserved but simulated r2, d2. 

}   # End of for loop over t (tt)

 ### Prepare results ###

 # Parameters

  theta.coef = (I+CC[T0:T1]) / (Npop - D[T0:T1] - D2);
  theta.coef = theta.coef + eps*(theta.coef == 0);
  thetagen = tthetagen / theta.coef;
#  thetagen = thetagen*(thetagen <= 1) + 1*(thetagen > 1);
 
 Ttheta = tthetagen[TT]
 Theta  = thetagen[TT]
 Rho    = rhogen[TT]
 Zeta   = zetagen[TT]
 Gamma  = gammagen[TT]
 Delta  = deltagen[TT]
 Eta    = etagen[TT]
 R0     = Theta/(Gamma+Delta)
 Duration = 1/(Gamma+Delta)  # Average duration of a disease

 # Counts #
 Iall = I[TT] + C[T] + R2[TT] + D2[TT];
 C = C[T];
 Rall = R[T] + R2[TT];
 R1 = R[T];
 O = O[T];
 U = U[T];
 R2 = R2[TT];
 Dall = D[T] + D2[TT];
 D1 = D[T]
 D2 = D2[TT]

 # Percentages #

Infection.rate              = Iall/Npop*100   
Testing.rate.infected       = C/Iall*100
Recovery.rate               = Rall/Iall*100
Mortality.rate              = Dall/Iall*100
Recovery.reporting.rate     = O/Rall*100

 Parameters = data.frame(Ttheta,Theta,Rho,Zeta,Gamma,Delta,Eta,R0,Duration)
 Counts = data.frame( Npop, Iall, C, Rall, R1, O, U, R2, Dall, D1, D2 );
 Percentages = data.frame(Infection.rate, Testing.rate.infected,
           Recovery.rate, Mortality.rate, Recovery.reporting.rate);

 return(data.frame(Country,Continent,Parameters,Counts,Percentages))
}   

SICROUD("US")
