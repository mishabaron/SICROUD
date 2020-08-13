# Data from https://data.humdata.org/dataset/novel-coronavirus-2019-ncov-cases

ConfirmedGlobal = read.csv(url("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_confirmed_global.csv&filename=time_series_covid19_confirmed_global.csv"));
DiedGlobal = read.csv(url("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_deaths_global.csv&filename=time_series_covid19_deaths_global.csv"));
RecoveredGlobal = read.csv(url("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_recovered_global.csv&filename=time_series_covid19_recovered_global.csv"));

States   = read.csv(url("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2Fnytimes%2Fcovid-19-data%2Fmaster%2Fus-counties.csv&filename=us-states.csv"))


# Data from https://data.humdata.org/dataset/novel-coronavirus-2019-ncov-cases
# Defines functions: conf, mort, counts, counts.state, counts.county


setwd("C:\\Users\\baron\\Documents\\Research\\Grants\\NSF-RAPID 2020\\Prep to research")
Population = read.csv("Population.csv");



########################### CONFIRMED CASES #################################

conf = function(Country){
CountryIndexC = which(ConfirmedGlobal$Country.Region == Country);
L = dim(ConfirmedGlobal)[2];   # Number of days. Data start from column 5 (01-22-2020)
ConfirmedMatrix = ConfirmedGlobal[CountryIndexC,]

if (nrow(ConfirmedMatrix)==1){ Conf <- ConfirmedMatrix[,5:L]
  } else { if (sum(is.na(ConfirmedMatrix[,1]))==0){
      Conf <- colSums(ConfirmedMatrix[,5:L])
  } else {k = which(is.na(ConfirmedMatrix[,1]));
      Conf <- ConfirmedMatrix[min(k),5:L]}}

T <<- L-4;  Days <<- 1:T;  Dates <<- as.Date(Days,origin="2020-01-21")

C <- as.numeric(Conf);

# Cumulative counts cannot decrease! #
for (t in 1:T){   C[t] = max(C[1:t]);  }

c <- c(0,C[2:(T)] - C[1:(T-1)]);

LastCHP <<- T-30;         # min(T-75+which.max(c[(T-75):T]), T-20); 

ConfDailyLastCHP = c[LastCHP:T];

DaysLastCHP <<- Days[LastCHP:T];
reg <<- lm( ConfDailyLastCHP ~ DaysLastCHP );

# Introduce days of the week
Weekday = rep("Thursday",T)     # Day 1 = Jan 22 = Thursday
Weekday[seq(2,T,7)] = "Friday"
Weekday[seq(3,T,7)] = "Saturday"
Weekday[seq(4,T,7)] = "Sunday"
Weekday[seq(5,T,7)] = "Monday"
Weekday[seq(6,T,7)] = "Tuesday"
Weekday[seq(7,T,7)] = "Wednesday"
WeekdayLastCHP = Weekday[LastCHP:T]

# Regression since day #LastCHP with weekly periodicity
reg2 = lm( ConfDailyLastCHP ~ DaysLastCHP + WeekdayLastCHP);
Cslope <<- coef(reg2)[2]

chat = predict(reg2)
chat <<- chat*(chat > 0)

# Predict the next h days...
h = 10;
DaysNew <<- seq(T+1,T+h,1);
DatesNew <<- as.Date(DaysNew,origin="2020-01-21")

WeekdayNew = rep("?",h); for (k in 1:7){WeekdayNew[seq(k,h,7)] = Weekday[T-7+k]}

cnew = predict(reg2, data.frame(DaysLastCHP=DaysNew, WeekdayLastCHP=WeekdayNew)); 
cnew <<- cnew*(cnew > 0)

 par(mfrow=c(1,1))

#Dates = as.Date(Days,origin="2020-01-21")
plot(Days,c,col="blue",type="b",lwd=3,main=Country)
#plot(Days,c,col="blue",type="b",lwd=3,main=Country,xaxt="n")
#axis.Date(1, at=seq(min(Dates), min(Dates)+T, 3))

abline(reg,col="red",lwd=1)
points(DaysLastCHP, chat, type="b", col="orange", lwd=3 )
points(DaysNew,cnew, type="b", col="orange", lwd=3 )
legend(0, max(ConfDailyLastCHP)/2, legend=c("Actual counts", "Predicted counts"),
        col=c("blue", "orange"), cex=1, lwd=3, bg='white') 

print(paste("Confirmed cases, last 10 days (slope =",round(Cslope,1),"recently):"))
print(c[(T-9):T])

Ctotal <<- C[T];

# print("Next 10 days forecast:")
# print(round(YhatNew))
return(c)
}



##################### CASUALTIES #########################################

mort = function(Country){
CountryIndexD = which(DiedGlobal$Country.Region == Country);
L = dim(DiedGlobal)[2];   # Number of days. Data start from column 5 (01-22-2020)
DiedMatrix = DiedGlobal[CountryIndexD,]

if (nrow(DiedMatrix)==1){ Died <- DiedMatrix[,5:L] 
  } else { if (sum(is.na(DiedMatrix[,1]))==0){
      Died <- colSums(DiedMatrix[,5:L])
  } else {k = which(is.na(DiedMatrix[,1]));
      Died <- DiedMatrix[min(k),5:L]}}

T = L-4; Days = 1:T;
D <- as.numeric(Died)

# Cumulative counts cannot decrease! #
for (t in 1:T){  D[t] = max(D[1:t]);   }

d <- c(0, D[2:(T)] - D[1:(T-1)]);
 
DLastCHP <- min(which.max(d),T-20); 
DiedDailyLastCHP = d[DLastCHP:T];

DaysDLastCHP = Days[DLastCHP:T];
reg = lm( DiedDailyLastCHP ~ DaysDLastCHP );

# Introduce days of the week
Weekday = rep("Thursday",T)     # Day 1 = Jan 22 = Thursday
Weekday[seq(2,T,7)] = "Friday"
Weekday[seq(3,T,7)] = "Saturday"
Weekday[seq(4,T,7)] = "Sunday"
Weekday[seq(5,T,7)] = "Monday"
Weekday[seq(6,T,7)] = "Tuesday"
Weekday[seq(7,T,7)] = "Wednesday"
WeekdayDLastCHP = Weekday[DLastCHP:T]

# Regression since day #LastCHP with weekly periodicity
reg2 = lm( DiedDailyLastCHP ~ DaysDLastCHP + WeekdayDLastCHP);
Dslope <<- coef(reg2)[2];
dhat = predict(reg2)
dhat = dhat*(dhat > 0)

# Predict the next h days...
h = 10;
DaysNew = seq(T+1,T+h,1);
WeekdayNew = rep("?",h); for (k in 1:7){WeekdayNew[seq(k,h,7)] = Weekday[T-7+k]}

dnew = predict(reg2, data.frame(DaysDLastCHP=DaysNew, WeekdayDLastCHP=WeekdayNew)); 
dnew <<- dnew *(dnew > 0)

 par(mfrow=c(1,1))
plot(Days,d,col="purple",type="b",lwd=3, main=Country)
abline(reg,col="red",lwd=1)
points(DaysDLastCHP, dhat, type="b", col="orange", lwd=3 )
points(DaysNew,dnew , type="b", col="orange", lwd=3 )
legend(0, max(d)/2, legend=c("Actual counts", "Predicted counts"),col=c("blue", "orange"), cex=1, lwd=3, bg='white') 

print(paste("Casualties, last 10 days (slope =",round(Dslope,1),"recently):"))
print(d[(T-9):T])

# print("Next 10 days forecast:")
# print(round(dnew ))

Dtotal <<- D[T];
return(d)
}


###################### RECOVERED ##################################

rec = function(Country){
CountryIndexR = which(RecoveredGlobal$Country.Region == Country);
L = dim(RecoveredGlobal)[2];   # Number of days. Data start from column 5 (01-22-2020)
RecoveredMatrix = RecoveredGlobal[CountryIndexR,]

if (nrow(RecoveredMatrix)==1){ Rec <- RecoveredMatrix[,5:L]
  } else { if (sum(is.na(RecoveredMatrix[,1]))==0){
      Rec <- colSums(RecoveredMatrix[,5:L])
  } else {k = which(is.na(RecoveredMatrix[,1]));
      Rec <- RecoveredMatrix[min(k),5:L]}}

T <<- L-4;  Days <<- 1:T;
R <- as.numeric(Rec);

# Cumulative counts cannot decrease! #
for (t in 1:T){R[t] = max(R[1:t]);}

r <- c(0,R[2:(T)] - R[1:(T-1)]);

plot(Days,r,col="green",type="b",lwd=3,main=Country)
legend(0, max(r)/2, legend=c("Reported recovered counts"),
        col=c("green"), cex=1, lwd=3, bg='white') 

Rtotal <<- R[T];

return(r)
}


###################### COMBINED COUNTS ###############################

counts = function(Country){
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

LastCHP <- T - 30;        # min(T-75+which.max(c[(T-75):T]),T-20); 
cLastCHP = c[LastCHP:T];
dLastCHP = d[LastCHP:T];

Days = 1:T; DaysLastCHP <- Days[LastCHP:T];
reg = lm( cLastCHP ~ DaysLastCHP );

# Introduce days of the week
Weekday = rep("Tuesday",T)     # Day 1 = Jan 21 = Tuesday
Weekday[seq(2,T,7)] = "Wednesday"
Weekday[seq(3,T,7)] = "Thursday"
Weekday[seq(4,T,7)] = "Friday"
Weekday[seq(5,T,7)] = "Saturday"
Weekday[seq(6,T,7)] = "Sunday"
Weekday[seq(7,T,7)] = "Monday"
WeekdayLastCHP = Weekday[LastCHP:T]

# Regression since day #LastCHP with weekly periodicity
reg2c = lm( cLastCHP ~ DaysLastCHP + WeekdayLastCHP);
Cslope <- coef(reg2c)[2]
chat = predict(reg2c)
chat <- chat*(chat > 0)

reg2d = lm( dLastCHP ~ DaysLastCHP + WeekdayLastCHP);
Dslope <- coef(reg2d)[2]
dhat = predict(reg2d)
dhat <- dhat*(dhat > 0)

# Predict the next h days...
h = 10;
DaysNew <- seq(T+1,T+h,1);
WeekdayNew = rep("?",h); for (k in 1:7){WeekdayNew[seq(k,h,7)] = Weekday[T-7+k]}

chatNew = predict(reg2c, data.frame(DaysLastCHP=DaysNew, WeekdayLastCHP=WeekdayNew)); 
chatNew <- chatNew *(chatNew > 0)
dhatNew = predict(reg2d, data.frame(DaysLastCHP=DaysNew, WeekdayLastCHP=WeekdayNew)); 
dhatNew <- dhatNew *(dhatNew > 0)

Ctotal <- C[T]; Dtotal <- D[T]; c.yesterday <- c[T]; d.yesterday <- d[T];
MortRate <- Dtotal/Ctotal;
Cpercent <- C[T]/sum(States$cases)

plot(Days,c,col="blue",type="b",lwd=3, ylab="Counts", xlim=c(0,T+10),
  main=paste(Country,"- ",Ctotal,"confirmed cases,",c[T],"yesterday, slope =",round(Cslope)))
abline(reg,col="red")
points( DaysLastCHP, chat,    type="b", col="orange", lwd=2 )
points( DaysNew,   chatNew, type="b", col="orange", lwd=2 )
points( Days,      d,       col="red",type="b",lwd=2)
points( Days,      c,       col="blue", type="b",lwd=3)
points( DaysNew,   dhatNew, type="b", col="orange", lwd=2 )

legend(0, max(c), 
  legend=c("Actual confirmed cases", "Actual fatalities", "Forecast"),
  col=c("blue", "red", "orange"), cex=1, lwd=3, bg='white') 


print(paste("Confirmed cases, last 10 days (slope =",round(Cslope,1),"recently):"))
print(c[(T-9):T])

print(paste("Casualties, last 10 days (slope =",round(Dslope,1),"recently):"))
print(d[(T-9):T])

MortRate = Dtotal/Ctotal;

print(paste("Official data:",Ctotal,"confirmed cases;", Dtotal,"casualties; mortality rate =",round(100*MortRate,1),"%"))
}

### Quiet Counts, printing suppressed, only plot, and no regression ###########

counts.quiet = function(Country){
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

Ctotal <- C[T]; Dtotal <- D[T]; c.yesterday <- c[T]; d.yesterday <- d[T];
MortRate <- Dtotal/Ctotal;
Cpercent <- C[T]/sum(States$cases)

par(mfrow=c(1,1))
plot(Days,c,col="blue",type="b",lwd=3, ylab="Counts", xlim=c(0,T),
  main=paste(Country,"- ",Ctotal,"confirmed cases,",c[T],"yesterday"))
points( Days,      d,       col="red",type="b",lwd=2)
points( Days,      c,       col="blue", type="b",lwd=3)

legend(0, max(c), 
  legend=c("Confirmed cases", "Casualties"),
  col=c("blue", "red"), cex=1, lwd=3, bg='white') 
}

##################### STATES ##############################################

# From https://data.humdata.org/dataset/nyt-covid-19-data?
# FIPS = Federal Information Processing Standards, https://www.census.gov/quickfacts/fact/note/US/fips

# States   = read.csv(url("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2Fnytimes%2Fcovid-19-data%2Fmaster%2Fus-counties.csv&filename=us-states.csv"))

Date = as.numeric(States$date)  # In date format, such as 2020-06-13
Day1 = min(Date)
T = as.numeric(max(Date)-Day1+1) # Number of days so far. Day 1 = Jan 21, 2020

setwd("C:\\Users\\baron\\Documents\\Research\\Grants\\NSF-RAPID 2020\\Prep to research");
StatePop = read.csv("PopulationStates.csv")


counts.state = function(State){
C = rep(0,T); D = rep(0,T); # Confirmed and Died cumulative
c = rep(0,T); d = rep(0,T); # Confirmed and Died daily

for (t in 1:T){
  Indices = which(States$state == State & Date == Day1 + t - 1);
  C[t] = sum(States$cases[Indices]);
  D[t] = sum(States$deaths[Indices]);
}
  c[2:T] = C[2:T]-C[1:(T-1)];
  d[2:T] = D[2:T]-D[1:(T-1)];
 
LastCHP <- T - 30;        # min(T-75+which.max(c[(T-75):T]),T-20); 
cLastCHP = c[LastCHP:T];
dLastCHP = d[LastCHP:T];

Days = 1:T; DaysLastCHP <- Days[LastCHP:T];
reg = lm( cLastCHP ~ DaysLastCHP );

# Introduce days of the week
Weekday = rep("Tuesday",T)     # Day 1 = Jan 21 = Tuesday
Weekday[seq(2,T,7)] = "Wednesday"
Weekday[seq(3,T,7)] = "Thursday"
Weekday[seq(4,T,7)] = "Friday"
Weekday[seq(5,T,7)] = "Saturday"
Weekday[seq(6,T,7)] = "Sunday"
Weekday[seq(7,T,7)] = "Monday"
WeekdayLastCHP = Weekday[LastCHP:T]

# Regression since day #LastCHP with weekly periodicity
reg2c = lm( cLastCHP ~ DaysLastCHP + WeekdayLastCHP);
Cslope <- coef(reg2c)[2]
chat = predict(reg2c)
chat <- chat*(chat > 0)

reg2d = lm( dLastCHP ~ DaysLastCHP + WeekdayLastCHP);
Dslope <- coef(reg2d)[2]
dhat = predict(reg2d)
dhat <- dhat*(dhat > 0)

# Predict the next h days...
h = 10;
DaysNew <- seq(T+1,T+h,1);
WeekdayNew = rep("?",h); for (k in 1:7){WeekdayNew[seq(k,h,7)] = Weekday[T-7+k]}

chatNew = predict(reg2c, data.frame(DaysLastCHP=DaysNew, WeekdayLastCHP=WeekdayNew)); 
chatNew <- chatNew *(chatNew > 0)
dhatNew = predict(reg2d, data.frame(DaysLastCHP=DaysNew, WeekdayLastCHP=WeekdayNew)); 
dhatNew <- dhatNew *(dhatNew > 0)

Ctotal <- C[T]; Dtotal <- D[T]; c.yesterday <- c[T]; d.yesterday <- d[T];
MortRate <- Dtotal/Ctotal;
Cpercent <- C[T]/sum(States$cases)

# Index = which(StatePop$State==State);

par(mfrow=c(1,1))
plot(Days,c,col="blue",type="b",lwd=3, ylab="Counts", xlim=c(0,T+10),
  main=paste(State," - ",Ctotal,"confirmed cases,",c[T],"yesterday, slope =",round(Cslope)))
abline(reg,col="red")
points( DaysLastCHP, chat,    type="b", col="orange", lwd=2 )
points( DaysNew,   chatNew, type="b", col="orange", lwd=2 )
points( Days,      d,       col="red",type="b",lwd=2)
points( Days,      c,       col="blue", type="b",lwd=3)
points( DaysNew,   dhatNew, type="b", col="orange", lwd=2 )

legend(0, max(c), 
  legend=c("Actual confirmed cases", "Actual fatalities", "Forecast"),
  col=c("blue", "red", "orange"), cex=1, lwd=3, bg='white') 


print(paste("Confirmed cases, last 10 days (slope =",round(Cslope,1),"recently):"))
print(c[(T-9):T])

print(paste("Casualties, last 10 days (slope =",round(Dslope,1),"recently):"))
print(d[(T-9):T])

MortRate = Dtotal/Ctotal;

print(paste("Official data:",Ctotal,"confirmed cases;", Dtotal,"casualties; mortality rate =",round(100*MortRate,1),"%"))
}


##################### COUNTIES ##############################################

counts.county = function(State,County){
C = rep(0,T); D = rep(0,T); # Confirmed and Died cumulative
c = rep(0,T); d = rep(0,T); # Confirmed and Died daily

for (t in 1:T){
  Indices = which(States$state == State & States$county == County & Date == Day1 + t - 1);
  C[t] = sum(States$cases[Indices]);
  D[t] = sum(States$deaths[Indices]);
}
  c[2:T] = C[2:T]-C[1:(T-1)];
  d[2:T] = D[2:T]-D[1:(T-1)];
 
LastCHP <- T - 30;        # min(T-75+which.max(c[(T-75):T]),T-20); 
cLastCHP = c[LastCHP:T];
dLastCHP = d[LastCHP:T];

Days = 1:T; DaysLastCHP <- Days[LastCHP:T];
reg = lm( cLastCHP ~ DaysLastCHP );

# Introduce days of the week
Weekday = rep("Tuesday",T)     # Day 1 = Jan 21 = Tuesday
Weekday[seq(2,T,7)] = "Wednesday"
Weekday[seq(3,T,7)] = "Thursday"
Weekday[seq(4,T,7)] = "Friday"
Weekday[seq(5,T,7)] = "Saturday"
Weekday[seq(6,T,7)] = "Sunday"
Weekday[seq(7,T,7)] = "Monday"
WeekdayLastCHP = Weekday[LastCHP:T]

# Regression since day #LastCHP with weekly periodicity
reg2c = lm( cLastCHP ~ DaysLastCHP + WeekdayLastCHP);
Cslope <- coef(reg2c)[2]
chat = predict(reg2c)
chat <- chat*(chat > 0)

reg2d = lm( dLastCHP ~ DaysLastCHP + WeekdayLastCHP);
Dslope <- coef(reg2d)[2]
dhat = predict(reg2d)
dhat <- dhat*(dhat > 0)

# Predict the next h days...
h = 10;
DaysNew <- seq(T+1,T+h,1);
WeekdayNew = rep("?",h); for (k in 1:7){WeekdayNew[seq(k,h,7)] = Weekday[T-7+k]}

chatNew = predict(reg2c, data.frame(DaysLastCHP=DaysNew, WeekdayLastCHP=WeekdayNew)); 
chatNew <- chatNew *(chatNew > 0)
dhatNew = predict(reg2d, data.frame(DaysLastCHP=DaysNew, WeekdayLastCHP=WeekdayNew)); 
dhatNew <- dhatNew *(dhatNew > 0)

Ctotal <- C[T]; Dtotal <- D[T]; c.yesterday <- c[T]; d.yesterday <- d[T];
MortRate <- Dtotal/Ctotal;
Cpercent <- C[T]/sum(States$cases)

# Index = which(StatePop$State==State);

par(mfrow=c(1,1))
plot(Days,c,col="blue",type="b",lwd=3, ylab="Counts", xlim=c(0,T+10),
  main=paste(County," County - ",Ctotal,"confirmed cases,",c[T],"yesterday, slope =",round(Cslope)))
abline(reg,col="red")
points( DaysLastCHP, chat,    type="b", col="orange", lwd=2 )
points( DaysNew,   chatNew, type="b", col="orange", lwd=2 )
points( Days,      d,       col="red",type="b",lwd=2)
points( Days,      c,       col="blue", type="b",lwd=3)
points( DaysNew,   dhatNew, type="b", col="orange", lwd=2 )

legend(0, max(c), 
  legend=c("Actual confirmed cases", "Actual fatalities", "Forecast"),
  col=c("blue", "red", "orange"), cex=1, lwd=3, bg='white') 


print(paste("Confirmed cases, last 10 days (slope =",round(Cslope,1),"recently):"))
print(c[(T-9):T])

print(paste("Casualties, last 10 days (slope =",round(Dslope,1),"recently):"))
print(d[(T-9):T])

MortRate = Dtotal/Ctotal;

print(paste("Official data:",Ctotal,"confirmed cases;", Dtotal,"casualties; mortality rate =",round(100*MortRate,1),"%"))
}

counts.county("Maryland","Montgomery")
counts.state("Maryland")
counts("US")


