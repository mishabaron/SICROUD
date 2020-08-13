# SICROUD
Modeling and estimation for COVID-19 under-reported counts

The project includes two R codes - humdataload.R and SICROUD.R

### humdataload.R ###

The R program "humdataload.R" reads the most recently 
updated COVID-19 counts data from The Humanitarian Data Exchange (HDX)
site https://data.humdata.org. These data contain confirmed infected, 
recovered, and perished counts, by country and day. US dataset contains
counts for all US states and counties.

The data are augmented with population information in CSV files
"Population.csv" and "PopulationStates.csv".

Download these CSV files and insert the proper path.

The code "humdataload.R" contains several functions:

conf = function(Country)   # Obtains and plots confirmed cases as a time plot
mort = function(Country)   # Same, for perished cases
rec = function(Country)    # Same, for recovered cases

counts = function(Country) # Combination of conf, mort, and rec. Graphs all the
                           # counts in a combined time plot

counts.state = function(State)  		# Same, for any US state
counts.county = function(State,County)	# Same, for any county within a State

The inputs Country, State, and County are of character type. 
For example,

counts("US");
counts.state("Maryland");
counts.county("Maryland","Montgomery");


### SICROUD.R ###

SICROUD.R is a function of Country, SICROUD = function(Country), where 
Country is of character type. For example, SICROUD("US");

SICROUD.R fits the SICROUD model, estimates all its parameters and
unobserved counts, and prints the summary output.

It requires to run the code "humdataload.R" first.

The output consists of:
- the estimated infection rate
- the estimated testing rate
- the estimated recovery rate
- the estimated mortality rate
- the estimated reporting rate
- the estimated frequency of infected untested individuals
- the estimated basic reproduction number
- the estimated total number of disease transmissions
- the estimated total number of recoveries
- the estimated total number of casualties
- the estimated number of unconfirmed transmissions
- the estimated number of unreported recoveries

Examples:
> SICROUD("US")
> SICROUD("Sweden")
> x = SICROUD("Ireland")

