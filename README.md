# Modeling Hunter Gatherer/Farmer Interactions


## About the simulation

This simulation models the interaction between Hunter Gatherers and Farmers in Neolthic Europe.
It represents the expansion of farming throughout Europe during this time period.
It is a Non-WF model featuring spatial competition, reproduction and learning. It is a 2D spatial model that can be displayed on a map.

### About the map:

The map is a square approximation of Europe which by default is 4,000 km^2.
This means that 1km on the screen (which is a 1x1) box is 0.00025.
The simulation can be run either on the map of Europe or on a simple black box.

### Color schemes:

#### Phenotypic coloring
Hunter gatherers (HGs) and farmers can be represented as colors by phenotype: The general defaults are listed below.
Another color can also be turned on to represent first generation offspring of a farmer and a HG. By default this is off.
```
Red: Hunter Gatherer
Blue: Farmer
Green: Hunter Gatherer who became a farmer by learning (Phenotypically now a farmer)
```

#### Genotypic coloring
The genotypic coloring scheme shows the proportion of farmer ancestry in the population rather than the phenotype of HG vs farmer.
The coloring is a spectrum from black to blue. The gradation is as follows:
```
Blue: Farmer
 |
Black: HG
```

## Code:

#### Initialize function

```
initialize()
{
	initializeSLiMModelType("nonWF");
	initializeSLiMOptions(dimensionality="xy", periodicity="");
```

This bit of code begins the model and creates the non-WF model and the xy dimensionality of the map.

##### Parameters

```
  // ---------------------------------------------------
	//  PARAMETERS --> Initialize Constant Params
	// ---------------------------------------------------
```

*Next begins the set up where parameters for the model will be set. Some of these will not need to be touched but some should be tuned to your liking.*

###### Parameters for how the population grows, interacts and moves

```
  // Carrying Capacities and Pop Sizes:
	// ***********************************
	defineConstant("SN", 1000); // Starting number of individuals
	
	defineConstant("HGK", 400); // carrying capacity for HGs for density dependent scaling
	defineConstant("FK", 2000); // carrying capacity for farmers for density dependent scaling
	defineConstant("K", 1300); // general k for age fitness scaling
 ```
These parameters discribe the starting number of individuals and the carrying capacity for HGs and farmers as well as a number used for the fitness scaling by age.
These parameters should be adjusted to your liking.

```
  // Movement and interaction Distances (ENTER IN SQUARE KILOMETERS):
	// ***********************************
	defineConstant("S_km", 200); // spatial competition distance
	defineConstant("MD_km", 100); // Mating distance
	defineConstant("LD_km", 100); // Learning distance
	defineConstant("OMD_km", 100); // Offspring movement away from parents. Right now OMD is the same for both populations but we could make it unique
	defineConstant("FMD_km", 300); // How far farmers diffuse away from their location
	defineConstant("HGMD_km", 400); // How far HGs diffuse away from their location
```

These parameters are all distance related and will be entered in km. This is for simplicity for the map.
The actual map will use the fraction of the 1x1 square that is equivalent to 1km. 

```
	// Learning, death and mating rate params:
	// ***********************************
	defineConstant("L", 0.005); // Learning rate 
	defineConstant("LP", 0.6); // Learning pecentage = the ratio of farmers to HGs required in an area for an individual HG to learn from a farmer 
	defineConstant("M", 0.1); // Mating rate
	
	// Age related mortality table
	defineConstant("age_scale", c(0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0));
 ```
 
 These parameters are all about rates. This is the likelihood that farmers will teach a HG how to farm. The first param L is just simply how often this interaction happens.
 The second (LP) is the ratio of farmers to HGs in an area that is required for a HG to learn. 
 For example, if a given area is majority HG it may be less likely that the minority farmers teach. 
 This of course can be changed as you see fit. By default here we see the param is 60%.
 
 Lastly, the age related mortality table is a life table. It is implimented similarly to SLiM recipe 16.2 Age structure (a life table model) by Benjamin C. Haller and Philipp W. Messer. It is decribed in the [manual](http://benhaller.com/slim/SLiM_Manual.pdf) in a succinct and useful way:
 
*"the addition of the defined constant [age_scale] ... is our life table. 
It gives the probability of mortality for each age; newly generated juveniles have a mortality of 0.7 (i.e., 70%), then the mortality drops to zero for several years, and then it ramps gradually upward with increasing age until it reaches 1.0 [at which all individuals of this age] will die. Note that this is only the age-related mortality; density-dependence will also cause mortality, as we will see below, but that will be additional to this age-related mortality, which would occur even in
a population that was not limited by its density." -Benjamin C. Haller and Philipp W. Messer (SLiM Manual)*


The ages of individuals correspond to the indices of the life table.
 
 
###### Parameters for user preferences on how the model will look and run

```
	// ---------------------------------------------------
	// RUN TIME PREFERENCE PARAMETERS
	// ---------------------------------------------------
```

First, we will look at coloring schemes. These are discribed at the top of the markdown.

*Color schemes*
```
	// Determine Coloring Schemes:
	// ***********************************	
	defineConstant("Color_scheme", 1); // Parameter of 0 = Genomic Coloring, 1 = Behavoiral Coloring
	defineConstant("Color_option1", 0); // Parameter of 0 = No special color for 1st generation 'hybrid' offspring, 1 = Coloring for 1st generation 'hybrid' offspring
	defineConstant("Color_option2", 1); // Parameter of 0 = No special color for 1st generation HG individuals who have learned, 1 = Coloring for 1st generation HG individuals who have learned
```
Color_scheme is a binary choice for if you want genomic coloring gradient or phenotypic behavorial coloring

The next two are options that determine if you want special colors for specific conditons like those who have recently learned or first generation "hybrids"

*Preferences regarding the map the individuals move around on*

```
	
	// Map Prefs:
	// ***********************************	
	defineConstant("map", 1); // 0 = black square, 1 = map of europe
	defineConstant("mapfile", "C:/PATH/GOES/HERE/TO/MAP/europe.png"); // Path to URL
	defineConstant("map_size_km2", 4000);
```

The first parameter determines if you would like a black square for simplicity or a map of europe

The next is a path in YOUR file system that the map file can be found in. The map is avaliable for download in the repo. This path points to where the map is in your system so it MUST be changed to fit your file stucture.
