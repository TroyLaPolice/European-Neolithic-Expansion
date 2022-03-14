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

## Phenotypic coloring
Hunter gatherers (HGs) and farmers can be represented as colors by phenotype: The general defaults are listed below.
Another color can also be turned on to represent first generation offspring of a farmer and a HG.
```
Red: Hunter Gatherer
Blue: Farmer
Green: Hunter Gatherer who became a farmer by learning (Phenotypically now a farmer)
```

## Genotypic coloring
The genotypic coloring scheme shows the proportion of farmer ancestry in the population rather than the phenotype of HG vs farmer.
The coloring is a spectrum from black to blue. The gradation is as follows:
```
Blue: Farmer
 |
Black: HG
```

## Code:

## Initialize function

```
initialize()
{
	initializeSLiMModelType("nonWF");
	initializeSLiMOptions(dimensionality="xy", periodicity="");
```

This bit of code begins the model and creates the non-WF model and the xy dimensionality of the map.

### Parameters

```
  // ---------------------------------------------------
	//  PARAMETERS --> Initialize Constant Params
	// ---------------------------------------------------
```

Next begins the set up where parameters for the model will be set. Some of these will not need to be touched but some should be tuned to your liking.

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
  // Learning and mating rate params:
	// ***********************************
	defineConstant("L", 0.005); // Learning rate 
	defineConstant("LP", 0.6); // Learning pecentage = the ratio of farmers to HGs required in an area for an individual HG to learn from a farmer 
	defineConstant("M", 0.1); // Mating rate
 ```
 
 These parameters are all about rates. This is the likelihood that farmers will teach a HG how to farm. The first param L is just simply how often this interaction happens.
 The second (LP) is the ratio of farmers to HGs in an area that is required for a HG to learn. 
 For example, if a given area is majority HG it may be less likely that the minority farmers teach. 
 This of course can be changed as you see fit. By default here we see the param is 60%.
 
