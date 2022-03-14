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
	defineConstant("min_repro_age", 12); // Individuals MUST be OLDER than this age to reproduce
	
	// Age related mortality table
	defineConstant("age_scale", c(0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0));
 ```
 
 These parameters are all about rates for mating, learning and death. This is the likelihood that farmers will teach a HG how to farm. The first param L is just simply how often this interaction happens.
 The second (LP) is the ratio of farmers to HGs in an area that is required for a HG to learn. 
 For example, if a given area is majority HG it may be less likely that the minority farmers teach. 
 This of course can be changed as you see fit. By default here we see the param is 60%.
 
 Next is a minimum age required for reproduction. Inidividuals MUST BE OLDER than the provided age in this parameter for them to be able to reproduce. This prevents infants and small children from being able to reproduce which is unrealistic.
 
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

*First, we will look at coloring schemes. These are discribed at the top of the markdown.*

```
	// Determine Coloring Schemes:
	// ***********************************	
	defineConstant("Color_scheme", 1); // Parameter of 0 = Genomic Coloring, 1 = Behavoiral Coloring
	defineConstant("Color_option1", 0); // Parameter of 0 = No special color for 1st generation 'hybrid' offspring, 1 = Coloring for 1st generation 'hybrid' offspring
	defineConstant("Color_option2", 1); // Parameter of 0 = No special color for 1st generation HG individuals who have learned, 1 = Coloring for 1st generation HG individuals who have learned
```
Color_scheme is a binary choice for if you want genomic coloring gradient or phenotypic behavorial coloring

The next two are options that determine if you want special colors for specific conditons like those who have recently learned or first generation "hybrids"

*Preferences regarding the map on which the individuals move*

```
	
	// Map Prefs:
	// ***********************************	
	defineConstant("map", 1); // 0 = black square, 1 = map of europe
	defineConstant("mapfile", "C:/PATH/GOES/HERE/TO/MAP/europe.png"); // Path to URL
	defineConstant("map_size_km2", 4000);
```

The first parameter determines if you would like a black square for simplicity or a map of europe

The next is a path in YOUR file system that the map file can be found in. The map is avaliable for download in the repo. This path points to where the map is in your system so it MUST be changed to fit your file stucture.

*Next are several parameters that will not typically need to be changed but some can be if you want*

```
// ---------------------------------------------------
	//  PARAMETERS FOR MAP SCALING (NO NEED TO ADJUST)
	// ---------------------------------------------------
	
	defineConstant("S", (S_km / map_size_km2)); // spatial competition distance
	defineConstant("MD", (MD_km / map_size_km2)); // Mating distance
	defineConstant("LD", (LD_km / map_size_km2)); // Learning distance
	defineConstant("OMD", (OMD_km / map_size_km2)); // Offspring movement away from parents. Right now OMD is the same for both populations but we could make it unique
	defineConstant("FMD", (FMD_km / map_size_km2)); // How far farmers diffuse away from their location
	defineConstant("HGMD", (HGMD_km / map_size_km2)); // How far HGs diffuse away from their location
```

This block of code simply does the conversion from the values we input in km to make the model work in a 1x1 square.
This is handy becuase it avoids you having to do the math yourself. **This does not need to be adjusted**

*The next two blocks simply initialize the genetic component that is the marker for farmer ancesty and the interactions between individuals*

```
	
	// ----------------------------------------------------
	//  GENETIC COMPONENT --> Initialize Genomic Elements
	// ----------------------------------------------------
	initializeMutationType("m1", 0.5, "f", 0.0); // Tag farmer ancestry
	m1.convertToSubstitution = F;
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, 2999);
	initializeMutationRate(0.0);
	initializeRecombinationRate(0.01);
```

This introduces a marker mutation for farmers that we can see recombine with HGs as the simulation progresses. This is how genotypic coloring works.

```
	// ---------------------------------------------------
	//  INTERACTIONS --> Interaction Initialization
	// ---------------------------------------------------
	// spatial competition
	initializeInteractionType(1, "xy", reciprocal=T, maxDistance=S);
	
	// spatial mate choice
	initializeInteractionType(2, "xy", reciprocal=T, maxDistance=MD);
	
	// spatial mate choice
	initializeInteractionType(3, "xy", reciprocal=T, maxDistance=LD);
```

This initializes the interaction types between individuals:

*1) Spatial competion between nearby individuals*

*2) Mating*

*3) Learning*

These all take place within a certain distance range specified by parameters above.

#### First generation of the simulation- building inital population
```
1 early()
{
	if (map == 1)
	{
		// Set up map
		p1.setSpatialBounds(c(0.0, 0.0, 1.0, 1.0));
		mapImage = Image(mapfile);
		p1.defineSpatialMap("world", "xy", 1.0 - mapImage.floatK, valueRange=c(0.0, 1.0), colors=c("#ffffff", "#111111"));
		
		// start near a specific map location
		for (ind in p1.individuals)
		{
			do
				newPos = c(runif(1, 0, 1.0), runif(1, 0, 1.0));
			while (!p1.pointInBounds(newPos) | p1.spatialMapValue("world", newPos) == 0.0);
			ind.setSpatialPosition(newPos);
		}
	}
	
```

This section sets up how the image file of the map of Europe works if the user specifies they want to run a sim with the map.

It defines the map and its bounds.

```
	// Define z param in offspring (phenotype, 0 = HG, 1 = Farmer)
	// Make individuals near Turkey for farmers
	p1.individuals[p1.individuals.x > 0.7 & p1.individuals.y < 0.3].z = 1;
	
	// Tag genomic ancestry of farmers with marker mutations (m1)
	// Each marker mutation represents 1Mb
	indFarmers = p1.individuals[p1.individuals.z == 1];
	indFarmers.genomes.addNewMutation(m1, 0.0, 0:2999);
	
	// Add color to represent phenotype
	for (i in p1.individuals)
	{
		// -----------------------
		// Color Based on Phenotype
		// -----------------------
		value = i.countOfMutationsOfType(m1) / 6000;
		i.color = rgb2color(hsv2rgb(c(0.6, 1.0, value)));
		
		if (Color_scheme == 1)
		{
			// -----------------------
			// Color Based on Behavior
			// -----------------------
			// HGs are red, farmers are blue
			if (i.z == 0)
				i.color = "red";
			else
				i.color = "blue";
		}
	}
}
```

This portion sets up the z cordinate for indivuals and adds the marker mutation.

The z cordinate here represents the behavoral phenotype- farming vs HGing.

Individuals located near Turkey on the map begin as farmers and slowly expand throughout Europe.

This portion also sets up how the coloring schemes work. 

It sets it up both by phenotype (behavorial) coloring and genomic coloring based on user provided preference above.

#### Reproduction

```
first()
{
	// look for mates
	i2.evaluate();
}

reproduction()
{
	// ---------------------------------------------------
	//  MATING --> Individuals mate with those close by
	// ---------------------------------------------------
	// choose nearest neighbor as a mate, within the max distance
	if (individual.age > 12) // Reproductive age of the individual must be reached before mating
	{
		mate = i2.nearestNeighbors(individual, 1);
		if (mate.size())
		{
			if (mate.age > min_repro_age) // Reproductive age of the individual must be reached before mating
			{
				// Frequency of the interaction
				for (i in seqLen(rpois(1, M)))
				{
					// Only runs if a mate is nearby
					if (mate.size())
					{
						offspring = subpop.addCrossed(individual, mate);
						
						// -----------------------
						// Color Based on Phenotype
						// -----------------------
						value = offspring.countOfMutationsOfType(m1) / 6000;
						offspring.color = rgb2color(hsv2rgb(c(0.6, 1.0, value)));
						
						// Define z param in offspring (phenotype, 0 = HG, 1 = Farmer)
						offspring.z = 0;
						
						// If both parents are farmers the child is a farmer
						if (individual.z == 1 & mate.z == 1)
							offspring.z = 1;
						if (Color_scheme == 1)
						{
							// -----------------------
							// Color Based on Behavior
							// -----------------------
							// Add color to represent phenotype
							if (offspring.z == 0)
								offspring.color = "red";
							else
								offspring.color = "blue";
						}
						
						// If the one parent is a farmer and one parent is a HG,
						// 50% chance they become a farmer
						if (individual.z != mate.z)
						{
							offspring.z = rbinom(1, 1, 0.5);
							if (Color_option1 == 1)
							{
								if (offspring.z == 1)
									offspring.color = "purple"; // Color Based on Behavior
							}
						}
						if (map == 1)
						{
							// set offspring position if map = 1
							do
								pos = individual.spatialPosition + rnorm(2, 0, OMD);
							while (!p1.pointInBounds(pos) | p1.spatialMapValue("world", pos) == 0.0);
							offspring.setSpatialPosition(pos);
						}
						if (map != 1)
						{
							// set offspring position if map != 1
							do
								pos = individual.spatialPosition + rnorm(2, 0, OMD);
							while (!p1.pointInBounds(pos));
							offspring.setSpatialPosition(pos);
						}
					}
				}
			}
		}
	}
}
```

First the simulation looks for possible mates nearby and then the reproduction function is run.
This reproduction function runs for each individual each generation.
In the reproduction function the phenotype of new offspring is tagged with the Z coordinate of the individual. We also see this in the function where the individuals are initialized at the start of the simulation.
