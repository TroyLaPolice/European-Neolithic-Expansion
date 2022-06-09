# Modeling Hunter Gatherer/Farmer Interactions


## About the simulation

This simulation models the interaction between Hunter Gatherers and Farmers in Neolithic Europe.
It represents the expansion of farming throughout Europe during this time period.
It is a Non-WF model featuring spatial competition, reproduction and learning. It is a 2D spatial model that can be displayed on a map.

### About the map:

The map is a square approximation of Europe which by default is 3,700 km^2.

### Color schemes:

#### Behavioral coloring
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
 to
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
	defineConstant("SN", 10000); // Starting number of individuals
	defineConstant("HGK", 0.064); // carrying capacity for HGs (ENTER IN INDIVIDUALS PER KM2) for density dependent scaling 
	defineConstant("FK", 1.28); // carrying capacity for farmers (ENTER IN INDIVIDUALS PER KM2) for density dependent scaling
	
 ```
These parameters describe the starting number of individuals in the sim and the carrying capacity (given in individuals per km^2) for HGs and farmers.
These parameters should be adjusted to your liking.

```
	// Movement and interaction Distances (ENTER IN KILOMETERS):
	// ***********************************
	defineConstant("S", 50); // spatial competition distance
	defineConstant("MD", 30); // Mating distance
	defineConstant("LD", 20); // Learning distance
	defineConstant("OMD", 20); // Offspring movement away from parents. Right now OMD is the same for both populations but we could make it unique
	defineConstant("FMD", 10); // How far farmers diffuse away from their location
	defineConstant("HGMD", 70); // How far HGs diffuse away from their location
	
```

These parameters are all distance related and will be entered in km.

```
	// Learning, death and mating rate params:
	// ***********************************
	defineConstant("L", 0.009); // Learning rate 
	defineConstant("LP", 0.6); // Learning percentage = the ratio of farmers to HGs required in an area for an individual HG to learn from a farmer 
	defineConstant("M", 0.1); // Mating rate
	defineConstant("min_repro_age", 0); // Individuals MUST be OLDER than this age to reproduce
	
	// Age related mortality table
	defineConstant("age_scale", c(0.211180124, 0.211180124, 0.211180124, 0.211180124, 0.211180124, 0.251968504, 0.251968504, 0.251968504, 0.251968504, 0.251968504, 0.105263158, 0.105263158, 0.105263158, 0.105263158, 0.105263158, 0.164705882, 0.164705882, 0.164705882, 0.164705882, 0.164705882, 0.164705882, 0.253521127, 0.253521127, 0.253521127, 0.253521127, 0.253521127, 0.301886792, 0.301886792, 0.301886792, 0.301886792, 0.301886792, 0.378378378, 0.378378378, 0.378378378, 0.378378378, 0.378378378, 0.47826087, 0.47826087, 0.47826087, 0.47826087, 0.47826087, 0.583333333, 0.583333333, 0.583333333, 0.583333333, 0.583333333, 0.6, 0.6, 0.6, 0.6, 0.6, 1.0));
 ```
 
 These parameters are all about rates for mating, learning and death. This is the likelihood that farmers will teach a HG how to farm. The first param L is just simply how often this interaction happens.
 The second (LP) is the ratio of farmers to HGs in an area that is required for a HG to learn. 
 For example, if a given area is majority HG it may be less likely that the minority farmers teach. 
 This of course can be changed as you see fit. By default here we see the param is 60%.
 
 Next is a minimum age required for reproduction. Individuals MUST BE OLDER than the provided age in this parameter for them to be able to reproduce. This prevents infants and small children from being able to reproduce which is unrealistic.
 
 Lastly, the age related mortality table is a life table. This data comes from "age at death" studies. It is implemented similarly to SLiM recipe 16.2 Age structure (a life table model) by Benjamin C. Haller and Philipp W. Messer. It is described in the [manual](http://benhaller.com/slim/SLiM_Manual.pdf) in a succinct and useful way:
 
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

*First, we will look at coloring schemes. These are described at the top of the markdown.*

```
	// Determine Coloring Schemes:
	// ***********************************	
	defineConstant("Color_scheme", 1); // Parameter of 0 = Genomic Coloring, 1 = Behavioral Coloring
	defineConstant("Color_option1", 0); // Parameter of 0 = No special color for 1st generation 'hybrid' offspring, 1 = Coloring for 1st generation 'hybrid' offspring
	defineConstant("Color_option2", 1); // Parameter of 0 = No special color for 1st generation HG individuals who have learned, 1 = Coloring for 1st generation HG individuals who have learned
```
Color_scheme is a binary choice for if you want genomic coloring gradient or phenotypic behavioral coloring

The next two are options that determine if you want special colors for specific conditions like those who have recently learned or first generation "hybrids"

*Preferences regarding the map on which the individuals move*

```
	
	// Map Prefs:
	// ***********************************	

	defineConstant("mapfile", "C:/PATH/GOES/HERE/TO/MAP/europe.png"); // File Path to Map Image
	defineConstant("map_size_length_and_width", 3700);
```
The first parameter is a path in YOUR file system that the map file can be found in. The map is available for download in the repo. This path points to where the map is in your system so it MUST be changed to fit your file structure.

The next parameter is the length and width of the map. This will be used when building the map.

*Next are several parameters that **will not typically need to be changed** but some can be if you want*

```
	// ---------------------------------------------------
	//  PARAMETERS FOR MAP SCALING (NO NEED TO ADJUST)
	// ---------------------------------------------------
	defineConstant("map_size_km2", (map_size_length_and_width * map_size_length_and_width));
```

*The next two blocks simply initialize the genetic component that is the marker for farmer ancestry and the interactions between individuals*

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

*1) Spatial competition between nearby individuals*

*2) Mating*

*3) Learning*

These all take place within a certain distance range specified by parameters above.

#### First generation of the simulation- building initial population
```
1 early()
{
	// Set up map
	p1.setSpatialBounds(c(0.0, 0.0, map_size_length_and_width, map_size_length_and_width));
	mapImage = Image(mapfile);
	p1.defineSpatialMap("world", "xy", 1.0 - mapImage.floatK, valueRange=c(0.0, 1.0), colors=c("#ffffff", "#000000"));
	
	// start near a specific map location
	for (ind in p1.individuals)
	{
		do
			newPos = c(runif(1, 0, map_size_length_and_width), runif(1, 0, map_size_length_and_width));
		while (!p1.pointInBounds(newPos) | p1.spatialMapValue("world", newPos) == 0.0);
		ind.setSpatialPosition(newPos);
	}
	
```

This section sets up how the image file of the map of Europe works.

It defines the map and its bounds.

```
	// Define z param in offspring (phenotype, 0 = HG, 1 = Farmer)
	// Make individuals near anatolia and greece farmers
	p1.individuals[p1.individuals.x > (0.8 * map_size_length_and_width) & p1.individuals.y < (0.2 * map_size_length_and_width)].z = 1;
	
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

This portion sets up the z coordinate for individuals and adds the marker mutation.

The z coordinate here represents the behavioral phenotype- farming vs HGing.

Individuals located near Anatolia and Greece on the map begin as farmers and slowly expand throughout Europe.

This portion also sets up how the coloring schemes work. 

It sets it up both by phenotype (behavioral) coloring and genomic coloring based on user provided preference above.

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
	if (individual.age > min_repro_age) // Reproductive age of the individual must be reached before mating
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
						// Color Based on Genotype
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
						// they become a farmer
						if (individual.z != mate.z)
						{
							offspring.z = 1;
							if (Color_option1 == 1)
							{
								if (offspring.z == 1)
									offspring.color = "purple"; // Color Based on Behavior
							}
						}
						
						// set offspring position
						do
							pos = individual.spatialPosition + rnorm(2, 0, OMD);
						while (!p1.pointInBounds(pos) | p1.spatialMapValue("world", pos) == 0.0);
						offspring.setSpatialPosition(pos);
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

The offspring then moves away from its parents within a specified distance range. The code is slightly different given the map or square. This is because the individuals can't persist in the ocean.

Color is also specified here. It can be based on ancestry proportion or simply what the individual is behaviorally- HG or farmer.

#### Learning

```
late()
{
	// ---------------------------------------------------
	//  LEARNING --> HGs learn to farm from nearby farmers
	// ---------------------------------------------------
	i3.evaluate();
	for (individual in p1.individuals)
	{
		if (individual.z == 1)
			next;
		
		// Get a vector the nearest neighbors within the learning distance (D) 
		neighbors = i3.nearestNeighbors(individual);
		
		// Get ratio of HGs to farmers in the neighbors of the individual
		neighbor_freq = (sum(neighbors.z) / length(neighbors));
		
		// If the HG is surrounded by a certain ratio of farmers (LP) it has the ablity to convert to farming by learning
		if (neighbor_freq >= LP)
		{
			// choose nearest neighbor as a teacher, within the max distance
			teacher = i3.nearestNeighbors(individual, 1);
			
			// Frequency of the interaction
			for (i in seqLen(rpois(1, L)))
			{
				// Only runs if a potential teacher is nearby
				if (teacher.size())
				{
					// Binomial choice for teaching
					chance = rbinom(1, 1, 0.5);
					
					// If teacher is a farmer and individual is a HG
					if (teacher.z == 1 & chance == 1)
					{
						if (Color_option2 == 1)
						{
							// -----------------------
							// Color Based on Behavior
							// -----------------------
							// Change the HG to green to see the interaction
							// and change its phenotype (z coordinate) from 0 to 1
							// to represent the conversion to farmer
							// Only first generation HG -> F converts are green
							individual.color = "green";
						}
						individual.z = 2;
					}
				}
			}
		}
	}
	p1.individuals[p1.individuals.z == 2].z = 1;
}
```

Learning is implemented in a similar way to reproduction. Where it differs is that a specified ratio of farmers to HGs must be met for learning to happen.

Individual HGs can learn from farmers and their z coordinate changes but their ancestry stays the same.

#### Competition


First we evaluate the interaction type and define vectors of the two groups

```
early()
{
	i1.evaluate();
	
	// define vector of farmers and vector of HGs
	farmers = p1.individuals[p1.individuals.z == 1];
	HGs = p1.individuals[p1.individuals.z == 0];
```

This part handles competition between nearby individuals. This is density dependent. There can be different K's for HGs and farmers here which we will see later. First we count the number nearby competing individuals.

```
	// spatial competition provides density-dependent selection
	// Count number of neighbors within S for farmers and hunter gatherers
	farmers_num_in_s = i1.interactingNeighborCount(farmers);
	HG_num_in_s = i1.interactingNeighborCount(HGs);
```

This next part keeps individuals from living beyond realistic limits. Without this individuals in the sim can live hundreds of years because death it not dependent on age, only population density. (See life table above)

```
	// life table based individual mortality
	farmer_ages = farmers.age;
	HG_ages = HGs.age;
	
	// calculate chance of survival by refering to the mortality table by age and subtracting the chance of mortality from one
	farmer_survival = 1 - age_scale[farmer_ages];
	HG_survival = 1 - age_scale[HG_ages];
```

Finally we bring the two parts together.

```
	// density-dependence, factoring in individual mortality
	farmers.fitnessScaling = ((PI * (S^2) * FK) / (farmers_num_in_s + 1) * farmer_survival);
	HGs.fitnessScaling = ((PI * (S^2) * HGK) / (HG_num_in_s + 1) * HG_survival);
```


#### Movement of individuals

```
late()
{
	// move around
	for (ind in p1.individuals)
	{
		if (ind.z == 1)
		{
			// How far farmers diffuse away from their location			
			do
				newPos = ind.spatialPosition + runif(2, -FMD, FMD);
			while (!p1.pointInBounds(newPos) | p1.spatialMapValue("world", newPos) == 0.0);
			ind.setSpatialPosition(newPos);
		}
		if (ind.z == 0)
		{
			// How far HGs diffuse away from their location
			do
				newPos = ind.spatialPosition + runif(2, -HGMD, HGMD);
			while (!p1.pointInBounds(newPos) | p1.spatialMapValue("world", newPos) == 0.0);
			ind.setSpatialPosition(newPos);
		}
	}
}
```

The individuals cannot move to locations outside of the bounds of the map. They can potentially jump across the water (simulating water travel) to other land masses, assuming it is not beyond their movement range, but they cannot stay in the ocean.

The individuals can have different distances they can travel based on if they are a HG or a farmer. This is set up above in the parameters. This allows for simulation of HGs being more migratory and farmers being more localized around their farm.
