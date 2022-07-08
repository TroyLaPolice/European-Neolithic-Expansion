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

*First we need to set our working directory and if the user choses a special custom map they will need to provide the file name (otherwise leave this as 0)*

```
	//SET WORKING DIRECTORY AND CUSTOM MAP NAME IF DESIRED
	defineConstant("wd", "~/PATH/TO/FILE/HERE");
	defineConstant("custom_map_filename", 0); // Parameter = file name of map (as string) if user wants to use their own map and override the built in maps, else == 0
	defineConstant("output_name", "run_name_goes_here"); // To run default names make an empty string
```

*Next begins the set up where parameters for the model will be set. Some of these will not need to be touched but some should be tuned to your liking.*

###### Parameters for how the population grows, interacts and moves

```
  	// Carrying Capacities and Pop Sizes:
	// ***********************************
	defineConstant("SN", 8761600); // Starting number of individuals
	defineConstant("HGK", 0.64); // carrying capacity for HGs (ENTER IN INDIVIDUALS PER KM2) for density dependent scaling 
	defineConstant("FK", 1.28); // carrying capacity for farmers (ENTER IN INDIVIDUALS PER KM2) for density dependent scaling
	
 ```
The following parameters below describe the starting number of individuals in the sim for HGs and farmers and the distances they can travel.
These parameters should be adjusted to your liking.

The movement_distances parameter is a series of distances that are sampled from based on the probability that an individual travels that distance. The movement_distance_weights parameter is the vector of corresponding probabilities.

```
	// Movement and interaction Distances (ENTER IN KILOMETERS):
	// ***********************************
	defineConstant("S", 30); // spatial competition distance
	defineConstant("MD", 30); // Mating distance
	defineConstant("movement_distances", c(2.3, 7.3, 15, 25, 35, 45, 55, 100)); // Distances sampled from
	defineConstant("movement_distance_weights", c(0.42, 0.23, 0.16, 0.08, 0.07, 0.02,
0.01, 0.01)); // Weights for movement distance sampling
	defineConstant("LD", 10); // Learning distance
	defineConstant("northern_slowdown_effect", 2); // Number equals the effect of the slowdown in the north (i.e., how many times slower do they move
	defineConstant("northern_slowdown_distance", 0.3); // How far north on the map does the slow down begin (0.5 is mid point, scale is 0-1 with far South being 0 and North being 1)
	
```

These parameters are all distance related and will be entered in km.

```
	// Learning, death and mating rate params:
	// ***********************************
	defineConstant("L", 0.1); // Learning rate 
	defineConstant("LP", 0.6); // Learning percentage = the ratio of farmers to HGs required in an area for an individual HG to learn from a farmer 
	defineConstant("HGM", 0.1); // HG fertility rate
	defineConstant("FM", 0.1); // Farmer fertility rate
	defineConstant("IM", 0.01); // Interbreeding fertility rate
	defineConstant("min_repro_age", 0); // Individuals MUST be OLDER than this age to reproduce
	// Age related mortality table
	defineConstant("age_scale", c(0.211180124, 0.211180124, 0.211180124, 0.211180124, 0.211180124, 0.251968504, 0.251968504, 0.251968504, 0.251968504, 0.251968504, 0.105263158, 0.105263158, 0.105263158, 0.105263158, 0.105263158, 0.164705882, 0.164705882, 0.164705882, 0.164705882, 0.164705882, 0.164705882, 0.253521127, 0.253521127, 0.253521127, 0.253521127, 0.253521127, 0.301886792, 0.301886792, 0.301886792, 0.301886792, 0.301886792, 0.378378378, 0.378378378, 0.378378378, 0.378378378, 0.378378378, 0.47826087, 0.47826087, 0.47826087, 0.47826087, 0.47826087, 0.583333333, 0.583333333, 0.583333333, 0.583333333, 0.583333333, 0.6, 0.6, 0.6, 0.6, 0.6, 1.0));
 ```
 
 These parameters are all about rates for mating, learning and death. 
 
 The first perameter L refers to likelihood that farmers will teach a HG how to farm.
 
 The second (LP) is the ratio of farmers to HGs in an area that is required for a HG to learn. For example, if a given area is majority HG it may be less likely that the minority farmers teach. This of course can be changed as you see fit. By default here we see the param is 60%.
 
 The next block of parameters handles reproduction probabilities (i.e., the probability that a pair mates and produces an offspring) 
 
 HGM is this rate for a pairing of two hunter gatherers will produce an offspring <HGM * 100> percent of the time. 
 
 The other parameters are the same idea- FM is the mating probability for two farmers.
 
 IM is how often a farmer and hunter gatherer interbreed and produce an offspring (controls assortative mating)
 
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
	defineConstant("map_style", 5); // Parameter of 0 = no topography, 1 = super light topopgraphy, 2 = light topography, 3 = regular topography, 4 heavy topography, 5 = square, 6 = custom map
	
	defineConstant("water_crossings", 1); //Parameter of 0 = no water crossing paths, 1 = water crossing paths
	
	if (water_crossings == 1)
		defineConstant("file_extention", "_water.png");
	else
		defineConstant("file_extention", ".png");
```

*These parameters do not need to be changed! They handle what map to use. If you want to use a custom map, enter the file name at the begining of the script*

```
	defineConstant("mapfile_none", wd + "/EEA_map" + file_extention); // File Path to Map Image
	defineConstant("mapfile_topo_superlight", wd + "/EEA_map_topo_superlight" + file_extention); // File Path to Map Image
	defineConstant("mapfile_topo_light", wd + "/EEA_map_topo_light" + file_extention); // File Path to Map Image
	defineConstant("mapfile_topo_regular", wd + "/EEA_map_topo_regular" + file_extention); // File Path to Map Image
	defineConstant("mapfile_topo_heavy", wd + "/EEA_map_topo_heavy" + file_extention); // File Path to Map Image
	defineConstant("square", wd + "/square.png"); // File Path to Map Image
	defineConstant("custom_map", wd + custom_map_filename); // File Path to Map Image
```


These parameters handle the map size in km. If using provided EEA maps, length and width for map style images should be kept at a 1 to 1 aspect ratio (square) to avoid distortion. By default they are 3700 x 3700 km square maps

```
	defineConstant("map_size_length", 3700);
	defineConstant("map_size_width", 3700);
	
	// Map citation: https://www.eea.europa.eu/data-and-maps/figures/elevation-map-of-europe
```
The first parameter is a path in YOUR file system that the map file can be found in. The maps are available for download in the repo. This path points to where the map is in your system so it MUST be changed to fit your file structure.

The next parameter is the length and width of the map. This will be used when building the map.

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
	// Check user input for what style of topogrpahy they want on the map
	if (map_style == 0)
		mapImage = Image(mapfile_none); //none
	else if (map_style == 1)
		mapImage = Image(mapfile_topo_superlight); //superlight topo
	else if (map_style == 2)
		mapImage = Image(mapfile_topo_light); //light
	else if (map_style == 3)
		mapImage = Image(mapfile_topo_regular); //regular
	else if (map_style == 4)
		mapImage = Image(mapfile_topo_heavy); //heavy
	else if (map_style == 5)
		mapImage = Image(square); //square
	else if (map_style == 6)
		mapImage = Image(custom_map); //square
	
	// Set up map
	p1.setSpatialBounds(c(0.0, 0.0, map_size_width, map_size_length));
	p1.defineSpatialMap("map_object", "xy", 1.0 - mapImage.floatK, valueRange=c(0.0, 1.0), colors=c("#ffffff", "#000000"));
	
	// place individuals on the map
	for (ind in p1.individuals)
	{
		do
			newPos = c(runif(1, 0, map_size_width), runif(1, 0, map_size_length));
		while (!p1.pointInBounds(newPos) | p1.spatialMapValue("map_object", newPos) == 0.0);
		ind.setSpatialPosition(newPos);
	}
	
```

This section sets up how the image file of the map of Europe works.

It defines the map and its bounds.

```
	// Define z param in offspring (phenotype, 0 = HG, 1 = Farmer)
	// Make individuals near anatolia and greece farmers or near the edge if using the square
	//
	if (map_style == 5)
		p1.individuals[p1.individuals.x > (0.98 * map_size_width)].z = 1;
	else
		p1.individuals[p1.individuals.x > (0.72 * map_size_width) & p1.individuals.y < (0.09 * map_size_length)].z = 1;
	
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
				if (individual.z != mate.z)
					M = IM;
				if (individual.z == mate.z & individual.z == 1)
					M = FM;
				if (individual.z == mate.z & individual.z == 0)
					M = HGM;
				
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
						// offspring becomes a farmer
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
						if (map_style != 5 & individual.y > northern_slowdown_distance * map_size_length)
						{
							do
							{
								// This samples from a vector of movement distances based on the probability that they move this distance
								distance = sample(x = c(movement_distances), size = 1, replace = T, weights = c(movement_distance_weights));
								// Next we need to calculate the x and y coodinates
								radian_angle = runif(1, 0, 2*PI);
								coordiates = c(cos(radian_angle) * distance, sin(radian_angle) * distance) / northern_slowdown_effect;
								// Next we can reset the position
								pos = individual.spatialPosition + coordiates;
							}
							while (!p1.pointInBounds(pos) | p1.spatialMapValue("map_object", pos) == 0.0);
							offspring.setSpatialPosition(pos);
						}
						else
						{
							do
							{
								// This samples from a vector of movement distances based on the probability that they move this distance
								distance = sample(x = c(movement_distances), size = 1, replace = T, weights = c(movement_distance_weights));
								// Next we need to calculate the x and y coodinates
								radian_angle = runif(1, 0, 2*PI);
								coordiates = c(cos(radian_angle) * distance, sin(radian_angle) * distance);
								// Next we can reset the position
								pos = individual.spatialPosition + coordiates;
							}
							while (!p1.pointInBounds(pos) | p1.spatialMapValue("map_object", pos) == 0.0);
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
					// If teacher is a farmer and individual is a HG
					if (teacher.z == 1)
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
						individual.z = 1;
					}
				}
			}
		}
	}
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
	// Count number of neighbors within S for farmers
	farmers_num_in_s = i1.interactingNeighborCount(farmers);
	
	// Do the same for HGs if there are still HGs left
	if (length(HGs) != 0)
		HG_num_in_s = i1.interactingNeighborCount(HGs);
	
	// Life table based individual mortality, get vector of individual ages
	farmer_ages = farmers.age;
	
	// Do the same for HGs if there are still HGs left
	if (length(HGs) != 0)
		HG_ages = HGs.age;
```

This next part keeps individuals from living beyond realistic limits. Without this individuals in the sim can live hundreds of years because death it not dependent on age, only population density. (See life table above)

```
	// Modify mortality curve to account for population density around the indiviudal
	scaled_mortality_farmer = (farmers_num_in_s + 1) / (PI * (S^2) * FK + 1) * age_scale[farmer_ages];
	
	// Do the same if there are still HGs left
	if (length(HGs) != 0)
		scaled_mortality_HG = (HG_num_in_s + 1) / (PI * (S^2) * HGK + 1) * age_scale[HG_ages];
	
	// Set a maximum age and make sure there are no negative fittnesses
	scaled_mortality_farmer[farmer_ages == length(age_scale) - 1] = 1;
	scaled_mortality_farmer[scaled_mortality_farmer > 1] = 1;
	
	// Do the same for HGs if there are still HGs left
	if (length(HGs) != 0)
	{
		scaled_mortality_HG[HG_ages == length(age_scale) - 1] = 1;
		scaled_mortality_HG[scaled_mortality_HG>1] = 1;
	}
	
	// Calculate chance of survival by refering to the mortality table by age and subtracting the chance of mortality from one
	farmer_survival = 1 - scaled_mortality_farmer;
	
		
	// Do the same if there are still HGs left
	if (length(HGs) != 0)
		HG_survival = 1 - scaled_mortality_HG;
```

Finally we scale the individuals' fitness by the calculated value.

```
	// Scale the fitneess of the individual, density-dependent and factoring in individual age based mortality
	farmers.fitnessScaling = farmer_survival;
		
	// Do the same if there are still HGs left
	if (length(HGs) != 0)
		HGs.fitnessScaling = HG_survival;
}
```


#### Movement of individuals

```
late()
{
	// move around
	for (ind in p1.individuals)
	{
		// How far farmers diffuse away from their location
		if (ind.z == 1)
		{
			if (map_style != 5 & ind.y > northern_slowdown_distance * map_size_length)
			{
				do
				{
					// This samples from a vector of movement distances based on the probability that they move this distance
					distance = sample(x = c(movement_distances), size = 1, replace = T, weights = c(movement_distance_weights));
					// Next we need to calculate the x and y coodinates
					radian_angle = runif(1, 0, 2*PI);
					coordiates = c(cos(radian_angle) * distance, sin(radian_angle) * distance) / northern_slowdown_effect;
					// Next we can reset the position
					newPos = ind.spatialPosition + coordiates;
				}
				while (!p1.pointInBounds(newPos) | p1.spatialMapValue("map_object", newPos) == 0.0);
				ind.setSpatialPosition(newPos);
			}
			else
			{
				do
				{
					// This samples from a vector of movement distances based on the probability that they move this distance
					distance = sample(x = c(movement_distances), size = 1, replace = T, weights = c(movement_distance_weights));
					// Next we need to calculate the x and y coodinates
					radian_angle = runif(1, 0, 2*PI);
					coordiates = c(cos(radian_angle) * distance, sin(radian_angle) * distance);
					// Next we can reset the position
					newPos = ind.spatialPosition + coordiates;
				}
				while (!p1.pointInBounds(newPos) | p1.spatialMapValue("map_object", newPos) == 0.0);
				ind.setSpatialPosition(newPos);
			}
		}
		if (ind.z == 0)
		{
			// How far HGs diffuse away from their location
			do
			{
				// This samples from a vector of movement distances based on the probability that they move this distance
				distance = sample(x = c(movement_distances), size = 1, replace = T, weights = c(movement_distance_weights));
				// Next we need to calculate the x and y coodinates
				radian_angle = runif(1, 0, 2*PI);
				coordiates = c(cos(radian_angle) * distance, sin(radian_angle) * distance);
				// Next we can reset the position
				newPos = ind.spatialPosition + coordiates;
			}
			while (!p1.pointInBounds(newPos) | p1.spatialMapValue("map_object", newPos) == 0.0);
			ind.setSpatialPosition(newPos);
		}
	}
}
```

The individuals cannot move to locations outside of the bounds of the map. They can potentially jump across the water (simulating water travel) to other land masses, assuming it is not beyond their movement range, but they cannot stay in the ocean.

The individuals can have different distances they can travel based on if they are a HG or a farmer. This is set up above in the parameters. This allows for simulation of HGs being more migratory and farmers being more localized around their farm.

Of course if you chose to run the simulation with the simple black square the individuals can move anywhere within the given map size.

#### Finally, write output to files

```
1 early()
{
	// log runtime params
	param_string = paste(SN, HGK, FK, S, MD, LD, northern_slowdown_effect, northern_slowdown_distance, L, LP, HGM, FM, IM, min_repro_age, map_style, water_crossings, "[", age_scale, "]", "[", movement_distances, "]", "[", movement_distance_weights, "]");
	
	// File headings
	param_heading = paste("SN HGK FK S MD LD northern_slowdown_effect northern_slowdown_distance L LP HGM FM IM min_repro_age map_style water_crossings [ age_scale ]  [ movement_distances ]  [ movement_distance_weights ]");
	
	// Runtime params - write to file
	output_runtime_file_name = ("/sim_runtime_params_" + output_name + ".txt");
	writeFile(wd + output_runtime_file_name, param_heading, append=T);
	writeFile(wd + output_runtime_file_name, param_string, append=T);
	
	// Population stats headers - write to file
	stats_header_string = paste("Year", "PopulationSize", "TotalFarmers", "TotalHGs", "RatioFarmertoHG", "NewBirths");
	output_stats_file_name = ("/sim_pop_stats_per_year_" + output_name + ".txt");
	writeFile(wd + output_stats_file_name, stats_header_string, append=T);
	
	if (map_style == 5)
	{
		// Wave stats headers - write to file
		wave_stats_header_string = paste("Year", "PopulationSize", "TotalFarmers", "TotalHGs", "RatioFarmerToHG", "Farmers_in_Partition1", "Farmers_in_Partition2", "Farmers_in_Partition3", "Farmers_in_Partition4", "Farmers_in_Partition5", "Farmers_in_Partition6", "Farmers_in_Partition7", "Farmers_in_Partition8", "Farmers_in_Partition9", "Farmers_in_Partition10", "HGs_in_Partition1", "HGs_in_Partition2", "HGs_in_Partition3", "HGs_in_Partition4", "HGs_in_Partition5", "HGs_in_Partition6", "HGs_in_Partition7", "HGs_in_Partition8", "HGs_in_Partition9", "HGs_in_Partition10", "RatioFarmerToHG_Partition1", "RatioFarmerToHG_Partition2", "RatioFarmerToHG_Partition3", "RatioFarmerToHG_Partition4", "RatioFarmerToHG_Partition5", "RatioFarmerToHG_Partition6", "RatioFarmerToHG_Partition7", "RatioFarmerToHG_Partition8", "RatioFarmerToHG_Partition9", "RatioFarmerToHG_Partition10", "NewBirths");
		output_wave_stats_file_name = ("/sim_square_wave_stats_per_year_" + output_name + ".txt");
		writeFile(wd + output_wave_stats_file_name, wave_stats_header_string, append=T);
	}
}

1:8000 late()
{
	// provide feedback on progress for command line users
	year_counter = paste("Simulation Year: ", sim.generation);
	print(year_counter);
	
	if (sim.generation == 8000)
	{
		print("--------------------------------");
		print("Simulation Complete");
		print("--------------------------------");
	}
	
	// calculate num farmers
	num_farmers = sum(p1.individuals.z);
	
	//calculate num HGs
	num_HGs = (p1.individuals.length() - sum(p1.individuals.z));
	
	// calculate the ratio of farmers in the total population to file
	ratio = (sum(p1.individuals.z) / p1.individuals.length());
	
	// calculate population size statistics
	pop_size = p1.individuals.length();
	
	// calculate the number of new births this year
	new_births = length(p1.individuals.age == 0);
			
	// write outputs
	output_string = paste(sim.generation, pop_size, num_farmers, num_HGs, ratio, new_births);
	output_stats_file_name = ("/sim_pop_stats_per_year_" + output_name + ".txt");
	writeFile(wd + output_stats_file_name, output_string, append=T);
}

```
The following block of code runs of the sim is run on a square rather than a map and provides more detailed outputs regarding the wave progression
```
late()
{
	if (map_style == 5)
	{
		// define vector of farmers and vector of HGs
		farmers = p1.individuals[p1.individuals.z == 1];
		HGs = p1.individuals[p1.individuals.z == 0];
		
		// Split width into 10 equal parts
		partition_widths = map_size_width / 10;
		
		// set criteria for farmers in each partition
		farmers_partition1_bool = farmers.x <= partition_widths;
		farmers_partition2_bool = farmers.x <= 2 * partition_widths & farmers.x > partition_widths;
		farmers_partition3_bool = farmers.x <= 3 * partition_widths & farmers.x > 2 * partition_widths;
		farmers_partition4_bool = farmers.x <= 4 * partition_widths & farmers.x > 3 * partition_widths;
		farmers_partition5_bool = farmers.x <= 5 * partition_widths & farmers.x > 4 * partition_widths;
		farmers_partition6_bool = farmers.x <= 6 * partition_widths & farmers.x > 5 * partition_widths;
		farmers_partition7_bool = farmers.x <= 7 * partition_widths & farmers.x > 6 * partition_widths;
		farmers_partition8_bool = farmers.x <= 8 * partition_widths & farmers.x > 7 * partition_widths;
		farmers_partition9_bool = farmers.x <= 9 * partition_widths & farmers.x > 8 * partition_widths;
		farmers_partition10_bool = farmers.x <= 10 * partition_widths & farmers.x > 9 * partition_widths;
		
		// count farmers in each partition 
		farmers_partition1 = size(which(farmers_partition1_bool == T));
		farmers_partition2 = size(which(farmers_partition2_bool == T));
		farmers_partition3 = size(which(farmers_partition3_bool == T));
		farmers_partition4 = size(which(farmers_partition4_bool == T));
		farmers_partition5 = size(which(farmers_partition5_bool == T));
		farmers_partition6 = size(which(farmers_partition6_bool == T));
		farmers_partition7 = size(which(farmers_partition7_bool == T));
		farmers_partition8 = size(which(farmers_partition8_bool == T));
		farmers_partition9 = size(which(farmers_partition9_bool == T));
		farmers_partition10 = size(which(farmers_partition10_bool == T));
		
		// set criteria for HGs in each partition
		HGs_partition1_bool = HGs.x <= partition_widths;
		HGs_partition2_bool = HGs.x <= 2 * partition_widths & HGs.x > partition_widths;
		HGs_partition3_bool = HGs.x <= 3 * partition_widths & HGs.x > 2 * partition_widths;
		HGs_partition4_bool = HGs.x <= 4 * partition_widths & HGs.x > 3 * partition_widths;
		HGs_partition5_bool = HGs.x <= 5 * partition_widths & HGs.x > 4 * partition_widths;
		HGs_partition6_bool = HGs.x <= 6 * partition_widths & HGs.x > 5 * partition_widths;
		HGs_partition7_bool = HGs.x <= 7 * partition_widths & HGs.x > 6 * partition_widths;
		HGs_partition8_bool = HGs.x <= 8 * partition_widths & HGs.x > 7 * partition_widths;
		HGs_partition9_bool = HGs.x <= 9 * partition_widths & HGs.x > 8 * partition_widths;
		HGs_partition10_bool = HGs.x <= 10 * partition_widths & HGs.x > 9 * partition_widths;
		
		// count HGs in each partition 
		HGs_partition1 = size(which(HGs_partition1_bool == T));
		HGs_partition2 = size(which(HGs_partition2_bool == T));
		HGs_partition3 = size(which(HGs_partition3_bool == T));
		HGs_partition4 = size(which(HGs_partition4_bool == T));
		HGs_partition5 = size(which(HGs_partition5_bool == T));
		HGs_partition6 = size(which(HGs_partition6_bool == T));
		HGs_partition7 = size(which(HGs_partition7_bool == T));
		HGs_partition8 = size(which(HGs_partition8_bool == T));
		HGs_partition9 = size(which(HGs_partition9_bool == T));
		HGs_partition10 = size(which(HGs_partition10_bool == T));
		
		// calculate the ratio of farmers each partition
		ratio1 = (farmers_partition1 / (farmers_partition1 + HGs_partition1));
		ratio2 = (farmers_partition2 / (farmers_partition2 + HGs_partition2));
		ratio3 = (farmers_partition3 / (farmers_partition3 + HGs_partition3));
		ratio4 = (farmers_partition4 / (farmers_partition4 + HGs_partition4));
		ratio5 = (farmers_partition5 / (farmers_partition5 + HGs_partition5));
		ratio6 = (farmers_partition6 / (farmers_partition6 + HGs_partition6));
		ratio7 = (farmers_partition7 / (farmers_partition7 + HGs_partition7));
		ratio8 = (farmers_partition8 / (farmers_partition8 + HGs_partition8));
		ratio9 = (farmers_partition9 / (farmers_partition9 + HGs_partition9));
		ratio10 = (farmers_partition10 / (farmers_partition10 + HGs_partition10));
	
		// calculate total num farmers
		num_farmers = sum(p1.individuals.z);
	
		// calculate total num HGs
		num_HGs = (p1.individuals.length() - sum(p1.individuals.z));
	
		// calculate population size statistics
		pop_size = p1.individuals.length();
		
		// calculate the ratio of farmers in the total population to file
		ratio = (sum(p1.individuals.z) / p1.individuals.length());
		
		// calculate the number of new births this year
		new_births = length(p1.individuals.age == 0);
	
		// write outputs
		output_string = paste(sim.generation, pop_size, num_farmers, num_HGs, ratio, farmers_partition1, farmers_partition2, farmers_partition3, farmers_partition4, farmers_partition5, farmers_partition6, farmers_partition7, farmers_partition8, farmers_partition9, farmers_partition10, HGs_partition1, HGs_partition2, HGs_partition3, HGs_partition4, HGs_partition5, HGs_partition6, HGs_partition7, HGs_partition8, HGs_partition9, HGs_partition10, ratio1, ratio2, ratio3, ratio4, ratio5, ratio6, ratio7, ratio8, ratio9, ratio10, new_births);
		
		// output to file
		output_stats_file_name = ("/sim_square_wave_stats_per_year_" + output_name + ".txt");
		writeFile(wd + output_stats_file_name, output_string, append=T);
		
	}
}
```
