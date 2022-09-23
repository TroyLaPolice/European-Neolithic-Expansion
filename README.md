# Modeling Hunter Gatherer/Farmer Interactions


## About the simulation

This simulation models the interaction between Hunter Gatherers and Farmers in Neolithic Europe.
It represents the expansion of farming throughout Europe during this time period.
It is a Non-WF model featuring spatial competition, reproduction and learning. It is a 2D spatial model that can be displayed on a map.

### About the map:

The map is a square approximation of Europe which, by default is 3,700 km^2.

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

###### These can be altered via the command line by using the -d flag and then setting the parameter equal to a value

i.e. -d HGK=0.60

```
  	// ---------------------------------------------------
	//  PARAMETERS --> Initialize Constant Params
	// ---------------------------------------------------
```

*First we need to set our working directory and if the user choses a special custom map they will need to provide the file name (otherwise leave this as 0)*

```
	//SET WORKING DIRECTORY AND CUSTOM MAP NAME IF DESIRED
	if (!exists("wd"))
		defineConstant("wd", "~/PATH/TO/FILE/HERE");
	if (!exists("custom_map_filename"))
		defineConstant("custom_map_filename", 0); // Parameter = file name of map (as string) if user wants to use their own map and override the built in maps, else == 0
	if (!exists("output_name"))
		defineConstant("output_name", "run_name_goes_here"); // To run default names make an empty string
```

*Next begins the set up where parameters for the model will be set. Some of these will not need to be touched but some should be tuned to your liking.*

###### Parameters for how the population grows, interacts and moves

```
  	// Carrying Capacities and Pop Sizes:
	// ***********************************
	if (!exists("SN"))
		defineConstant("SN", 876160); // Starting number of individuals
	if (!exists("HGK"))
		defineConstant("HGK", 0.064); // carrying capacity for HGs (ENTER IN INDIVIDUALS PER KM2) for density dependent scaling 
	if (!exists("FK"))
		defineConstant("FK", 0.128); // carrying capacity for farmers (ENTER IN INDIVIDUALS PER KM2) for density dependent scaling
	
 ```
The following parameters below describe the starting number of individuals in the sim for HGs and farmers and the distances they can travel.
These parameters should be adjusted to your liking.

The movement_distances parameter is a series of distances that are sampled from based on the probability that an individual travels that distance. The movement_distance_weights parameter is the vector of corresponding probabilities.

```
	// Movement and interaction and competition:
	// ***********************************
	if (!exists("S"))
		defineConstant("S", 30); // spatial competition distance (ENTER IN KILOMETERS)
	if (!exists("C"))
		defineConstant("C", 1); // If C = 1 that means there will be competition between individuals of different phenotypes (ie HGs compete with Farmers for space) If C = 0 then competition only happens within their like groups
	if (!exists("SDX"))
		defineConstant("SDX", 30); // Movement standard diviation (sigma) for distribution of distances sampled from in the x direction (ENTER IN KILOMETERS)
	if (!exists("SDY"))
		defineConstant("SDY", 30); // Movement standard diviation (sigma) for distribution of distances sampled from in the y direction (ENTER IN KILOMETERS)
		
```
 This next block contains parameters about rates for mating, learning and death. 
 	
```
/ Learning, death and mating rate params:
	// ***********************************
	if (!exists("L"))
		defineConstant("L", 0.1); // Learning rate - product of constant number of teachers times probability of learning per contact
	if (!exists("gamma"))
		defineConstant("gamma", 1); // preference of farmer teacher over HG
	if (!exists("LD"))
	defineConstant("LD", 30); // Learning distance (ENTER IN KILOMETERS)
	if (!exists("MD"))
		defineConstant("MD", 30); // Mating distance (ENTER IN KILOMETERS)	
	if (!exists("MP"))
		defineConstant("MP", 0.5); // Level of assortative mating. Probability that mates will choose behaviorally similar mates. 1 = full assortative, 0.5 = no assortative, 0 = individuals only mate with individuals of opposite phenotypes
	if (!exists("min_repro_age"))
		defineConstant("min_repro_age", 11); // Individuals MUST be OLDER than this age to reproduce
	
	// Age related mortality table
	defineConstant("age_scale", c(0.211180124, 0.211180124, 0.211180124, 0.211180124, 0.211180124, 0.251968504, 0.251968504, 0.251968504, 0.251968504, 0.251968504, 0.105263158, 0.105263158, 0.105263158, 0.105263158, 0.105263158, 0.164705882, 0.164705882, 0.164705882, 0.164705882, 0.164705882, 0.164705882, 0.253521127, 0.253521127, 0.253521127, 0.253521127, 0.253521127, 0.301886792, 0.301886792, 0.301886792, 0.301886792, 0.301886792, 0.378378378, 0.378378378, 0.378378378, 0.378378378, 0.378378378, 0.47826087, 0.47826087, 0.47826087, 0.47826087, 0.47826087, 0.583333333, 0.583333333, 0.583333333, 0.583333333, 0.583333333, 0.6, 0.6, 0.6, 0.6, 0.6, 1.0));
 ```
 
 The first perameter L refers to likelihood that farmers will teach a HG how to farm is is the product of constant number of teachers times probability of learning per contact.
 
 The second (gamma) is the preference of a HG to choose a farmer as a teacher over choosing a HG
 
 The next block of parameters handles reproduction probabilities (i.e., the probability that a pair mates and produces an offspring) 
 
 MP is the amount of assortative mating, i.e., the probability that mates will choose behaviorally similar mates. 1 = full assortative, 0.5 = no assortative, 0 = individuals only mate with individuals of opposite phenotypes
 
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
	if (!exists("map_style"))
		defineConstant("map_style", 5); // Parameter of 0 = no topography, 1 = super light topopgraphy, 2 = light topography, 3 = regular topography, 4 heavy topography, 5 = square, 6 = custom map
	if (!exists("num_partitions"))
		defineConstant("num_partitions", 20); // Number of partitions if using the square
	if (!exists("water_crossings"))
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
The maps are available for download in the repo. The path points to where the map is in your system so the name of the map MUST be changed to fit your file structure if you altered the name of the map.

The next parameter is the length and width of the map. This will be used when building the map.

These parameters handle the map size in km. If using provided EEA maps, length and width for map style images should be kept at a 1 to 1 aspect ratio (square) to avoid distortion. By default they are 3700 x 3700 km square maps

```
	if (!exists("map_size_length"))
		defineConstant("map_size_length", 3700);
	if (!exists("map_size_width"))
		defineConstant("map_size_width", 3700);;
	
	// Map citation: https://www.eea.europa.eu/data-and-maps/figures/elevation-map-of-europe
```


*The next two blocks initialize the genetic component that is the marker for farmer ancestry and the interactions between individuals*

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
	sim.addSubpop("p1", SN);
}

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
The next section is responsible for setting up the age distribution from the initial population. If the age stucture is not defined then there will be a mating depression and a population crash due to the unrealsitic distributuion of ages in the population. We want the population to be representitive from the first generation. There would never be a real life situation where inidividuals all spawned in at zero years old by themselves.

```
	// Set ages of created individuals (inital age dist)
	prob_of_surv_until = cumProduct(1 - age_scale); // Prob of survival until a certain age
	eq_age_distr = prob_of_surv_until / sum(prob_of_surv_until); // Equilibrium age distribution
	num_of_inds = SN; // get random ages for individuals sampled from age distribution
	eq_ages_for_sample = sample(x=0:(length(prob_of_surv_until) - 1), size=num_of_inds, replace=T, weights=eq_age_distr); // samples from pop to get ages
	p1.individuals.age = eq_ages_for_sample; // sets population ages
	eq_offspring_rate = 1 / sum(prob_of_surv_until[(min_repro_age):(length(prob_of_surv_until) - 1)]); // calculates how often individuals of repro age have a child
	defineConstant("M", eq_offspring_rate); // Define this value for use later
```

This section sets up how the image file of the map works.

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
2: first()
{
	// look for mates
	i2.evaluate(p1);
}

reproduction()
{
	// ---------------------------------------------------
	//  MATING --> Individuals mate with those close by
	// ---------------------------------------------------
	// choose nearest neighbor as a mate, within the max distance
	if (individual.age > min_repro_age) // Reproductive age of the individual must be reached before mating
	{
		mates = i2.nearestInteractingNeighbors(individual, p1.individuals.length());
		mates_ra = mates[mates.age > min_repro_age];
		if (mates_ra.size()) // Runs if there are possible mates near by
		{
			// Sets up a vector that checks if the individuals in the surrounding area are of like phenotype (i.e. if their z params are the same)
			mates_similarity_bool = mates_ra.z == individual.z;
			mates_similarity = asFloat(mates_similarity_bool);
			
			// Generates a vector of probabilites of the individual mating with a particular mate given their phenotype
			// If assortative mating is complete, HGs will only mate with other HGs and farmers only with farmers
			// If there is no assortative mating there will be an equal proability for all individuals
			mates_similarity[mates_similarity_bool == T] = MP;
			mates_similarity[mates_similarity_bool == F] = (1 - MP);
			
			// If there are individuals in the neary area that there is a chance the individual mates with it will run this section
			if (sum(mates_similarity) != 0)
			{
				// Samples an individual in the population based on the assortative mating probabilities. 
				// If no assortative mating all indivduals have an equal chance of getting selected
				mate = sample(mates_ra, 1, weights=mates_similarity);
				
				// Frequency of the interaction is based on the calculated reproduction value given by the mortality curve (below)
				for (i in seqLen(rpois(1, M)))
				{
					// Only runs if a mate is nearby
					if (mate.size())
					{
						// Generates an offspring
						offspring = subpop.addCrossed(individual, mate);
						
						// -----------------------
						// Color Based on Genotype
						// -----------------------
						value = offspring.countOfMutationsOfType(m1) / (sim.chromosome.lastPosition * 2);
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
						
						// Next we set the position of the offspring at the parent's location
						pos = individual.spatialPosition;
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
In the reproduction function the phenotype of new offspring is tagged with the Z coordinate of the individual. We also see this in the function where the individuals are initialized at the start of the simulation. The offspring is generated on the map at the parental location. 

Color is also specified here. It can be based on ancestry proportion or simply what the individual is behaviorally- HG or farmer.

Assortative mating is handled here as well. Based on the specified amount of assortative mating above (MP) the probabability of the individual choosing a similar mate is calculated and used as weights for the sampling function that selects a mate for the individual.

#### Learning

```
late()
{
	// ---------------------------------------------------
	//  LEARNING --> HGs learn to farm from nearby farmers
	// ---------------------------------------------------
	i3.evaluate(p1);
	for (individual in p1.individuals)
	{
		// Skip if the individual is a farmer and therefore doesn't learn farming as it is already practicing the behavior.
		if (individual.z == 1)
			next;
		
		// Get a vector the nearest neighbors within the learning distance (D) 
		neighbors = i3.nearestNeighbors(individual, p1.individuals.length());
		
		// Determine z for neighbors of the individual
		farmer_neighbors = sum(neighbors.z);
		HG_neighbors = (length(neighbors) - farmer_neighbors);
		
		// Calculate the probability that the individual will learn farming based on the behavorial characteristics of its neighbors and its teacher preference.
		LP = L * (farmer_neighbors) / (farmer_neighbors + gamma * HG_neighbors);
		
		// Only ruyn if the probability of learning is greater than zero
		if ((isNAN(LP) == F) & LP != 0.0)
		{	

			// Frequency of the interaction -> Individual learns at a probability defined above
			for (i in seqLen(rpois(1, LP)))
			{
				// Change to farmer
				individual.z = 1;
				
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
			}
		}
	}
}
```

Individual HGs can learn from farmers and their z coordinate changes but their ancestry stays the same.

Learning is based on the probability that the individual will learn farming based on the behavorial characteristics of its neighbors and its teacher preference.

HGs learn at a higher probability when surrounded by more farmers

#### Competition


First we evaluate the interaction type and define vectors of the two groups

```
late()
{
	i1.evaluate(p1);
	
	// define vector of farmers and vector of HGs
	farmers = p1.individuals[p1.individuals.z == 1];
	HGs = p1.individuals[p1.individuals.z == 0];
	
	// Life table based individual mortality, get vector of individual ages
	farmer_ages = farmers.age;
	
	// Do the same for HGs if there are still HGs left
	if (length(HGs) != 0)
		HG_ages = HGs.age;
		
```

This part handles competition between nearby individuals. This is density dependent. There can be different K's for HGs and farmers here which we will see later. First we count the number nearby competing individuals. This process is slightly different depending on if the parameter is selected where individuals compete with everyone or just their own group. It counts the number of nearby individuals and calculates the density.

```
	if (C == 0) // If individuals only complete within their own groups
	{
		// Count number of neighbors within S for farmers
		farmers_num_in_s = sapply(farmers, "sum(i1.nearestNeighbors(applyValue, count = farmers.length()).z==1.0);");
		
		// Do the same for HGs if there are still HGs left
		if (length(HGs) != 0)
			HG_num_in_s = sapply(HGs, "sum(i1.nearestNeighbors(applyValue, count = HGs.length()).z==0.0);");;
		
	}
	else // Individuals complete with everyone
	{
		// Count number of neighbors within S for everyone
		farmers_num_in_s = i1.interactingNeighborCount(farmers);
		HG_num_in_s = i1.interactingNeighborCount(HGs);
		
	}
	
	// Modify mortality curve to account for population density around the indiviudal
		dens_scal_farmer = (farmers_num_in_s + 1) / (PI * (S^2) * FK + 1); // Density of farmers in the given area
		dens = mean(dens_scal_farmer);
		defineGlobal("density", dens);
		
		// Do the same if there are still HGs left
		if (length(HGs) != 0)
			dens_scal_HG = (HG_num_in_s + 1) / (PI * (S^2) * HGK + 1); // Density of HGs in the given area
```

This next part keeps individuals from living beyond realistic limits. Without this individuals in the sim can live hundreds of years because death it not dependent on age, only population density. (See life table above)

```
	scaled_mortality_farmer = ((dens_scal_farmer - 1) * scal_fac + 1) * age_scale[farmer_ages]; // Scale by a growth factor if desired. If scaling factor = 0 there isn't any scaling, population grows unrestricted based on the other parameters that govern equalibrium.
	
	if (length(HGs) != 0)
		scaled_mortality_HG = ((dens_scal_HG - 1) * scal_fac + 1) * age_scale[HG_ages]; // Scale by a growth factor if desired. If scaling factor = 0 there isn't any scaling, population grows unrestricted based on the other parameters that govern equalibrium.
	
	// Set a maximum age and make sure there are no negative fittnesses
	scaled_mortality_farmer[farmer_ages == length(age_scale) - 1] = 1;
	scaled_mortality_farmer[scaled_mortality_farmer > 1] = 1;
	
	// Do the same for HGs if there are still HGs left
	if (length(HGs) != 0)
	{
		scaled_mortality_HG[HG_ages == length(age_scale) - 1] = 1;
		scaled_mortality_HG[scaled_mortality_HG > 1] = 1;
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
				
	expected_deaths = sum(scaled_mortality_farmer);
	defineGlobal("expected_deaths", expected_deaths);
}
```


#### Movement of individuals

```
late()
{
	// ---------------------------------------------------
	//  MOVEMENT --> How individuals move around
	// ---------------------------------------------------
	for (ind in p1.individuals)
	{
		// How far individuals diffuse away from their location
		// While the new position is off the map or in the ocean, a new position will be selected until the point is in bounds
		do
		{
			coordinate_matrix = rmvnorm(1, c(0, 0), matrix(c(SDX^2, 0, 0, SDY^2), nrow=2));
			coordinates = c(coordinate_matrix[0], coordinate_matrix[1]);
			
			// Next we can reset the position
			newPos = ind.spatialPosition + coordinates;
		}
		while (!p1.pointInBounds(newPos) | p1.spatialMapValue("map_object", newPos) == 0.0);
		ind.setSpatialPosition(newPos);
	}
}
```

The individuals cannot move to locations outside of the bounds of the map. They can potentially jump across the water (simulating water travel) to other land masses, assuming it is not beyond their movement range, but they cannot stay in the ocean. Of course if you chose to run the simulation with the simple black square the individuals can move anywhere within the given map size.

#### Birth rate function

There is a small function that runs prior to the rest of the code each year that calculates the number of individuals born, and the amount of reproductively aged individuals.

```
early()
{
	// calculate the number of new births this year
	new_births = sum(p1.individuals.age == 0);
	
	// calculate the number of reproductively aged individuals
	repro_age_inds = (sum(p1.individuals.age > min_repro_age));
	
	// calculate reproduction frequency per reproductive age individual
	repro_freq = (new_births / repro_age_inds);
	
	defineGlobal("new_births", new_births);
	defineGlobal("repro_age_inds", repro_age_inds);
	defineGlobal("repro_freq", repro_freq);
}
```

#### Finally, write output to files

This first early function only runs once and generates the headers for the files before the simulation adds the values
```
1 early()
{
	// log runtime params
	param_string = paste(SN, HGK, FK, S, C, SDX, SDY, L, gamma, LD, MD, MP, min_repro_age, scal_fac, map_style, num_partitions, water_crossings, map_size_length, map_size_width, "[", age_scale, "]");
	
	// File headings
	param_heading = paste("SN HGK FK S C SDX SDY L gamma LD MD MP min_repro_age scal_fac map_style num_partitions water_crossings map_size_length map_size_width [ age_scale ]");
	
	// Runtime params - write to file
	output_runtime_file_name = ("/sim_runtime_params_" + output_name + ".txt");
	writeFile(wd + output_runtime_file_name, param_heading, append=T);
	writeFile(wd + output_runtime_file_name, param_string, append=T);
	
	// Population stats headers - write to file
	stats_header_string = paste("Year", "PopulationSize", "TotalFarmers", "TotalHGs", "RatioFarmertoHG", "Farmer_Ancestry_All", "Farmer_Ancestry_Farmers", "Farmer_Ancestry_HGs", "Num_Repro_Age_Inds", "ReproFreq", "New_Births", "Density", "expected_deaths", sep=",");
	output_stats_file_name = ("/sim_pop_stats_per_year_" + output_name + ".txt");
	writeFile(wd + output_stats_file_name, stats_header_string, append=T);
	if (map_style == 5)
	{
		for (part in 1:num_partitions)
		{
			// Outputs for partition headers
			// *************************
			if (part == 1)
			{
				farmers_per_part = "Farmers_in_Partition1";
				HGs_per_part = "HGs_in_Partition1";
				all_per_part = "All_in_Partition1";
				farmer_ancestry_farmers = "Farmer_Ancestry_Partition1_Farmers";
				farmer_ancestry_HGs = "Farmer_Ancestry_Partition1_HGs";
				farmer_ancestry_all = "Farmer_Ancestry_Partition1_All";
				ratio_per_part = "RatioFarmerToHG_Partition1";
			}
			else
			{
				farmers_per_part = paste(farmers_per_part, ",", "Farmers_in_Partition", part, sep="");
				HGs_per_part = paste(HGs_per_part, ",", "HGs_in_Partition", part, sep="");
				all_per_part = paste(all_per_part, ",", "All_in_Partition", part, sep="");
				farmer_ancestry_farmers = paste(farmer_ancestry_farmers, ",", "Farmer_Ancestry_Partition", part, "_Farmers", sep="");
				farmer_ancestry_HGs = paste(farmer_ancestry_HGs, ",", "Farmer_Ancestry_Partition", part, "_HGs", sep="");
				farmer_ancestry_all = paste(farmer_ancestry_all, ",", "Farmer_Ancestry_Partition", part, "_All", sep="");
				ratio_per_part = paste(ratio_per_part, ",", "RatioFarmerToHG_Partition", part, sep="");
			}
		}
		
		// Final Loop Output
		loop_output = paste(all_per_part, farmers_per_part, HGs_per_part, ratio_per_part, farmer_ancestry_farmers, farmer_ancestry_HGs, farmer_ancestry_all, sep=",");
		
		// Wave stats headers - write to file
		wave_stats_header_string = paste("Year", "PopulationSize", "TotalFarmers", "TotalHGs", "RatioFarmerToHG", loop_output, "Farmer_Ancestry_All", "Farmer_Ancestry_All_Farmers", "Farmer_Ancestry_All_HGs", "Num_Repro_Age_Inds", "NewBirths", "ReproFreq", sep=",");
		output_wave_stats_file_name = ("/sim_square_wave_stats_per_year_" + output_name + ".txt");
		writeFile(wd + output_wave_stats_file_name, wave_stats_header_string, append=T);
	}
}

1:6000 late()
{
	// provide feedback on progress for command line users
	year_counter = paste("Simulation Year: ", sim.cycle);
	print(year_counter);
	
	// define vector of farmers and vector of HGs
	farmers = p1.individuals[p1.individuals.z == 1];
	HGs = p1.individuals[p1.individuals.z == 0];
	all_inds = p1.individuals;
	
	// calculate num farmers
	num_farmers = sum(p1.individuals.z);
	
	//calculate num HGs
	num_HGs = (p1.individuals.length() - sum(p1.individuals.z));
	
	// calculate the ratio of farmers in the total population to file
	ratio = (sum(p1.individuals.z) / p1.individuals.length());
	
	// calculate population size statistics
	pop_size = p1.individuals.length();
	
	// calculate the overall farming ancestry
	farmer_ancestry_all = mean(all_inds.countOfMutationsOfType(m1) / (sim.chromosome.lastPosition * 2 + 2));
	
	// calculate the farming ancestry in all farmers
	if (length(farmers) != 0)
		farmer_ancestry_farmers = mean(farmers.countOfMutationsOfType(m1) / (sim.chromosome.lastPosition * 2 + 2));
	else
		farmer_ancestry_farmers = 0.0;
	
	// calculate the farming ancestry in all HGs
	if (length(HGs) != 0)
		farmer_ancestry_HGs = mean(HGs.countOfMutationsOfType(m1) / (sim.chromosome.lastPosition * 2 + 2));
	else
		farmer_ancestry_HGs = 0.0;
	
	// write outputs
	output_string = paste(sim.cycle, pop_size, num_farmers, num_HGs, ratio, farmer_ancestry_all, farmer_ancestry_farmers, farmer_ancestry_HGs, repro_age_inds, repro_freq, new_births, density, expected_deaths, sep=",");
	output_stats_file_name = ("/sim_pop_stats_per_year_" + output_name + ".txt");
	writeFile(wd + output_stats_file_name, output_string, append=T);
	
	// log individual data
	//for (ind in p1.individuals)
	//{
	//individuals = paste(sim.cycle , ind.z, ind.x, ind.y);
	//writeFile(wd + "/sim_individuals.txt", individuals, append=T);
	//}
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
		all_inds = p1.individuals;
		
		// Split width into equal parts
		partition_widths = map_size_width / num_partitions;
		for (part in 1:num_partitions)
		{
			// Collect vector of farming individuals in particular partition
			if (part == 1)
				farmers_partition = farmers[farmers.x <= partition_widths];
			else if (part == 2)
				farmers_partition = farmers[farmers.x <= 2 * partition_widths & farmers.x > partition_widths];
			else
				farmers_partition = farmers[farmers.x <= part * partition_widths & farmers.x > (part - 1) * partition_widths];
			
			// Count farmers in partition
			num_farmers_partition = size(farmers_partition);
			if (length(farmers_partition) != 0)
				farmer_ancestry_partition_farmer = mean(farmers_partition.countOfMutationsOfType(m1) / (sim.chromosome.lastPosition * 2));
			else
				farmer_ancestry_partition_farmer = 0.0;
			
			// Do the same thing for HGs
			// *************************
			// Collect vector of HG individuals in particular partition
			if (part == 1)
				HGs_partition = HGs[HGs.x <= partition_widths];
			else if (part == 2)
				HGs_partition = HGs[HGs.x <= 2 * partition_widths & HGs.x > partition_widths];
			else
				HGs_partition = HGs[HGs.x <= part * partition_widths & HGs.x > (part - 1) * partition_widths];
			
			// Count HGs in partition
			num_HGs_partition = size(HGs_partition);
			if (length(HGs_partition) != 0)
				farmer_ancestry_partition_HG = mean(HGs_partition.countOfMutationsOfType(m1) / (sim.chromosome.lastPosition * 2));
			else
				farmer_ancestry_partition_HG = 0.0;
			
			// Do the same thing for all individuals
			// *************************
			// Collect vector of HG individuals in particular partition
			if (part == 1)
				all_inds_partition = all_inds[all_inds.x <= partition_widths];
			else if (part == 2)
				all_inds_partition = all_inds[all_inds.x <= 2 * partition_widths & all_inds.x > partition_widths];
			else
				all_inds_partition = all_inds[all_inds.x <= part * partition_widths & all_inds.x > (part - 1) * partition_widths];
			
			// Count all_inds in partition
			num_all_inds_partition = size(all_inds_partition);
			if (length(all_inds_partition) != 0)
				farmer_ancestry_partition_all = mean(all_inds_partition.countOfMutationsOfType(m1) / (sim.chromosome.lastPosition * 2));
			else
				farmer_ancestry_partition_all = 0.0;
			
			// calculate the ratio of farmers to HGs each partition
			// *************************
			ratio_part = (num_farmers_partition / num_all_inds_partition);
			
			// Outputs for partitions
			// *************************
			// Num individuals in each partition
			if (part == 1)
			{
				farmers_per_part = paste(num_farmers_partition);
				HGs_per_part = paste(num_HGs_partition);
				all_per_part = paste(num_all_inds_partition);
			}
			else
			{
				farmers_per_part = paste(farmers_per_part, num_farmers_partition, sep=",");
				HGs_per_part = paste(HGs_per_part, num_HGs_partition, sep=",");
				all_per_part = paste(all_per_part, num_all_inds_partition, sep=",");
			}
			
			// Ratio F:HGs in each partition
			if (part == 1)
				ratio_per_part = paste(ratio_part);
			else
				ratio_per_part = paste(ratio_per_part, ratio_part, sep=",");
			
			// Farmer ancestry:
			if (part == 1)
			{
				// In Farmers
				farmer_ancestry_farmers = paste(farmer_ancestry_partition_farmer);
				
				// In HGs
				farmer_ancestry_HGs = paste(farmer_ancestry_partition_HG);
				
				// In All
				farmer_ancestry_all = paste(farmer_ancestry_partition_all);
			}
			else
			{
				// In Farmers
				farmer_ancestry_farmers = paste(farmer_ancestry_farmers, farmer_ancestry_partition_farmer, sep=",");
				
				// In HGs
				farmer_ancestry_HGs = paste(farmer_ancestry_HGs, farmer_ancestry_partition_HG, sep=",");
				
				// In All
				farmer_ancestry_all = paste(farmer_ancestry_all, farmer_ancestry_partition_all, sep=",");
			}
		}
		
		// Final Loop Output
		loop_output = paste(all_per_part, farmers_per_part, HGs_per_part, ratio_per_part, farmer_ancestry_farmers, farmer_ancestry_HGs, farmer_ancestry_all, sep=",");
		
		// calculate total num farmers
		num_farmers = sum(p1.individuals.z);
		
		// calculate total num HGs
		num_HGs = (p1.individuals.length() - sum(p1.individuals.z));
		
		// calculate population size statistics
		pop_size = p1.individuals.length();
		
		// calculate the ratio of farmers in the total population to file
		ratio = (sum(p1.individuals.z) / p1.individuals.length());
		
		// calculate the overall farming ancestry
		farmer_ancestry_all = mean(all_inds.countOfMutationsOfType(m1) / (sim.chromosome.lastPosition * 2 + 2));
		
		/// calculate the farming ancestry in all farmers
		if (length(farmers) != 0)
			farmer_ancestry_farmers = mean(farmers.countOfMutationsOfType(m1) / (sim.chromosome.lastPosition * 2 + 2));
		else
			farmer_ancestry_farmers = 0.0;
		
		// calculate the farming ancestry in all HGs
		if (length(HGs) != 0)
			farmer_ancestry_HGs = mean(HGs.countOfMutationsOfType(m1) / (sim.chromosome.lastPosition * 2 + 2));
		else
			farmer_ancestry_HGs = 0.0;
		
		// write outputs
		output_string = paste(sim.cycle, pop_size, num_farmers, num_HGs, ratio, loop_output, farmer_ancestry_all, farmer_ancestry_farmers, farmer_ancestry_HGs, repro_age_inds, new_births, repro_freq, sep=",");
		
		// output to file
		output_stats_file_name = ("/sim_square_wave_stats_per_year_" + output_name + ".txt");
		writeFile(wd + output_stats_file_name, output_string, append=T);
	}
}
```
