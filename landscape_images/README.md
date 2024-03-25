*Preferences regarding the landscape map on which the individuals move*

```
	
	// Map Prefs: 
	// ***********************************	
	if (!exists("map_style"))
		defineConstant("map_style", 5); // Parameter of 0 = no topography, 1 = super light topography, 2 = light topography, 3 = regular topography, 4 heavy topography, 5 = square, 6 = custom map
	if (!exists("num_bins"))
		defineConstant("num_bins", 20); // Number of bins if using the square
	if (!exists("water_crossings"))
		defineConstant("water_crossings", 1); //Parameter of 0 = no water crossing paths, 1 = water crossing paths
	if (water_crossings == 1)
		defineConstant("file_extention", "_water.png");
	else
		defineConstant("file_extention", ".png");
```

"map_style" dictates the landscape type for the simulation. 
Parameter of 0 = no topography, 1 = super light topography, 2 = light topography, 3 = regular topography, 4 = heavy topography, 5 = square, 6 = custom map

"num_bins" is a parameter used only when running the sim on map style 5 (square). This parameter divides the landscape into discrete bins and allows for measurement of expansion speed and tracking simulation progress when run in a simple landscape model. This parameter controls the number of bins.

"water_crossings" provides the option that individuals have faux "land bridges" that simulate realistic routes that individuals traveled via waterways. This parameter only applies when run on the landscape map. The binary nature of our landscape model means that water travel may not be possible otherwise depending on moment range. This parameter provides the option to allow incremental movement over water.
 

*The following parameters below **do not need to be changed**! They handle what map to use. If you want to use a custom map, enter the file name at the beginning of the script*

```
	defineConstant("mapfile_none", wd + "/EEA_map" + file_extention); // File Path to Map Image
	defineConstant("mapfile_topo_superlight", wd + "/EEA_map_topo_superlight" + file_extention); // File Path to Map Image
	defineConstant("mapfile_topo_light", wd + "/EEA_map_topo_light" + file_extention); // File Path to Map Image
	defineConstant("mapfile_topo_regular", wd + "/EEA_map_topo_regular" + file_extention); // File Path to Map Image
	defineConstant("mapfile_topo_heavy", wd + "/EEA_map_topo_heavy" + file_extention); // File Path to Map Image
	defineConstant("square", wd + "/square.png"); // File Path to Map Image
	defineConstant("mapfile_island_removed", wd + "/EEA_map_island_removed" + file_extention); // File Path to Map Image
	defineConstant("custom_map", wd + custom_map_filename); // File Path to Map Image
```
***The maps are available for download in the repo.** The path points to where the map is in your system. If you altered the name of the map the name of the map MUST be changed to fit your file structure. It is recommended that you DO NOT change the map names. If you do not alter them, setting the WD at the top of the script is sufficient.*

*The next parameter set is the length and width of the map. This will be used when building the map.*

The following parameters handle the map size in km. If using provided EEA maps from the repo, length and width for map style images should be kept at a 1 to 1 aspect ratio (square) to avoid distortion. By default they are 3700 x 3700 km square maps

```
	if (!exists("map_size_length"))
		defineConstant("map_size_length", 3700);
	if (!exists("map_size_width"))
		defineConstant("map_size_width", 3700);;
	
	// Map citation: https://www.eea.europa.eu/ds_resolveuid/558D91E1-3DB0-4639-9F70-2012CC4453A5
```
