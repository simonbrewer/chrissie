# chrissie

Test of hydra model for Lake Chrissie region

1. Run extractCLclim.R in cru/cl - crops out an approximate region for cliamte data
2. Run getClimDailyGrid.R - calculates climate indicies for this region
3. Run regGrid.R - regrids the climate indicies to the finer DEM grid
4. Run testHydra.R