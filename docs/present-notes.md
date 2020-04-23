# Note for video presentation

## Background/motivation

The motivation behind this project is that global climate change, and in particular sea level rise, represents a growing threat, and is therefore a key target for scientific research.
The West Antarctic Ice Sheet, due to it's large size and relative instability, will likely dominate sea level rise in the future.
The SMB, or total snow accumulation in a year, is an important influence on the ice sheet sea level contribution, and is therefore the focus of our project, where we investigate the spatial and temporal correlations of annual snow accumulation with various climatic and topographic variables over the central West Antarctic Ice Sheet.

## Data/methods

To do this, we principally rely on 3 data sets.
The first dataset consists of annual accumulation estimates, consisting of accumulation time series spaced every 25 meters along a radar flightline.
You can see this flightline here in the greater context of Antarctica.
These time series serve as our primary dependent variable in our analysis.

Our second data set consists of surface topography attributes across Antarctica at 8-m resolution, which we use as predictors in our model.
The final dataset consists of various climatic variables with 30x30 km grid cells, of which we are most interested in near surface wind vectors to use as additional predictors of accumulation.

A significant portion of our project consisted of combining all these data into a consistent and tidy format.
Getting all our data into the same coordinate reference system, resolution, temporal span, and geographic projection, and extracting the data only at our points of interest, was non-trivial.
This is particularly challenging in Antarctica, where the proximity to the South Pole causes our traditional notions of cardinal directions to break down.

After cleaning and formatting, we used weighted least-squares regression to explore temporal changes in accumulation rates, and multilinear regression to quantify the spatial controls on mean accumulation rate and on the temporal trend in accumulation rates.

So taking a look at some of the key results from this, we start with...

## Results/visualizations/conclusions

This figure shows the mean accumulation rate over the period 1979-2008 in mm w.e.
You can see this plot zooms into our region of interest, with this same triangular flightline present.
The key observation for this figure is the low accumulation rates to the southwest (remember directions are weird in Antarctica), with increasing accumulation rates moving north and east.

The next results figure shows the change in accumulation rate from 1979-2008, expressed as a percent change per year from the mean.
The main observation here is the rough East-West gradient in trends, with strong decreases in accumulation rates of more than 2-3% per year to the west and relatively stable rates as you move eastward.

Finally we have a map of wind speed (on the left) and dominant wind direction (on the right).
These show generally higher wind speeds to the north and southwest.
Wind direction for the most part is consistently east by southeast, but a significant change occurs in the southwesternmost region, where then dominant wind direction transitions to north by northeast.

The multilinear spatial regression led to some interesting results worth mentioning.
All of our predictors were significant in predicting mean accumulation rate and temporal trends in accumulation, but the two most influential were the wind direction and the northing component.
We predict regions more towards the interior to have lower accumulation rates, as well as regions with more northerly wind flow.
Together our predictors can explain ~78% of the spatial variability in mean accumulation rates.

The most important controlling spatial variables on the temporal accumulation trend are again the wind direction and northing location, with the mean accumulation rate also influential.
The modeling results suggest anomalously dry western locations with more southerly wind flow are the most susceptible to accumulation losses.
Our selected predictors explain ~60% of the spatial variability in the temporal trends in accumulation in the region.
