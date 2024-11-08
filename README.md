# EMPRESS (Earth Modeling PRioritized Event Selection Scheme)
---


## Overview:
EMPRESS is a workflow designed to optimize event selection for regional seismic modeling, particularly for full-waveform inversions that are limited by the number of earthquakes that can be modeled. EMPRESS includes several interactive filtering and weighting steps to optimize various event attributes and contributions, including phase measurement quality, total number of raypaths, regional sampling, and azimuthal coverage. The workflow is executed in a Jupyter notebook (`EMPRESS.ipynb`) and includes visualization tools to show how various filters/weights affect overall coverage of a region of interest.


## EMPRESS package contents:

1. `EMPRESS.ipynb` - the Jupyter notebook containing the interactive EMPRESS workflow
2. `box_sampling.py` - a Python script that calculates the total box sampling of each raypath in an input file
3. `box_sampling_input.py` - an input file with variables to define the input dataset file and attributes, as well as the modeling region.
4. `geo_math.py` and `refmodels.py` - libraries of necessary functions for `EMPRESS.ipynb` and `box_sampling.py`
5. `/EMPRESS/events_lists` - a subdirectory where optimized lists of events will be stored


## Usage notes:
### **Step 1**:
Define the lateral (latitude/longitude) and radial (depth) extents of the region of interest, referred to as the "box". For this step, the box extent is first defined in `box_sampling_input.py` and then `EMPRESS.ipynb` is used to vizualize the resulting box. The extent variables can be iteratively modified and visualized until an appropriate region is selected. Note that you must use the -180/180 longitude convention and the -90/90 latitude convention. To do this, define the following variables in `box_sampling_input.py`:

- `width` - longitudinal extent of the box in degrees
- `height` - latitudinal extent of the box in degrees
- `south` - latitude that defines the southernmost extent of the box in degrees 
- `west` - longitude that defines the westernmost extent of the box in degrees 
- `box_depth_min` - depth that defines the top of the box in km
- `box_depth_max` - depth that defines the bottom/base of the box in km

and then run `PART 1` of `EMPRESS.ipynb`. If the box doesn't cover the target region appropriately, update `box_sampling_input.py` and rerun `PART 1`. Do this until a satisfactory box has been defined.


### **Step 2**:
Run `box_sampling.py` to compute the ray theoretical percentage of box sampling for each phase measurement in the input dataset, typically a dataset of travel time measurements. The input dataset must be a CSV of event-station pairs, defined with the following headers:

- `EVENT_ID` - a unique ID for each event, typically a datetime for the earthquake origin time
- `EQ_LAT` - latitude of the event
- `EQ_LON` - longitude of the event
- `EQ_MAG` - magnitude of the event
- `EQ_DEP` - depth of the event
- `STA_LAT` - latitude of the station
- `STA_LON` - longitude of the station
- `PHASE` - phase measured

Note that the naming convention for phases in `PHASE` must be the same as in the TauP toolkit (Crotwell et al., 1999; https://github.com/crotwell/TauP). That is, phase names must not include integers; indicate multiples as longer strings rather than with numbers. For example, `S3` should be instead listed as `SSS`. Use `m` at the end of the entire phase name to indicate a major arc path (e.g., `ScSScSScSm`) and/or `s` at the beginning of the entire phase name to indicate a depth phase (e.g., `sS`).

Once the CSV is properly formatted and uploaded to the `/EMPRESS` directory, define the following in `box_sampling_input.py`:

- `dataset_input_name` - full name of the input dataset file (e.g., `input_file.csv`)
- `columns_to_keep` - any additional headers to be retained to be used as filters/weights later (e.g., information about measurement quality). These must be numerical values.
- `dataset_save_name` - full name of the output dataset file, which will include the ray theoretical percentage of box sampling for each phase measurement (e.g., `output_file.csv`)

The output dataset file will contain the initial information from the input dataset as well as:

- `AZIMUTH` - the azimuth of each raypath between the event and the station
- `MAX_DEPTH_KM` - the bottoming depth of each raypath between the event and the station
- `TOTAL_LENGTH_KM` - the total length in km of each raypath between the event and the station
- `LENGTH_IN_BOX_KM` - the length of the portion of each raypath that samples the box
- `%_IN_BOX` - the ray theoretical percentage of the portion of each raypath that samples the box

In addition to the output dataset file, `box_sampling.py` also creates the subdirectory `/EMPRESS/event_files/`. This directory is populated with CSVs for each event that contain the entry and exit coordinates of each raypath segment that samples the box.

Note that `box_sampling.py` can take several hours or more to run depending on the computer and the size of the input dataset.


### **Step 3**:
Use `EMPRESS.ipynb` to interactively apply filters and weights to the dataset and develop an optimized list of events. The user may indicate whether or not to save the filters/weights that they have applied to come up with a given list of events. `EMPRESS.ipynb` output includes:

- `/EMPRESS/event_lists/{list_name}_list.csv` - the actual list of events, in order of priority
- `/EMPRESS/event_lists/{list_name}_params.csv` - the filters/weights applied to develop the corresponding list of events




