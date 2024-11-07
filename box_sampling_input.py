## STEP 1: define box extent
# lateral dimensions/width (in longitudinal direction) and height (in latitudinal direction)
width = 170 # box width (longitude), must be an integer less than 180
height = 110 # box height (latitude), must be an integer

# south (latitude) and west (longitude) boundary values of the box
south = -70 # southern box extent (latitude), must be an integer
west = 120 # western box extent (longitude), must be an integer

## define box dimensions:
# radial dimensions/box top and bottom:
box_depth_min = 2891 - 1000
box_depth_max = 2891



## STEP 2: define input dataset and its attributes
# define file names, columns to keep from the input dataset
dataset_input_name = 'hongyu_data.csv'
columns_to_keep = ['COMPREHENSIVE_WEIGHT']
dataset_save_name = 'dataset_with_path_length.csv'