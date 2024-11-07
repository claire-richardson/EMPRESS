import os
import shutil
import geo_math
import numpy as np
import pandas as pd
import box_sampling_input
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from matplotlib.dates import date2num
ref_model = TauPyModel(model = 'prem')

## DEFINE FUNCTIONS
def find_nearest(array, value, type):
    if type == 'absolute':
        array_diff = abs(array - value)
        id = np.where(array_diff == array_diff.min())[0][0]
        return [id]
    
    elif type == 'bounds':
        array_diff = array - value
        smallest_negative_value = array_diff[np.where(array_diff <= 0)].max()
        smallest_positive_value = array_diff[np.where(array_diff >= 0)].min()
        negative_id = np.where(array_diff == smallest_negative_value)[0][0]
        positive_id = np.where(array_diff == smallest_positive_value)[0][0]
        return [negative_id, positive_id]

def dist_at_depth(z1, dist1, z2, dist2, z_mid):
    R = (z_mid - z1) / (z2 - z1)
    return round(dist1 + ((dist2 - dist1) * R), 10)

def vector_mag(vector):
    return np.sqrt((vector[0] ** 2) + (vector[1] ** 2) + (vector[2] ** 2))

def vector_norm(vector):
    mag = vector_mag(vector)
    x = vector[0] / mag
    y = vector[1] / mag
    z = vector[2] / mag
    return [x, y, z]

def check_lon_extent(intersection_point, starting_point, array):
    if starting_point == 'eq':
        intersection_dist = p1_dist + intersection_point[2]
    elif starting_point == 'sta':
        intersection_dist = p2_dist - intersection_point[2]
    
    nearest_int_ids = find_nearest(array[:, 2], intersection_dist, 'bounds')                                    
    int_dist_1 = array[nearest_int_ids[0], 2]
    int_dist_2 = array[nearest_int_ids[1], 2]
    int_depth_1 = array[nearest_int_ids[0], 3]
    int_depth_2 = array[nearest_int_ids[1], 3]
    
    intersection_depth = geo_math.new_depth(int_depth_1, int_depth_2, int_dist_1, int_dist_2, intersection_dist)

    intersection_lats.append(intersection_point[0])
    intersection_lons.append(intersection_point[1])
    intersection_dists.append(intersection_dist)
    intersection_depths.append(intersection_depth)

def find_gc_intersections(start_vector, end_vector, box_lat1, box_lon1, box_lat2, box_lon2):
    intersect_coords = []
    a21 = geo_math.coord2cart(box_lat1, box_lon1)
    a22 = geo_math.coord2cart(box_lat2, box_lon2)
    N1 = np.cross(a11, a12)
    N2 = np.cross(a21, a22)
    L = np.cross(N1, N2)
    I1 = vector_norm(L)
    I2 = [-I1[0], -I1[1], -I1[2]]
    I1_pt = geo_math.cart2coord(I1[0], I1[1], I1[2])
    I2_pt = geo_math.cart2coord(I2[0], I2[1], I2[2])

    # first, test if I1 lies on the first arc (the raypath)
    theta_a11_i1 = np.degrees(np.arccos(np.dot(start_vector, I1) / (vector_mag(start_vector) * vector_mag(I1))))
    theta_a12_i1 = np.degrees(np.arccos(np.dot(end_vector, I1) / (vector_mag(end_vector) * vector_mag(I1))))
    theta_a13_i1 = np.degrees(np.arccos(np.dot(start_vector, end_vector) / (vector_mag(start_vector) * vector_mag(end_vector))))
    
    # if it does, check if it lies on the second arc as well
    a1_i1 = round(theta_a11_i1 + theta_a12_i1, 1)
    if a1_i1 == round(theta_a13_i1, 1):
        theta_a21_i1 = np.degrees(np.arccos(np.dot(a21, I1) / (vector_mag(a21) * vector_mag(I1))))
        theta_a22_i1 = np.degrees(np.arccos(np.dot(a22, I1) / (vector_mag(a22) * vector_mag(I1))))
        theta_a23_i1 = np.degrees(np.arccos(np.dot(a21, a22) / (vector_mag(a21) * vector_mag(a22))))

        # if the first intersection point (I1) also lies on the second arc, it's an intersection point
        a2_i1 = round(theta_a21_i1 + theta_a22_i1, 1)
        if a2_i1 == round(theta_a23_i1, 1):
            intersect_coords.append([round(I1_pt[0], 5), round(I1_pt[1], 5)])

    # now, repeat the tests for I2
    theta_a11_i2 = np.degrees(np.arccos(np.dot(start_vector, I2) / (vector_mag(start_vector) * vector_mag(I2))))
    theta_a12_i2 = np.degrees(np.arccos(np.dot(end_vector, I2) / (vector_mag(end_vector) * vector_mag(I2))))
    theta_a13_i2 = np.degrees(np.arccos(np.dot(start_vector, end_vector) / (vector_mag(start_vector) * vector_mag(end_vector))))
    
    # if I2 falls on the first arc, check if it lies on the second arc as well
    a1_i2 = round(theta_a11_i2 + theta_a12_i2, 1)
    if a1_i2 == round(theta_a13_i2, 1):
        theta_a21_i2 = np.degrees(np.arccos(np.dot(a21, I2) / (vector_mag(a21) * vector_mag(I2))))
        theta_a22_i2 = np.degrees(np.arccos(np.dot(a22, I2) / (vector_mag(a22) * vector_mag(I2))))
        theta_a23_i2 = np.degrees(np.arccos(np.dot(a21, a22) / (vector_mag(a21) * vector_mag(a22))))

        # if the second intersection point (I2) also lies on the second arc, it's also an intersection point
        a2_i2 = round(theta_a21_i2 + theta_a22_i2, 1)
        if a2_i2 == round(theta_a23_i2, 1):
            intersect_coords.append([round(I2_pt[0], 5), round(I2_pt[1], 5)])

    return intersect_coords

def taup_raypath(eq_depth, eq_lat, eq_lon, sta_lat, sta_lon, phase):
    df_pathfile = pd.DataFrame()
    phase_name = phase.split('_')[0]
    if 'm' in phase_name:
        phase2 = phase_name.replace('m', '')
        arrivals = ref_model.get_ray_paths_geo(eq_depth, eq_lat, eq_lon, sta_lat, sta_lon, [phase2])
        if not arrivals:
            pass
        else:
            try:
                for arrival in arrivals:
                    arr_dist = arrival.purist_distance
                    if arr_dist > 180.:
                        pathfile = arrival.path
                        df_pathfile = pd.DataFrame(pathfile)
                        df_pathfile['dist'] = df_pathfile['dist'] * (180. / np.pi)
                        df_pathfile = df_pathfile[['lat', 'lon', 'dist', 'depth']]
                        break
            except:
                pass
    else:
        arrivals = ref_model.get_ray_paths_geo(eq_depth, eq_lat, eq_lon, sta_lat, sta_lon, [phase_name])
        if not arrivals:
            pass
        else:            
            arrival = arrivals[0]
            pathfile = arrival.path
            df_pathfile = pd.DataFrame(pathfile)
            df_pathfile['dist'] = df_pathfile['dist'] * (180. / np.pi)
            df_pathfile = df_pathfile[['lat', 'lon', 'dist', 'depth']]
    
    return np.array(df_pathfile)

def raypath_length(raypath_array):
    arr_dist = raypath_array[1:, 2] - raypath_array[0:-1, 2] # epicentral distance of each segment
    arr_dist = np.insert(arr_dist, 0, 0.)
    arr_dist = np.radians(arr_dist)
    arr_radius = 6371 - raypath_array[:, 3] # radius of each point    
    arr_length = np.sqrt(np.power(arr_radius[:-1], 2) + np.power(arr_radius[1:], 2) - (2 * arr_radius[:-1] * arr_radius[1:] * np.cos(arr_dist[1:])))
    arr_length = np.insert(arr_length, 0, 0.)
    mantle_ids = np.where(raypath_array[:,3] == 2891)[0]
    if (len(mantle_ids) == 2) and (mantle_ids[1] - mantle_ids[0] == 1):
        d1 = raypath_array[mantle_ids[0], 2]
        d2 = raypath_array[mantle_ids[1], 2]
        length = (2 * np.pi * (6371 - 2891)) / 360 * (d2 - d1)
        arr_length[mantle_ids[1]] = length
    else:
        pass
    return arr_length

## START CODE:
try:
    os.mkdir(f'./event_files')
except:
    shutil.rmtree(f'./event_files')
    os.mkdir(f'./event_files')

cols = ['EVENT_ID', 'EQ_LAT', 'EQ_LON', 'EQ_MAG', 'EQ_DEP', 'STA_LAT', 'STA_LON', 'PHASE']

added_col_ids = []
for col_name in box_sampling_input.columns_to_keep:
    added_col_ids.append(len(cols))
    cols.append(col_name)
df_data = pd.read_csv(box_sampling_input.dataset_input_name)
df_data = df_data[cols]
df_data['AZIMUTH'] = 0.
# df_data['MEAN_AZ_IN_BOX'] = 0.
# df_data['STD_AZ_IN_BOX'] = 0.
df_data['MAX_DEPTH_KM'] = 0.
df_data['TOTAL_LENGTH_KM'] = 0.
df_data['LENGTH_IN_BOX_KM'] = 0.
df_data['%_IN_BOX'] = 0.

w = box_sampling_input.west
s = box_sampling_input.south
e = w + box_sampling_input.width
if e >= 180.:
    e -= 360

n = s + box_sampling_input.height
if n >= 90.:
    raise Exception(f'--  Northern extent too large, n = {n}; adjust `height` variable in `box_sampling_input.py`')

box_depth_min = box_sampling_input.box_depth_min
box_depth_max = box_sampling_input.box_depth_max

df_data_ar = np.array(df_data)
idx_to_drop = []

# for every line in the dataset:
for line in range(len(df_data)):
    try:
        eq_id = df_data_ar[line, 0]
        eq_lat = df_data_ar[line, 1]
        eq_lon = df_data_ar[line, 2]
        eq_dep = df_data_ar[line, 4]
        sta_lat = df_data_ar[line, 5]
        sta_lon = df_data_ar[line, 6]
        phase = df_data_ar[line, 7]
        length_in_box = 0.
        # az_in_box = []

        path_intersections = np.array([[eq_lat, eq_lon, 0., eq_dep]])
        path_az = geo_math.azimuth(eq_lat, eq_lon, sta_lat, sta_lon)

        if os.path.exists(f'./event_files/{eq_id}.csv'):
            pass
        else:
            with open(f'./event_files/{eq_id}.csv', 'w') as fout:
                fout.write('STA_LAT,STA_LON,PHASE,ENTRY_LAT,ENTRY_LON,ENTRY_DIST,ENTRY_DEP,EXIT_LAT,EXIT_LON,EXIT_DIST,EXIT_DEP,SEG_LENGTH\n')

        ## FIRST, FIND BOTTOMING DEPTH
        arr_raypath = taup_raypath(eq_dep, eq_lat, eq_lon, sta_lat, sta_lon, phase)
        if len(arr_raypath) == 0:
            idx_to_drop.append(line)
            pass
        else:
            max_depth = arr_raypath.T[3].max()
            total_length = raypath_length(arr_raypath).sum()
        
            # if the path bottoms above the box, the whole path can just be omitted.
            if max_depth < box_depth_min:
                pass
            
            # otherwise, loop through the subsegments that dip into the box (this will be only one segment for paths without reflections)
            else:
                # find all of the indices where depth == 0. this will split the path up into segments if it has reflections.
                zero_depth_idx = np.where(arr_raypath.T[3] == 0)[0]
                
                # then loop through each of those segments to narrow them even further to just the bit that falls in the depth range.
                idx_start = 0
                for zd_idx in zero_depth_idx:
                    # find the indices of the points of the current leg that dip into the depth range
                    arr_sub_raypath = arr_raypath[idx_start:zd_idx] # array of the entire current reflection segment
                    top_depth_idx = np.where(arr_sub_raypath.T[3] > box_depth_min)[0] # indices of the current reflection segment that are deeper than the box top
                    bottom_depth_idx = np.where(arr_sub_raypath.T[3] > box_depth_max)[0] # indices of the current reflection segment that are deeper than the box bottom                
                    
                    # if the leg doesn't actually dip into the depth range (i.e., depth phases), pass
                    if len(top_depth_idx) == 0:
                        pairs = 0
    
                    # if the leg only dips below the top of the box but stays above its base
                    elif (len(bottom_depth_idx) == 0) and (len(top_depth_idx) > 0):
                        pairs = 1
                        arr_sub_raypath_segment = arr_sub_raypath[top_depth_idx[0] - 1: top_depth_idx[-1] + 2]
                        entry_depth = box_depth_min
                        exit_depth = box_depth_min
                        
                    # if the leg dips below the top of the box AND the base of the box
                    elif (len(bottom_depth_idx) > 0) and (len(top_depth_idx) > 0):
                        pairs = 2
                    
                    for pair in range(pairs):
                        pair += 1
                        if (pair == 1) and (pairs == 2):
                            arr_sub_raypath_segment = arr_sub_raypath[top_depth_idx[0] - 1:bottom_depth_idx[0] + 1]
                            entry_depth = box_depth_min
                            exit_depth = box_depth_max
                            
                        elif pair == 2:
                            arr_sub_raypath_segment = arr_sub_raypath[bottom_depth_idx[-1]:top_depth_idx[-1] + 2]
                            entry_depth = box_depth_max
                            exit_depth = box_depth_min
    
                        # find the entry point
                        # slice to only include the portion of the path that dips into the depth range
                        # in the case where it's the first reflection segment and the earthquake is already in the box:
                        if (eq_dep > box_depth_min) and (zd_idx == zero_depth_idx[0]) and (pair == 1):
                            entry_lat = eq_lat
                            entry_lon = eq_lon
                            entry_dist = 0.
                            entry_depth = eq_dep
    
                        # otherwise, the entry point is later along the path
                        else:                        
                            # add new points for the exact top of the depth range
                            lat1 = arr_sub_raypath_segment[0, 0]
                            lon1 = arr_sub_raypath_segment[0, 1]
                            dist1 = arr_sub_raypath_segment[0, 2]
                            depth1 = arr_sub_raypath_segment[0, 3]
    
                            lat2 = arr_sub_raypath_segment[1, 0]
                            lon2 = arr_sub_raypath_segment[1, 1]
                            dist2 = arr_sub_raypath_segment[1, 2]
                            depth2 = arr_sub_raypath_segment[1, 3]
                
                            new_dist = dist_at_depth(depth1, dist1, depth2, dist2, entry_depth)
                            new_coords = geo_math.GCP_point(lat1, lon1, lat2, lon2, (dist2 - dist1), (new_dist - dist1))
    
                            entry_lat = new_coords[0]
                            entry_lon = new_coords[1]
                            entry_dist = new_dist
    
                        arr_sub_raypath_segment[0, 0] = entry_lat
                        arr_sub_raypath_segment[0, 1] = entry_lon
                        arr_sub_raypath_segment[0, 2] = entry_dist
                        arr_sub_raypath_segment[0, 3] = entry_depth
                        
                        # now, find the exit point
                        # if the exit point is actually the station, e.g., if the box covers the whole-mantle
                        if (zd_idx == zero_depth_idx[-1]) and (box_depth_min == 0.) and (pair == pairs):
                            exit_lat = sta_lat
                            exit_lon = sta_lon
                            exit_dist = arr_raypath[-1, 2]
    
                        else:        
                            lat1 = arr_sub_raypath_segment[-2, 0]
                            lon1 = arr_sub_raypath_segment[-2, 1]
                            dist1 = arr_sub_raypath_segment[-2, 2]
                            depth1 = arr_sub_raypath_segment[-2, 3]
                        
                            lat2 = arr_sub_raypath_segment[-1, 0]
                            lon2 = arr_sub_raypath_segment[-1, 1]
                            dist2 = arr_sub_raypath_segment[-1, 2]
                            depth2 = arr_sub_raypath_segment[-1, 3]
                
                            new_dist = dist_at_depth(depth1, dist1, depth2, dist2, exit_depth)
                            new_coords = geo_math.GCP_point(lat1, lon1, lat2, lon2, (dist2 - dist1), (new_dist - dist1))
    
                            exit_lat = new_coords[0]
                            exit_lon = new_coords[1]
                            exit_dist = new_dist
    
                        arr_sub_raypath_segment[-1, 0] = exit_lat
                        arr_sub_raypath_segment[-1, 1] = exit_lon
                        arr_sub_raypath_segment[-1, 2] = exit_dist
                        arr_sub_raypath_segment[-1, 3] = exit_depth
    
                        pair_dist = exit_dist - entry_dist
                        
                        if pair_dist < 180:
                            lats = [entry_lat, exit_lat]
                            lons = [entry_lon, exit_lon]
                            dists = [entry_dist, exit_dist]
                            depths = [entry_depth, exit_depth]
                    
                        elif pair_dist > 180:
                            nearest_dist_id = find_nearest(arr_sub_raypath_segment[:, 2], (entry_dist + (pair_dist / 2)), 'absolute')
                            lats = [entry_lat, arr_sub_raypath_segment[nearest_dist_id, 0], exit_lat]
                            lons = [entry_lon, arr_sub_raypath_segment[nearest_dist_id, 1], exit_lon]
                            dists = [entry_dist, arr_sub_raypath_segment[nearest_dist_id, 2], exit_dist]
                            depths = [entry_depth, arr_sub_raypath_segment[nearest_dist_id, 3], exit_depth]
    
                        intersection_lats = lats
                        intersection_lons = lons
                        intersection_dists =  dists
                        intersection_depths = depths
    
                        no_intersection_len = len(intersection_lats)
    
                        ## NOW MOVE ON TO LATERAL BOUNDARIES
                        ## FIRST, FIND THE NORTH/SOUTH BOUNDARIES:
                        for coord in range(len(lats) - 1):
                            p1_lat = lats[coord]
                            p1_lon = lons[coord]
                            p1_dist = dists[coord]
                            p1_depth = depths[coord]
                            p2_lat = lats[coord + 1]
                            p2_lon = lons[coord + 1]
                            p2_dist = dists[coord + 1]
                            p2_depth = depths[coord + 1]
                            az = geo_math.azimuth(p1_lat, p1_lon, p2_lat, p2_lon)
                            baz = geo_math.azimuth(p2_lat, p2_lon, p1_lat, p1_lon)
    
                            inflection = geo_math.inflection_finder(az, baz, p1_lat, p1_lon, p2_lat, p2_lon, p1_dist, p2_dist)
                            # if there's no inflection, there's at most one intersecting point on either or both the north and/or south boundaries.
                            if inflection[0] == False:
                                # if the entry and exit points are on either side of the north bound, find the intersection point
                                if ((p1_lat > n) and (p2_lat < n)) or ((p1_lat < n) and (p2_lat > n)):
                                    int = geo_math.known_lat(az, p1_lat, p1_lon, n)
                                    # check if the intersection point is between the east/west bounds and keep it if it is:
                                    if (int[1] < e) or (int[1] > w):
                                        check_lon_extent(int, 'eq', arr_sub_raypath_segment)
    
                                # if the entry and exit points are on either side of the south bound, find the intersection point
                                if ((p1_lat < s) and (p2_lat > s)) or ((p1_lat > s) and (p2_lat < s)):
                                    int = geo_math.known_lat(az, p1_lat, p1_lon, s)
                                    # check if the intersection point is between the east/west bounds and keep it if it is:
                                    if (int[1] < e) or (int[1] > w):
                                        check_lon_extent(int, 'eq', arr_sub_raypath_segment)
    
                            # if there's an inflection, there are potentially two intersecting points on one of the boundaries, and potentially one on the other:
                            elif inflection[0] == True:
                                # if the inflection point is above the northern boundary and the earthquake is below the northern boundary, find the intersection point
                                if (inflection[2] > n) and (p1_lat < n):
                                    int = geo_math.known_lat(az, p1_lat, p1_lon, n)
                                    # check if the intersection point is between the east/west bounds and keep it if it is:
                                    check_lon_extent(int, 'eq', arr_sub_raypath_segment)
                        
                                # if the inflection point is above the northern boundary and the station is below the northern boundary, find the intersection point
                                if (inflection[2] > n) and (p2_lat < n):
                                    int = geo_math.known_lat(baz, p2_lat, p2_lon, n)
                                    # check if the intersection point is between the east/west bounds and keep it if it is:
                                    check_lon_extent(int, 'sta', arr_sub_raypath_segment)
                                
                                # if the inflection point is in the N. hemisphere and the earthquake or station is below the southern boundary, find the intersection point
                                if ((inflection[2] > 0) and (p2_lat < s)) or ((inflection[2] > 0.) and (p1_lat < s)):
                                    int = geo_math.known_lat(az, p1_lat, p1_lon, s)
                                    # check if the intersection point is between the east/west bounds and keep it if it is:
                                    check_lon_extent(int, 'eq', arr_sub_raypath_segment)
                                
                                # if the inflection point is below the southern boundary and the earthquake is above the southern boundary, find the intersection point
                                if (inflection[2] < s) and (p1_lat > s):
                                    int = geo_math.known_lat(az, p1_lat, p1_lon, s)
                                    # check if the intersection point is between the east/west bounds and keep it if it is:
                                    check_lon_extent(int, 'eq', arr_sub_raypath_segment)
                    
                                # if the inflection point is below the southern boundary and the station is above the southern boundary, find the intersection point
                                if (inflection[2] < s) and (p2_lat > s):
                                    int = geo_math.known_lat(baz, p2_lat, p2_lon, s)
                                    # check if the intersection point is between the east/west bounds and keep it if it is:
                                    check_lon_extent(int, 'sta', arr_sub_raypath_segment)
                    
                                # if the inflection point is in the S. hemisphere and the earthquake or station is above the northern boundary, find the intersection point
                                if ((inflection[2] < 0.) and (p2_lat > n)) or ((inflection[2] < 0.) and (p1_lat > n)):
                                    int = geo_math.known_lat(az, p1_lat, p1_lon, n)
                                    # check if the intersection point is between the east/west bounds and keep it if it is:
                                    check_lon_extent(int, 'eq', arr_sub_raypath_segment)
                        
                        ## SECOND, FIND THE EAST/WEST BOUNDARIES:
                        a11 = geo_math.coord2cart(entry_lat, entry_lon)
                        a12 = geo_math.coord2cart(exit_lat, exit_lon)
                        e_ints = find_gc_intersections(a11, a12, n, e, s, e)
                        if len(e_ints) > 0.:
                            for i in e_ints:
                                intersection_dist = geo_math.GCP_length(eq_lat, eq_lon, i[0], i[1])
                                if (intersection_dist >= arr_sub_raypath_segment[-1, 2]) or (intersection_dist <= arr_sub_raypath_segment[0, 2]):
                                    pass
                                else:
                                    intersection_lats.append(i[0])
                                    intersection_lons.append(i[1])
                                    intersection_dists.append(intersection_dist)
                                    
                                    nearest_int_ids = find_nearest(arr_sub_raypath_segment[:, 2], intersection_dist, 'bounds')                                    
                                    int_dist_1 = arr_sub_raypath_segment[nearest_int_ids[0], 2]
                                    int_dist_2 = arr_sub_raypath_segment[nearest_int_ids[1], 2]
                                    int_depth_1 = arr_sub_raypath_segment[nearest_int_ids[0], 3]
                                    int_depth_2 = arr_sub_raypath_segment[nearest_int_ids[1], 3]
                                    
                                    intersection_depth = geo_math.new_depth(int_depth_1, int_depth_2, int_dist_1, int_dist_2, intersection_dist)
                                    intersection_depths.append(intersection_depth)
    
                        w_ints = find_gc_intersections(a11, a12, n, w, s, w)
                        for i in w_ints:
                            if len(i) > 0:
                                intersection_dist = geo_math.GCP_length(eq_lat, eq_lon, i[0], i[1])
                                if (intersection_dist >= arr_sub_raypath_segment[-1, 2]) or (intersection_dist <= arr_sub_raypath_segment[0, 2]):
                                    pass
                                else:
                                    intersection_lats.append(i[0])
                                    intersection_lons.append(i[1])
                                    intersection_dists.append(intersection_dist)
                                    
                                    nearest_int_ids = find_nearest(arr_sub_raypath_segment[:, 2], intersection_dist, 'bounds')                                    
                                    int_dist_1 = arr_sub_raypath_segment[nearest_int_ids[0], 2]
                                    int_dist_2 = arr_sub_raypath_segment[nearest_int_ids[1], 2]
                                    int_depth_1 = arr_sub_raypath_segment[nearest_int_ids[0], 3]
                                    int_depth_2 = arr_sub_raypath_segment[nearest_int_ids[1], 3]
                                    
                                    intersection_depth = geo_math.new_depth(int_depth_1, int_depth_2, int_dist_1, int_dist_2, intersection_dist)
                                    intersection_depths.append(intersection_depth)
                        
                        # by now, if there are still no boundaries detected, the segment is either entirely inside or entirely outside of the box
                        if len(intersection_lats) == no_intersection_len:
                            # find out if it's entirely inside the box:
                            if (n > entry_lat > s) and ((w < entry_lon) or (e > entry_lon)):
                                int_1 = [entry_lat, entry_lon, entry_dist, entry_depth]
                                int_2 = [exit_lat, exit_lon, exit_dist, exit_depth]
                                path_intersections = np.vstack((path_intersections, int_1, int_2))
    
                                segment_length = raypath_length(arr_sub_raypath_segment).sum()
                                length_in_box += segment_length
                                # az_in_box.append(geo_math.azimuth(p1_lat, p1_lon, p2_lat, p2_lon))

                                with open(f'./event_files/{eq_id}.csv', 'a') as fout:
                                    fout.write(f'{sta_lat},{sta_lon},{phase},{p1_lat},{p1_lon},{p1_dist},{p1_depth},{p2_lat},{p2_lon},{p2_dist},{p2_depth},{segment_length}\n')

                        # if boundaries have been found, find out which segments are in the box and how long they are:
                        elif len(intersection_lats) > no_intersection_len:
                            intersection_rounded_lats = [round(val, 2) for val in intersection_lats]
                            intersection_rounded_lons = [round(val, 2) for val in intersection_lons]
                            intersection_rounded_dists = [round(val, 2) for val in intersection_dists]
                            intersection_rounded_depths  = [round(val, 2) for val in intersection_depths]
                            
                            unsorted_rounded_ar = np.array([intersection_rounded_lats, intersection_rounded_lons, intersection_rounded_dists, intersection_rounded_depths]).T
                            unsorted_ar = np.array([intersection_lats, intersection_lons, intersection_dists, intersection_depths]).T
                            sorted_rounded_ar = unsorted_rounded_ar[unsorted_rounded_ar[:, 2].argsort()]
                            sorted_ar = unsorted_ar[unsorted_ar[:, 2].argsort()]
                            unique_ids = np.unique(sorted_rounded_ar[:, 2], return_index = True)[1]
                            intersection_ar = sorted_ar[unique_ids]
                            
                            # find the mid point of the current segment
                            for i in range(len(intersection_ar) - 1):
                                p1_lat = intersection_ar[i, 0]
                                p1_lon = intersection_ar[i, 1]
                                p1_dist = intersection_ar[i, 2]
                                p1_depth = intersection_ar[i, 3]
                            
                                p2_lat = intersection_ar[i + 1, 0]
                                p2_lon = intersection_ar[i + 1, 1]
                                p2_dist = intersection_ar[i + 1, 2]
                                p2_depth = intersection_ar[i + 1, 3]
                            
                            # find the distance of the segment and its midpoint
                            mid_dist = p2_dist - p1_dist
                            mid_pt = geo_math.GCP_point(p1_lat, p1_lon, p2_lat, p2_lon, mid_dist, (mid_dist / 2))
    
                            # if the midpoint is in the box, add to cumulative dist
                            if (n > mid_pt[0] > s) and ((w < mid_pt[1]) or (e > mid_pt[1])):
                                
                                nearest_start_ids = find_nearest(arr_sub_raypath_segment[:, 2], p1_dist, 'bounds')
                                nearest_end_ids = find_nearest(arr_sub_raypath_segment[:, 2], p2_dist, 'bounds')
                                intersecting_segment_array = arr_sub_raypath_segment[nearest_start_ids[0]:nearest_end_ids[-1] + 1]
                                
                                intersecting_segment_array[0, 0] = p1_lat
                                intersecting_segment_array[0, 1] = p1_lon
                                intersecting_segment_array[0, 2] = p1_dist
                                intersecting_segment_array[0, 3] = p1_depth
                                
                                intersecting_segment_array[-1, 0] = p2_lat
                                intersecting_segment_array[-1, 1] = p2_lon
                                intersecting_segment_array[-1, 2] = p2_dist
                                intersecting_segment_array[-1, 3] = p2_depth
    
                                int_1 = [p1_lat, p1_lon, p1_dist, p1_depth]
                                int_2 = [p2_lat, p2_lon, p2_dist, p2_depth]
                                path_intersections = np.vstack((path_intersections, int_1, int_2))
    
                                segment_length = raypath_length(intersecting_segment_array).sum()
                                length_in_box += segment_length
                                # az_in_box.append(geo_math.azimuth(p1_lat, p1_lon, p2_lat, p2_lon))

                                with open(f'./event_files/{eq_id}.csv', 'a') as fout:
                                    fout.write(f'{sta_lat},{sta_lon},{phase},{p1_lat},{p1_lon},{p1_dist},{p1_depth},{p2_lat},{p2_lon},{p2_dist},{p2_depth},{segment_length}\n')
                                
                    idx_start = zd_idx

        percent_in_box = (length_in_box / total_length) * 100
        # if len(az_in_box) > 0:
            # mean_az_in_box = np.mean(np.array(az_in_box))
            # std_az_in_box = np.std(np.array(az_in_box))
        else:
            # mean_az_in_box = np.nan
            # std_az_in_box = np.nan
        df_data.loc[line, 'AZIMUTH'] = path_az
        # df_data.loc[line, 'MEAN_AZ_IN_BOX'] = mean_az_in_box
        # df_data.loc[line, 'STD_AZ_IN_BOX'] = std_az_in_box
        df_data.loc[line, 'MAX_DEPTH_KM'] = max_depth
        df_data.loc[line, 'TOTAL_LENGTH_KM'] = total_length
        df_data.loc[line, 'LENGTH_IN_BOX_KM'] = length_in_box
        df_data.loc[line, '%_IN_BOX'] = percent_in_box

    except:
        pass

df_geo_data = df_data.drop(idx_to_drop).reset_index(drop = True)
df_geo_data.to_csv(box_sampling_input.dataset_save_name, index = False)


