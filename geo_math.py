##########
# 
# Assorted geographical functions
# written by Claire Richardson (crricha5@asu.edu)
#
# READ ME:
#
# This module contains various functions for geographical applications on a spherical body.
# These are written with Earth in mind and some functions are written to be Earth specific,
# but can be directly used or adapted for other spherical applications (e.g., other
# planetary bodies). Natively, these function use four input variables from `box_sampling_input.py`,
# however these can be modified by the user for convenience/if desired. They are listed at
# the top of the page, beneath the import statements. The spherical trig used in many of
# these functions was sourced from my favorite †extbook of all time:
#
# "Spherical Trigonometry: For the Use of Colleges and Schools" by Isaac Todhunter (1886).
# link to open source PDF from Project Gutenberg: https://www.gutenberg.org/ebooks/19770
#
# IMPORTANT!!!!: THESE FUNCTIONS WERE WRITTEN ASSUMING A -180/180 LONGITUDE CONVENTION!!!!
# Some will work with 0/360, but others will not. I recommend converting your coordinates
# to -180/180 before using these functions if they aren't already to ensure proper
# performace.
#
# These functions are certainly not exhaustive of the possible geographical scenarios
# one may find themselves needing to solve. This module arose from a specific research
# project I was working on during my PhD (the development of the SITRUS package, for which 
# this module is necessary). The configurations of these functions, i.e., my choice of
# arguments and and outputs, is consistent with what was convenient for the development of
# that package. Still, I hope that they are useful and convenient for others, as I went
# through much pain to learn and parse the trig that constitutes many of them. Some are
# fairly ugly, which partly falls out of the fact that I've attemped to inclulde all of
# the boundary cases I could think of.
#
# Finally, there are surely redundancies with other available tools, but I wrote many of
# these out of frustration of not being able to find what I needed. Similarly, I wrote the
# docstrings with an emphasis on clarity at the expense of verbosity, also based on past
# frustrations with too-sparse documentation. All of this being said, this is some of the
# first-ever code that I wrote. I used these to teach myself Python, and programming in
# general, so if there are any issues PLEASE let me know! You can email me at
# crricha5@asu.edu. Thanks very much :)
#
# the functions in this module include:
# 1.  path_length: compute the distance in km between two points at depth
# 2.  GCP_length: compute the epicentral distance bewteen two points
# 3.  coord2cart: convert lat/lon to cartesian coordinates
# 4.  cart2coord: convert cartesian coordinates to lat/lon
# 5.  GCP_point: find the coordinates of a point between two others along a
#         great circle (GC)
# 6.  azimuth: find the azimuth between two points
# 7.  new_depth: find the depth of a point between two others along a GC
# 8.  midpt_depth: find the depth of the midpoint between two given depths
# 9.  pt_from_dist: find the lat/lon of a point at a given epicentral distance from a
#         starting point
# 10. known_lat: find the lat/lon/distance of a point whose lat is known, given the
#         azimuth and distance from a starting point
# 11. known_lon: find the lat/lon/distance of a point whose lon is known, given the
#         azimuth and distance from a starting point
# 12. inflection_finder: find whether there is an inflection along the GC between two
#         points, and if there is, find its lat/lon/distance from the starting point
#
##########

import math
import box_sampling_input

total_radius = 6371 # radius of the Earth (or other planetary body) in kilometers
CMB_depth = 2891 # depth of the core-mantle boundary in kilometers
cdp = 10 # added to avoid rounding errors and inconsistencies. default: 10
rdp = 5 # added to avoid rounding errors and inconsistencies. default: 5

## law of cosines with differing depths:
def path_length(z1, z2, phi1, phi2):
    '''
    Function to calculate the length of a 1D ray segment between two points at the same or different depths in a sphere using the law of cosines.
    Can be applied to any two points along a great circle, i.e., the starting point can be either at zero distance or at some positive non-zero distance from the starting point if the extent of a smaller sub-segment is desired.
    ======
    Inputs:
    - `z1`: depth of point 1 (start point of segment) [format: float or int; unit: km]
    - `z2`: depth of point 2 (endpoint of segment) [format: float or int; unit: km]
    - `phi1`: distance along the entire path segment of point 1 (can be zero) [format: float or int; unit: degrees]
    - `phi2`: distance along the entire path segment of point 2 [format: float or int; unit: degrees]
    ======     
    Outputs:
    - path length: the length of the path segment between points 1 & 2 [format: float; unit: km]
    ======
    '''
    if round(z1, rdp) == round(z2, rdp) == CMB_depth:
        length = ((2. * math.pi * (total_radius - CMB_depth)) / 360.) * (phi2 - phi1)
    else:    
        r1 = total_radius - z1
        r2 = total_radius - z2
        phi_deg = phi2 - phi1
        phi_rad = math.radians(phi_deg)
        length = math.sqrt((r1**2) + (r2**2) - (2 * r1 * r2 * math.cos(phi_rad)))
    return float(length)
        
## length in degrees along a path, calculated from the haversine formula:
def GCP_length(lat1, lon1, lat2, lon2):
    '''
    Function to calculate the epicentral distance between two points along a great circle using the haversine formula.
    ======
    Inputs:
    - `lat1`: latitude of point 1 [format: float or int; unit: degrees]
    - `lon1`: longitude of point 1 [format: float or int; unit: degrees]
    - `lat2`: latitude of point 2 [format: float or int; unit: degrees]
    - `lon2`: longitude of point 2 [format: float or int; unit: degrees]
    ======
    Outputs:
    - epicentral distance: epicentral distance between points 1 & 2 [format: float; unit: degrees]
    ======
    '''
    lon1 += 180.
    lon2 += 180.
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)
    a = (math.sin((lat2_rad - lat1_rad) / 2) ** 2) + math.cos(lat1_rad) * math.cos(lat2_rad) * (math.sin((lon2_rad - lon1_rad) / 2) ** 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return round(math.degrees(c), cdp)

## convert lat/lon to cartesian:
def coord2cart(latitude, longitude):
    '''
    Function to convert geographical coordinates to cartesian coordinates (x, y, z) on a unit sphere.
    ======
    Inputs:
    - `latitude`: latitude of the coordinates you wish to convert
    - `longitude`: longitude of the coordinates you wish to convert
    ======
    Outputs:
    - unit cartesian coordinates corresponding to `latitude` and `longitude` [format: list of three floats [x, y, z]; unit: [unitless, unitless, unitless]]
    ======
    '''
    lat_rad = math.radians(latitude)
    lon_rad = math.radians(longitude)
    x = round(math.cos(lon_rad) * math.cos(lat_rad), cdp)
    y = round(math.sin(lon_rad) * math.cos(lat_rad), cdp)
    z = round(math.sin(lat_rad), cdp)
    return [float(x), float(y), float(z)]

## convert cartesian to lat/lon:
def cart2coord(x, y, z):
    '''
    Function to convert unit cartesian coordinates (x, y, z) on a unit sphere to geographical coordinates
    ======
    Inputs:
    - `x`: coordinate along the x-axis of the unit sphere (value must be from -1 to 1) [format: float or int; unit: none]
    - `y`: coordinate along the y-axis of the unit sphere (value must be from -1 to 1) [format: float or int; unit: none]
    - `z`: coordinate along the z-axis of the unit sphere (value must be from -1 to 1) [format: float or int; unit: none]
    ======
    Outputs:
    - geographical coordinates corresponding to (`x`, `y`, `z`) [format: list of two floats [lat, lon]; unit: [degrees, degrees]]
    ======
    '''
    lat = round(math.degrees(math.asin(z / 1.)), cdp)
    lon = round(math.degrees(math.atan2(y, x)), cdp)
    if lon >= 180.:
        lon -= 360.
    if lon < -180:
        lon += 360
    return [float(lat), float(lon)]

## find the coordinates of a point between two other points on the same great circle path
def GCP_point(lat1, lon1, lat2, lon2, total_epidist, int_epidist):
    '''
    Function to find the coordinates of a point along a great circle path between two other points on the GCP, given a target intermediate distance from point 1.
    ======
    Inputs:
    - `lat1`: latitude of point 1 [format: float or int; unit: degrees]
    - `lon1`: longitude of point 1 [format: float or int; unit: degrees]
    - `lat2`: latitude of point 2 [format: float or int; unit: degrees]
    - `lon2`: longitude of point 2 [format: float or int; unit: degrees]
    - `total_epidist`: total epicentral distance between points 1 & 2 [format: float or int; unit: degrees]
    - `int_epidist`: intermediate distance between point 1 and the point whose coordinates are to be determined [format: float or int; unit: degrees]
    ======
    Outputs:
    - geographical coordinats at the point `int_epidist` degrees from point 1 [format: list of two floats [lat, lon]; unit: [degrees, degrees]]
    ======
    '''
    # first, convert bound points to cartesian:
    cart1 = coord2cart(lat1, lon1)
    cart2 = coord2cart(lat2, lon2)
    x1, y1, z1 = cart1[0], cart1[1], cart1[2]
    x2, y2, z2 = cart2[0], cart2[1], cart2[2]
    # second, calculate the third point in cartesian:
    epidist_rad = math.radians(total_epidist)
    phi_rad = math.radians(int_epidist)
    x3 = ((x1 * math.sin(epidist_rad - phi_rad)) + (x2 * math.sin(phi_rad))) / (math.sin(epidist_rad))
    y3 = ((y1 * math.sin(epidist_rad - phi_rad)) + (y2 * math.sin(phi_rad))) / (math.sin(epidist_rad))
    z3 = ((z1 * math.sin(epidist_rad - phi_rad)) + (z2 * math.sin(phi_rad))) / (math.sin(epidist_rad))
    # third, convert the third point to lat/lon:
    new_coord = cart2coord(x3, y3, z3)
    lat3 = new_coord[0]
    lon3 = new_coord[1]
    return [lat3, lon3]
    
## calculate azimuth of a line segment between geographical coordinates
def azimuth(lat1, lon1, lat2, lon2):
    '''
    Function to compute the azimuth of a great circle path between two points.
    ======
    Inputs:
    - `lat1`: latitude of point 1 [format: float or int; unit: degrees]
    - `lon1`: longitude of point 1 [format: float or int; unit: degrees]
    - `lat2`: latitude of point 2 [format: float or int; unit: degrees]
    - `lon2`: longitude of point 2 [format: float or int; unit: degrees]
    ======
    Outputs:
    - azimuth between points 1 & 2 [format: float; unit: degrees between 0 and 360]
    ======
    '''
    if round(lat1, rdp) == 90. and round(lat2, rdp) < round(lat1, rdp):
        az = 180.
    elif round(lat1, rdp) == -90. and round(lat2, rdp) > round(lat1, rdp):
        az = 0.
    else:
        delta_lat = lat2 - lat1
        delta_lon = lon2 - lon1
        lat1_rad = math.radians(lat1)
        lat2_rad = math.radians(lat2)
        lon1_rad = math.radians(lon1)
        lon2_rad = math.radians(lon2)
        delta_lat_rad = math.radians(delta_lat)
        delta_lon_rad = math.radians(delta_lon)
        az_rad = math.atan2((math.sin(delta_lon_rad) * math.cos(lat2_rad)), (math.cos(lat1_rad) * math.sin(lat2_rad) - math.sin(lat1_rad) * math.cos(lat2_rad) * math.cos(delta_lon_rad)))
        az = math.degrees(az_rad)
    if round(az, rdp) < 0.:
        az = 360. - abs(az)
    if round(az, rdp) == 360.:
        az = 0.
    return round(az, cdp)

## calculate the depth of a new point between two other points, given the epicentral distance to the new point
def new_depth(depth1, depth2, dist1, dist2, new_dist):
    '''
    Function to compute the depth of a point between two other points at depth in a sphere along a great circle path.
    Finds the linear depth, i.e., doesn't take ray parameter into account.
    ======
    Inputs:
    - `depth1`: depth of point 1 (start point) [format: float or int; unit: km/m/cm/etc.]
    - `depth2`: depth of point 2 (end point) [format: float or int; unit: km/m/cm/etc.]
    - `dist1`: epicentral distance of point 1 along the greater raypath (can be 0) [format: float or int; unit: degrees]
    - `dist2`: epicentral distance of point 2 along the greater raypath [format: float or int; unit: degrees]
    - `new_dist`: epicentral distance of the intermediate point of interest, relative to the greater raypath [format: float or int; unit: degrees]
    ======
    Outputs:
    - depth of the intermediate point along a line between points 1 & 2 [format: float; unit: km/m/cm/etc.]
    ======
    '''
    if round(depth1, rdp) == round(depth2, rdp):
        return depth1
    else:        
        ratio = (new_dist - dist1) / (dist2 - dist1)
        z = ((depth2 - depth1) * ratio) + depth1
        return round(z, cdp)

## calculate the depth of the midpoint between two points
def midpt_depth(z1, z2):
    '''
    Function to find the linear mid-depth between two other depths.
    ======
    Inputs:
    - `z1`: depth of point 1 [format: float or int; unit: km/m/cm/etc.]
    - `z2`: depth of point 2 [format: float or int; unit: km/m/cm/etc.]
    ======
    Outputs:
    - mid-depth between depths 1 and 2 [format: float; unit: km/m/cm/etc.]
    '''
    if round(z1, rdp) == round(z2, rdp):
        return round(z1, cdp)
    else:
        z3 = ((z2 - z1) / 2) + z1
        return round(z3, cdp)

###################
# SPHERICAL TRIG: #
###################

## Function to compute the latitude and longitude of a desired point on a great circle, given one point and the distance from it
def pt_from_dist(azimuth, lat1, lon1, dist):
    '''
    Function to compute the latitude and longitude of a desired point some distance from a given point along a great circle.
    ======
    Inputs:
    - `azimuth`: azimuth along which to find the new point [format: float or int; unit: degrees]
    - `lat1`: latitude of the starting point [format: float or int; unit: degrees]
    - `lon1`: longitude of the starting point [format: float or int; unit: degrees]
    - `dist`: epicentral distance along the great circle defined by `azimuth` from the starting point to find the new point [format: float or int; unit: degrees]
    ======
    Outputs:
    - geographical coordinats at the point `dist` degrees from the starting point [format: list of two floats [lat, lon]; unit: [degrees, degrees]]
    ======
    '''
    azimuth = float(azimuth)
    lat1 = float(lat1)
    lon1 = float(lon1)
    
    if round(azimuth, rdp) == 360. or round(azimuth, rdp) == 0.:
        new_lat = lat1 + dist
        new_lon = lon1
        if new_lat > 90.:
            diff = new_lat - 90.
            new_lat = 90. - diff
            if lon1 < 0.:
                new_lon = lon1 + 180.
            elif lon1 >= 0.:
                new_lon = lon1 - 180.
                
    elif round(azimuth, rdp) == 180.:
        new_lat = lat1 - dist
        new_lon = lon1
        if new_lat < -90.:
            diff = -90. - new_lat
            new_lat = -90. + diff
            if lon1 < 0.:
                new_lon = lon1 + 180.
            elif lon1 >= 0.:
                new_lon = lon1 - 180.
                
    # if path is eastbound:
    elif 0. < round(azimuth, rdp) < 180.:
        B_ang = math.radians(azimuth)
        c_side = math.radians(90. - lat1)
        a_side = math.radians(dist)
        x = math.atan((math.cos(0.5 * (a_side - c_side)) / math.cos(0.5 * (a_side + c_side))) * (1 / math.tan(0.5 * B_ang)))
        y = math.atan((math.sin(0.5 * (a_side - c_side)) / math.sin(0.5 * (a_side + c_side))) * (1 / math.tan(0.5 * B_ang)))
        A_ang = x + y
        if A_ang < 0.:
            A_ang += math.pi
        b_side = math.acos((math.cos(c_side) * math.cos(a_side)) + (math.sin(c_side) * math.sin(a_side) * math.cos(B_ang)))
        new_lat = 90. - math.degrees(b_side)
        new_lon = lon1 + math.degrees(A_ang)
        
    # if the path is westbound:
    elif 180. < round(azimuth, rdp) < 360.:
        B_ang = math.radians(360. - azimuth)
        c_side = math.radians(90. - lat1)
        a_side = math.radians(dist)
        x = math.atan((math.cos(0.5 * (a_side - c_side)) / math.cos(0.5 * (a_side + c_side))) * (1 / math.tan(0.5 * B_ang)))
        y = math.atan((math.sin(0.5 * (a_side - c_side)) / math.sin(0.5 * (a_side + c_side))) * (1 / math.tan(0.5 * B_ang)))
        A_ang = x + y
        if A_ang < 0.:
            A_ang += math.pi
        b_side = math.acos((math.cos(c_side) * math.cos(a_side)) + (math.sin(c_side) * math.sin(a_side) * math.cos(B_ang)))
        new_lat = 90. - math.degrees(b_side)
        new_lon = lon1 - math.degrees(A_ang)
        
    if round(new_lon, rdp) >= 180.:
        new_lon -= 360.
    if round(new_lon, rdp) < -180:
        new_lon += 360

    return [round(new_lat, cdp), round(new_lon, cdp)]


def known_lat(az, lat1, lon1, bound_lat):
    '''
    Function to compute the geographical coordinates and epicentral distance of a point of interest along the great circle path from another known point, given the azimuth and latitude of the point of interest.
    ======
    Inputs:
    - `az`: azimuth along which the starting point and the desired point lie [format: float or int; unit: degrees]
    - `lat1`: latitude of the starting point [format: float or int; unit: degrees]
    - `lon1`: longitude of the starting point [format: float or int; unit: degrees]
    - `bound_lat`: the latitude of the point of interest [format: float or int; unit: degrees]
    ======
    Outputs:
    - geographical coordinates of the point of interest and the epicentral ditance from the starting point to the point of interest [format: list of three floats [lat, lon, dist]; unit: [degrees, degrees, degrees]]
    ======
    '''
    az = float(az)
    lat1 = float(lat1)
    lon1 = float(lon1)
    bound_lat = float(bound_lat)
    b_side = math.radians(90. - bound_lat)
    c_side = math.radians(90. - lat1)
    if 0. < round(az, rdp) < 90.:
        B_ang = math.radians(az)
        C_ang = math.pi - (math.asin((math.sin(c_side) / math.sin(b_side)) * math.sin(B_ang)))
        A_ang = 2 * math.atan((math.cos(0.5 * (b_side - c_side)) / math.cos(0.5 * (b_side + c_side))) * (1 / (math.tan(0.5 * (B_ang + C_ang)))))
        a_side = 2 * math.atan((math.cos(0.5 * (B_ang + C_ang)) / math.cos(0.5 * (B_ang - C_ang))) * (math.tan(0.5 * (b_side + c_side))))
        bound_dist = math.degrees(a_side)
        bound_lon = lon1 + math.degrees(A_ang)

    elif 90. < round(az, rdp) < 180.:
        B_ang = math.radians(az)
        C_ang = math.asin((math.sin(c_side) / math.sin(b_side)) * math.sin(B_ang))
        A_ang = 2 * math.atan((math.cos(0.5 * (b_side - c_side)) / math.cos(0.5 * (b_side + c_side))) * (1 / (math.tan(0.5 * (B_ang + C_ang)))))
        a_side = 2 * math.atan((math.cos(0.5 * (B_ang + C_ang)) / math.cos(0.5 * (B_ang - C_ang))) * (math.tan(0.5 * (b_side + c_side))))
        bound_dist = math.degrees(a_side)
        bound_lon = lon1 + math.degrees(A_ang)

    elif 180. < round(az, rdp) < 270.:
        B_ang = math.radians(360. - az)
        C_ang = math.asin((math.sin(c_side) / math.sin(b_side)) * math.sin(B_ang))
        A_ang = 2 * math.atan((math.cos(0.5 * (b_side - c_side)) / math.cos(0.5 * (b_side + c_side))) * (1 / (math.tan(0.5 * (B_ang + C_ang)))))
        a_side = 2 * math.atan((math.cos(0.5 * (B_ang + C_ang)) / math.cos(0.5 * (B_ang - C_ang))) * (math.tan(0.5 * (b_side + c_side))))
        bound_dist = math.degrees(a_side)
        bound_lon = lon1 - math.degrees(A_ang)

    elif 270. < round(az, rdp) < 360.:
        B_ang = math.radians(360. - az)
        C_ang = math.pi - (math.asin((math.sin(c_side) / math.sin(b_side)) * math.sin(B_ang)))
        A_ang = 2 * math.atan((math.cos(0.5 * (b_side - c_side)) / math.cos(0.5 * (b_side + c_side))) * (1 / (math.tan(0.5 * (B_ang + C_ang)))))
        a_side = 2 * math.atan((math.cos(0.5 * (B_ang + C_ang)) / math.cos(0.5 * (B_ang - C_ang))) * (math.tan(0.5 * (b_side + c_side))))           
        bound_dist = math.degrees(a_side)
        bound_lon = lon1 - math.degrees(A_ang)

    elif round(az, rdp) == 0. or round(az, rdp) == 360.:
        bound_lon = lon1
        bound_dist = bound_lat - lat1

    elif round(az, rdp) == 90.:
        B_ang = math.radians(az)
        a_side = math.acos(math.cos(b_side) / math.cos(c_side))
        A_ang = math.acos(math.tan(c_side) / math.tan(b_side))
        bound_lon = lon1 + math.degrees(A_ang)
        bound_dist = math.degrees(a_side)

    elif round(az, rdp) == 180.:
        bound_lon = lon1
        bound_dist = lat1 - bound_lat

    elif round(az, rdp) == 270.:
        B_ang = math.radians(360. - az)
        a_side = math.acos(math.cos(b_side) / math.cos(c_side))
        A_ang = math.acos(math.tan(c_side) / math.tan(b_side))
        bound_lon = lon1 - math.degrees(A_ang)
        bound_dist = math.degrees(a_side)

    if round(bound_lon, rdp) >= 180.:
        bound_lon -= 360.
    if round(bound_lon, rdp) < -180:
        bound_lon += 360
    return [round(bound_lat, cdp), round(bound_lon, cdp), round(bound_dist, cdp)]
    
def known_lon(az, lat1, lon1, bound_lon):
    '''
    Function to compute the geographical coordinates and epicentral distance of a point of interest along the great circle path from another known point, given the azimuth and longitude of the point of interest.
    ======
    Inputs:
    - `az`: azimuth along which the starting point and the desired point lie [format: float or int; unit: degrees]
    - `lat1`: latitude of the starting point [format: float or int; unit: degrees]
    - `lon1`: longitude of the starting point [format: float or int; unit: degrees]
    - `bound_lon`: the longitude of the point of interest [format: float or int; unit: degrees]
    ======
    Outputs:
    - geographical coordinates of the point of interest and the epicentral ditance from the starting point to the point of interest [format: list of three floats [lat, lon, dist]; unit: [degrees, degrees, degrees]]
    ======
    '''
    az = float(az)
    lat1 = float(lat1)
    lon1 = float(lon1)
    bound_lon = float(bound_lon)
    c_side = math.radians(90. - lat1)
    if 0. < round(az, rdp) < 90.:
        B_ang = math.radians(az)
        if round(bound_lon, rdp) < 0. and round(lon1, rdp) > 0.:
            diff1 = bound_lon + 180.
            diff2 = 180. - lon1
            A_ang = math.radians(diff1 + diff2)
        else:
            A_ang = math.radians(bound_lon - lon1)
        x = math.atan((math.cos(0.5 * (A_ang - B_ang)) / math.cos(0.5 * (A_ang + B_ang))) * math.tan(0.5 * c_side))
        y = math.atan((math.sin(0.5 * (A_ang - B_ang)) / math.sin(0.5 * (A_ang + B_ang))) * math.tan(0.5 * c_side))
        a_side = x + y
        b_side = x - y
        bound_dist = math.degrees(a_side)
        bound_lat = 90 - math.degrees(b_side)
    elif 90. < round(az, rdp) < 180.:
        B_ang = math.radians(az)
        if round(bound_lon, rdp) < 0. and round(lon1, rdp) > 0.:
            diff1 = bound_lon + 180.
            diff2 = 180. - lon1
            A_ang = math.radians(diff1 + diff2)
        else:
            A_ang = math.radians(bound_lon - lon1)
        x = math.atan((math.cos(0.5 * (A_ang - B_ang)) / math.cos(0.5 * (A_ang + B_ang))) * math.tan(0.5 * c_side))
        y = math.atan((math.sin(0.5 * (A_ang - B_ang)) / math.sin(0.5 * (A_ang + B_ang))) * math.tan(0.5 * c_side))
        a_side = x + y
        b_side = x - y
        if a_side < 0.:
            a_side += math.pi
        if b_side < 0.:
            b_side = math.pi + b_side
        bound_dist = math.degrees(a_side)
        bound_lat = 90. - math.degrees(b_side)
    elif 180. < round(az, rdp) < 270.:
        B_ang = math.radians(360. - az)
        if round(bound_lon, rdp) > 0. and round(lon1, rdp) < 0.:
            diff1 = lon1 + 180.
            diff2 = 180. - bound_lon
            A_ang = math.radians(diff1 + diff2)
        else:
            A_ang = math.radians(lon1 - bound_lon)
        x = math.atan((math.cos(0.5 * (A_ang - B_ang)) / math.cos(0.5 * (A_ang + B_ang))) * math.tan(0.5 * c_side))
        y = math.atan((math.sin(0.5 * (A_ang - B_ang)) / math.sin(0.5 * (A_ang + B_ang))) * math.tan(0.5 * c_side))
        a_side = x + y
        b_side = x - y
        if a_side < 0.:
            a_side += math.pi
        if b_side < 0.:
            b_side = math.pi + b_side
        bound_dist = math.degrees(a_side)
        bound_lat = 90. - math.degrees(b_side)
    elif 270. < round(az, rdp) < 360.:
        B_ang = math.radians(360. - az)
        if round(bound_lon, rdp) > 0. and round(lon1, rdp) < 0.:
            diff1 = lon1 + 180.
            diff2 = 180. - bound_lon
            A_ang = math.radians(diff1 + diff2)
        else:
            A_ang = math.radians(lon1 - bound_lon)
        x = math.atan((math.cos(0.5 * (A_ang - B_ang)) / math.cos(0.5 * (A_ang + B_ang))) * math.tan(0.5 * c_side))
        y = math.atan((math.sin(0.5 * (A_ang - B_ang)) / math.sin(0.5 * (A_ang + B_ang))) * math.tan(0.5 * c_side))
        a_side = x + y
        b_side = x - y
        bound_dist = math.degrees(a_side)
        bound_lat = 90 - math.degrees(b_side)
    elif round(az, rdp) == 90.:
        B_ang = math.radians(az)
        if round(bound_lon, rdp) < 0. and round(lon1, rdp) > 0.:
            diff1 = bound_lon + 180.
            diff2 = 180. - lon1
            A_ang = math.radians(diff1 + diff2)
        else:
            A_ang = math.radians(bound_lon - lon1)
        a_side = math.atan(math.tan(A_ang) * math.sin(c_side))
        b_side = math.atan(math.tan(c_side) / math.cos(A_ang))
        if a_side < 0.:
            a_side += math.pi
        if b_side < 0.:
            b_side = math.pi + b_side
        bound_dist = math.degrees(a_side)
        bound_lat = 90. - math.degrees(b_side)
    elif round(az, rdp) == 270.:
        B_ang = math.radians(360. - az)
        if round(bound_lon, rdp) > 0. and round(lon1, rdp) < 0.:
            diff1 = lon1 + 180.
            diff2 = 180. - bound_lon
            A_ang = math.radians(diff1 + diff2)
        else:
            A_ang = math.radians(lon1 - bound_lon)
        a_side = math.atan(math.tan(A_ang) * math.sin(c_side))
        b_side = math.atan(math.tan(c_side) / math.cos(A_ang))
        if a_side < 0.:
            a_side += math.pi
        if b_side < 0.:
            b_side = math.pi + b_side
        bound_dist = math.degrees(a_side)
        bound_lat = 90. - math.degrees(b_side)
    if round(bound_lon, rdp) >= 180.:
        bound_lon -= 360.
    if round(bound_lon, rdp) < -180:
        bound_lon += 360
    return [round(bound_lat, cdp), round(bound_lon, cdp), round(bound_dist, cdp)]


def inflection_finder(az, baz, lat1, lon1, lat2, lon2, dist1, dist2):
    '''
    Function to find whether there is an inflection between two points on a geographical sphere (i.e., whether a raypath traveling along a great circle changes heading from north to south or vice/versa along its route. this is mostly applicable for raypaths near the poles).
    ======
    Inputs:
    - `az`: azimuth along which the raypath points, from point 1 to point 2 [format: float or int; unit: degrees]
    - `baz`: back azimuth along which the raypath lies, from point 2 to point 1 [format: float or int; unit: degrees]
    - `lat1`: latitude of point 1 [format: float or int; unit: degrees]
    - `lon1`: longitude of point 1 [format: float or int; unit: degrees]
    - `lat2`: latitude of point 2 [format: float or int; unit: degrees]
    - `lon2`: longitude of point 2 [format: float or int; unit: degrees]
    - `dist1`: epicentral distance along the greater raypath that point 1 lies (can be 0) [format: float or int; unit: degrees]
    - `dist2`: epicentral distance along the greater raypath that point 2 lies [format: float or int; unit: degrees]
    ======
    Outputs:
    - if there *is* an inflection point:
        - returns True, the epicentral distance along the greater raypath to the inflection point, and the latitude and longitude of the inflection point [format: list of a boolean and three floats [True, inflection_dist, inflection_lat, inflection_lon]; unit: [None, degrees, degrees, degrees]]
    - if there *is not* an inflection point:
        - returns False [format: list of one boolean value [False]; unit: [None]]
    ======
    '''
    az = float(az)
    baz = float(baz)
    # find whether there is an inflection and, if there is, find which quadrant the azimuth is in
    if 0. < round(az, rdp) < 90.:
        quadrant = 'I'
        if 270. < round(baz, rdp) < 360.:
            inflection = True
        else:
            inflection = False
    elif 90. < round(az, rdp) < 180.:
        quadrant = 'II'
        if 180. < round(baz, rdp) < 270.:
            inflection = True
        else:
            inflection = False
    elif 180. < round(az, rdp) < 270.:
        quadrant = 'III'
        if 90. < round(baz, rdp) < 180.:
            inflection = True
        else:
            inflection = False
    elif 270. < round(az, rdp) < 360.:
        quadrant = 'IV'
        if 0. < round(baz, rdp) < 90.:
            inflection = True
        else:
            inflection = False
    elif round(az, rdp) == 0.:
        delta_dist = dist2 - dist1
        dist2pole = GCP_length(lat1, lon1, 90., lon1)
        if delta_dist > dist2pole:
            inflection = True
            quadrant = 'north'
        else:
            inflection = False
            quadrant = None
    elif round(az, rdp) == 180.:
        delta_dist = dist2 - dist1
        dist2pole = GCP_length(lat1, lon1, -90., lon1)
        if delta_dist > dist2pole:
            inflection = True
            quadrant = 'south'
        else:
            inflection = False
            quadrant = None
    elif round(az, rdp) == 90. or round(az, rdp) == 270.:
        inflection = False
        quadrant = None
    # if there is an inflection point, find the info for the point
    if inflection == True:
        c_side = math.radians(90. - lat1)
        # calculate info for infleciton point (depth, dist, lat, lon, shell, block, type)
        if quadrant == 'north':
            inflection_lat = 90.
            inflection_lon = lon1
            inflection_dist = dist1 + dist2pole
        elif quadrant == 'south':
            inflection_lat = -90.
            inflection_lon = lon1
            inflection_dist = dist1 + dist2pole
        elif quadrant == 'I':
            B_ang = math.radians(az)
            C_ang = math.radians(360. - baz)
            inflection_colat = math.degrees(math.asin(math.sin(B_ang) * math.sin(c_side)))
            inflection_ang = math.atan(math.cos(B_ang) * math.tan(c_side))
            inflection_dist = dist1 + math.degrees(inflection_ang)
            inflection_A_prime = math.degrees(math.asin(math.sin(inflection_ang) / math.sin(c_side)))
            inflection_lat = 90. - inflection_colat
            inflection_lon = lon1 + inflection_A_prime

        elif quadrant == 'II':
            B_ang = math.radians(az)
            C_ang = math.radians(360. - baz)
            inflection_colat = 180. - math.degrees(math.asin(math.sin(B_ang) * math.sin(c_side)))
            inflection_ang = math.atan(math.cos(B_ang) * math.tan(c_side))
            inflection_dist = dist1 + math.degrees(inflection_ang)
            inflection_A_prime = math.degrees(math.asin(math.sin(inflection_ang) / math.sin(c_side)))
            inflection_lat = 90. - inflection_colat
            inflection_lon = lon1 + inflection_A_prime

        elif quadrant == 'III':
            B_ang = math.radians(360. - az)
            C_ang = math.radians(baz)
            inflection_colat = 180. - math.degrees(math.asin(math.sin(B_ang) * math.sin(c_side)))
            inflection_ang = math.atan(math.cos(B_ang) * math.tan(c_side))
            inflection_dist = dist1 + math.degrees(inflection_ang)
            inflection_A_prime = math.degrees(math.asin(math.sin(inflection_ang) / math.sin(c_side)))
            inflection_lat = 90. - inflection_colat
            inflection_lon = lon1 - inflection_A_prime

        elif quadrant == 'IV':
            B_ang = math.radians(360. - az)
            C_ang = math.radians(baz)
            inflection_colat = math.degrees(math.asin(math.sin(B_ang) * math.sin(c_side)))
            inflection_ang = math.atan(math.cos(B_ang) * math.tan(c_side))
            inflection_dist = dist1 + math.degrees(inflection_ang)
            inflection_A_prime = math.degrees(math.asin(math.sin(inflection_ang) / math.sin(c_side)))
            inflection_lat = 90. - inflection_colat
            inflection_lon = lon1 - inflection_A_prime
                
        if round(inflection_lon, rdp) >= 180.:
            inflection_lon -= 360.
        if round(inflection_lon, rdp) < -180:
            inflection_lon += 360
        
        return [inflection, round(inflection_dist, cdp), round(inflection_lat, cdp), round(inflection_lon, cdp)]
    # if there isn't an inflection point, just return False
    else:
        return [inflection]