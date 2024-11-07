##########
#
# Functions to call PREM values at given depths
# 
# READ ME:
# 
# This module contains three functions, each of which contain the polynomials for the elastic
# parameters (density, S-wave velocity, P-wave velocity) defined by the Preliminary Reference
# Earth Model (PREM; Dziewonski and Anderson, 1981).
#
# Functions in this module include:
# 1.  prem_density: get the density of PREM at a given depth
# 2.  prem_vel: get the S-wave or P-wave velocity of PREM at a given depth
#
##########

## function to evaluate PREM density at a given depth:
def prem_density(depth): # units = g / cm^3
    '''
    Function to compute the density of PREM at a given depth
    ======
    Inputs:
    - `depth`: depth at which to find PREM's density
    ======
    Outputs:
    - PREM's density at `depth` [format: float; unit: g / cm^3]
    '''
    a = 6371
    r = a - depth
    x = r / a
    # inner core: (0 - 1221.5 km)
    if 0 < r < 1221.5:
        return 13.0885 - (8.8381 * (x**2))
    # outer core: (1221.5 - 3480.0)
    elif 1221.5 < r < 3480:
        return 12.5815 - (1.2638 * x) - (3.6426 * (x**2)) - (5.5281 * (x**3))
    # lower mantle (3480.0 - 5701)
    elif 3480 <= r < 5701:
        return 7.9565 - (6.4761 * x) + (5.5283 * (x**2)) - (3.0807 * (x**3))
    # transition zone (5701 - 5771)
    elif 5701 <= r < 5771:
        return 5.3197 - (1.4836 * x)
    # transition zone (5771 - 5971)
    elif 5771 <= r < 5971:
        return 11.2494 - (8.0298 * x)
    # transition zone (5971 - 6151)
    elif 5971 <= r < 6151:
        return 7.1089 - (3.8045 * x)
    # LVZ & LID (6151 - 6346.6)
    elif 6151 <= r < 6346.6:
        return 2.6910 + (0.6924 * x)
    # crust (6346.6 - 6356)
    elif 6346.6 <= r < 6356:
        return 2.900
    # crust (6356 - 6371)
    elif 6356 <= r <= 6371:
        return 2.600


## function to evaluate PREM s - or p-wave velocity at a given depth:
def prem_vel(wave_type, depth):
    '''
    Function to compute the P-wave velocity of PREM at a given depth
    ======
    Inputs:
    - `wave_type`: polarity of the wave for which to find the PREM velocity (either 'S' or 'P) [format: string; unit: unitless]
    - `depth`: depth at which to find PREM's P-wave velocity
    ======
    Outputs:
    - PREM's P-wave velocity at `depth` [format: float; unit: km / sec]
    '''
    a = 6371
    r = a - depth
    x = r / a
    if wave_type == 'P':
        # inner core: (0 - 1221.5 km)
        if 0 < r < 1221.5:
            return 11.2622 - (6.3640 * (x**2))
        # outer core: (1221.5 - 3480.0)
        elif 1221.5 < r < 3480:
            return  11.0487 - (4.0362 * x) + (4.8023 * (x**2)) - (13.5732 * (x**3))
        # lower mantle (3480.0 - 3630.0)
        elif 3480 <= r < 3630:
            return  15.3891 - (5.3181 * x) + (5.5242 * (x**2)) - (2.5514 * (x**3))
        # lower mantle (3630.0 - 5600.0)
        elif 3630 <= r < 5600:
            return 24.9520 - (40.4673 * x) + (51.4832 * (x**2)) - (26.6419 * (x**3))
        # lower mantle (5600 - 5701)
        elif 5600 <= r < 5701:
            return 29.2766 - (23.6027 * x) + (5.5242 * (x**2)) - (2.5514 * (x**3))
        # transition zone (5701 - 5771)
        elif 5701 <= r < 5771:
            return 19.0957 - (9.8672 * x)
        # transition zone (5771 - 5971)
        elif 5771 <= r < 5971:
            return 39.7027 - (32.6166 * x)
        # transition zone (5971 - 6151)
        elif 5971 <= r < 6151:
            return 20.3926 - (12.2569 * x)
        # LVZ & LID (6151 - 6346.6)
        elif 6151 <= r < 6346.6:
            return 4.1875 + (3.9382 * x)
        # crust (6346.6 - 6356)
        elif 6346.6 <= r < 6356:
            return 6.800
        # crust (6356 - 6371)
        elif 6356 <= r <= 6371:
            return 5.800

    elif wave_type == 'S':
        # inner core: (0 - 1221.5 km)
        if 0 < r < 1221.5:
            return 3.6678 - (4.4475 * (x**2))
        # outer core: (1221.5 - 3480.0)
        elif 1221.5 < r < 3480:
            return 0
        # lower mantle (depth: 2891 - 2741)
        elif 3480 <= r < 3630:
            return 6.9254 + (1.4672 * x) - (2.0834 * (x**2)) + (0.9783 * (x**3))
        # lower mantle (depth: 2741 - 771)
        elif 3630 <= r < 5600:
            return 11.1671 - (13.7818 * x) + (17.4575 * (x**2)) - (9.2777 * (x**3))
        # lower mantle (depth: 771 - 670)
        elif 5600 <= r < 5701:
            return 22.3459 - (17.2473 * x) - (2.0834 * (x**2)) + (0.9783 * (x**3))
        # transition zone (depth: 670 - 600)
        elif 5701 <= r < 5771:
            return 9.9839 - (4.9324 * x)
        # transition zone (depth: 600 - 400)
        elif 5771 <= r < 5971:
            return 22.3512 - (18.5856 * x)
        # transition zone (depth: 400 - 220)
        elif 5971 <= r < 6151:
            return 8.9496 - (4.4597 * x)
        # LVZ & LID (depth: 220 - 24.4)
        elif 6151 <= r < 6346.6:
            return 2.1519 + (2.3481 * x)
        # crust (depth: 24.4 - 15)
        elif 6346.6 <= r < 6356:
            return 3.900
            # crust (depth: 15 - 0)
        elif 6356 <= r <= 6371:
            return 3.200
        