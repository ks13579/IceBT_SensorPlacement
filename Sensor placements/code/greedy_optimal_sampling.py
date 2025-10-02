import numpy as np
import random
from scipy.interpolate import CubicSpline
from scipy.stats import qmc

from sampling_error_calculation import do_sensor_analysis

#Function to generate initial sensor set
def get_inner_sensors(sensor_min_depth,sensor_max_depth):
    
    sample_qmch=np.load('halton_samples_sensors.npy',allow_pickle=True)

    l_bounds = [sensor_min_depth]
    u_bounds = [sensor_max_depth]
    sample_scaled = qmc.scale(sample_qmch, l_bounds, u_bounds)
    
    
    sample_select=[]
    
    while(len(sample_select) != 2):
        si=random.randrange(100000)
        if (sample_scaled[si] != sensor_min_depth and sample_scaled[si] != sensor_max_depth):
            sample_select.append(sample_scaled[si])
        if(len(sample_select) == 2 and sample_select[0] == sample_select[1]):
            sample_select=[]
            print("--inner points same--")
            continue
    return sample_select

#Function to generate optimal sensor locations for the given number of sensors.
def generate_optimized_sensor_locations_candiate_set(T_end, sensor_count, sensor_min_depth, sensor_max_depth, z_grid, used_profiles,sensor_std,error_grid):

    status='okay'
    #    generate candidates set with candidates being all grid points in z direction (since analysis function
    #    currently only returns the errors at the grid points)
    #    candidate_set = z_grid

    initial_sensor_count = 4
    print('initial_sensor_count',initial_sensor_count)
    
    #    set initial sensor set with just 4 sensors as starting point
    #    out of these 4 sensor one is always placed at the top 
    #    and another one is always placed at the bottom of the borehole
    sensors = np.linspace(sensor_min_depth,sensor_max_depth, num=initial_sensor_count, endpoint=True)
    sensor_select=get_inner_sensors(sensor_min_depth,sensor_max_depth)
    sensors[1]=min(sensor_select)
    sensors[2]=max(sensor_select)
    
    print(sensors)
     
    for i in range(sensor_count - initial_sensor_count):
        #    make sure that initially chosen sensors are not considered
    
        #    run analysis with current choice of sensors
        sensors = np.sort(sensors)
        accumulated_error, errors = do_sensor_analysis(T_end,sensors, sensor_std, z_grid, error_grid, used_profiles)

        max_indices = np.argmax(accumulated_error)

        #    add new sensor 
        sensors = np.append(sensors, [error_grid[max_indices]])
        sensors = np.sort(sensors)
        
        #    checks for duplicate sensor locations 
        #    and break the loop to prevents adding more sensors to the same location
        sn_values, counts = np.unique(sensors, return_counts=True)
        if any(counts>1):
            status='duplicate_%f_for_%d_sensors'%(error_grid[max_indices],sensor_count)
            break
   
    return sensors,status
