import numpy as np
import dill
import random
from scipy.interpolate import CubicSpline
from scipy.stats import qmc

from sampling_error_calculation import do_sensor_analysis
from greedy_optimal_sampling import generate_optimized_sensor_locations_candiate_set



random.seed(9001)
    
# array containing final temperature profiles in borehole
T_end2006 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA2_6_simulation_10_ext1000ws.npy')
# array containing surface temperature profiles

T_top2006 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA2_6_surfaceTem_10_ext1000ws.npy')


t_grid = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_sample_times_fwd.npy')
z_grid = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_sample_depths_fwd.npy')

N_z = z_grid.size
N_t = t_grid.size

used_profiles = []
idx = 0
for pa in [2]:
    for beta in [0.6]:
        for k in range(10):
            used_profiles.append(idx)
            idx = idx + 1
borehole_simulations=np.zeros((z_grid.size,10))
borehole_simulations[:,0:10]=T_end2006.T[:,:10]


sensor_min_depth = 0
sensor_max_depth = 200
sensor_count = 20
sensor_std = 0

error_grid_size = 100000
error_grid = np.linspace(sensor_min_depth,sensor_max_depth, num=error_grid_size, endpoint=True)
error_grid_delta = error_grid[1]-error_grid[0]

go_status='okay'

if go_status=='okay':
    print('GRIP_res_%dsensor_10derr_go_avg'%(20))
    no_of_cases=1000
    sensor_depths_go_1000=[]
    
    for n in range(no_of_cases):
        sensor_depths_go_1000.append(0)

    for rep in range(no_of_cases):        
        sensor_interim, status=generate_optimized_sensor_locations_candiate_set(borehole_simulations, sensor_count, sensor_min_depth, sensor_max_depth, z_grid, used_profiles,sensor_std,error_grid)
        if status != 'okay':
            break
        else:
            sensor_depths_go_1000[rep]=sensor_interim
    

    if status == 'okay':
        sensor_depths_go_1000=np.array(sensor_depths_go_1000)
        sensor_depths_go_avg=np.mean(sensor_depths_go_1000,axis=0)
        accumulated_error_go_avg, errors_go_avg= do_sensor_analysis(borehole_simulations,sensor_depths_go_avg, sensor_std, z_grid, error_grid, used_profiles)

    else:
        go_status=status
        print(go_status)
        
if go_status=='okay':
    np.save('../output/Discussion/Advection_impact/GRIP_adv1000_sensor_depths2006_go_avg',sensor_depths_go_avg)
    np.save('../output/Discussion/Advection_impact/GRIP_adv1000_accumulated_error_2006_go_avg',accumulated_error_go_avg)
        