import numpy as np
import dill
import random
from scipy.interpolate import CubicSpline
from scipy.stats import qmc

from sampling_error_calculation import do_sensor_analysis
from greedy_optimal_sampling import generate_optimized_sensor_locations_candiate_set




sensor_min_depth = 0
sensor_max_depth = 200
sensor_count = 20
sensor_std = 0

error_grid_size = 100000
error_grid = np.linspace(sensor_min_depth,sensor_max_depth, num=error_grid_size, endpoint=True)
error_grid_delta = error_grid[1]-error_grid[0]

random.seed(9001)
    
# array containing final temperature profiles in borehole
T_end2006 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA2_6_simulation_1000.npy')

# array containing surface temperature profiles
T_top2006 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA2_6_surfaceTem_1000.npy')

t_grid = np.load('../data_borehole_simulations/EDML_simulations/EDML_sample_times_fwd.npy')
z_grid = np.load('../data_borehole_simulations/EDML_simulations/EDML_sample_depths_fwd.npy')

N_z = z_grid.size
N_t = t_grid.size

sensor_depths2006_go_inc_profs={}
for i in range(1000):
    sensor_depths2006_go_inc_profs[i]=np.zeros((0,0))

for i in range(1000):
    random.seed(9001)
    sensor_depths_go_20=[]
    used_profiles = []
    idx = 0
    for k in range(1):
            used_profiles.append(idx)
            idx = idx + 1

    for n in range(20):
        sensor_depths_go_20.append(0)

    for rep in range(20):
        sensor_interim, status=generate_optimized_sensor_locations_candiate_set(np.reshape(T_end2006[i], (1,np.product(T_end2006[i].shape))).T , sensor_count, sensor_min_depth, sensor_max_depth, z_grid, used_profiles,sensor_std,error_grid)
        if status != 'okay':
            break
        else:
            sensor_depths_go_20[rep]=sensor_interim
    sensor_depths_go_20=np.array(sensor_depths_go_20)
    sensor_depths_go_avg=np.mean(sensor_depths_go_20,axis=0)
    accumulated_error_go_avg, errors_go_avg= do_sensor_analysis(np.reshape(T_end2006[i], (1,np.product(T_end2006[i].shape))).T,sensor_depths_go_avg, sensor_std, z_grid, error_grid, used_profiles)
    sensor_depths2006_go_inc_profs[i]=np.array(list(zip(sensor_depths_go_avg, accumulated_error_go_avg)))
    
    
    
with open('../output/Discussion/Past_climate_impact/EDML_sensor_depths2006_go_avg20_each_prof.pkl', 'wb') as f1:
    dill.dump(sensor_depths2006_go_inc_profs,f1)