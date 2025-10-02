import numpy as np
import random
from scipy.interpolate import CubicSpline
from scipy.stats import qmc

from sampling_error_calculation import do_sensor_analysis
from greedy_optimal_sampling import generate_optimized_sensor_locations_candiate_set



rm_sensors=[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

for sn in rm_sensors:
    
    go_status='okay'
    
    random.seed(9001)
    
    # array containing final temperature profiles in borehole
    T_end1000 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA10_simulation_10.npy')
    T_end1006 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA1_6_simulation_10.npy')
    T_end1001 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA11_simulation_10.npy')
    T_end2000 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA20_simulation_10.npy')
    T_end2006 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA2_6_simulation_10.npy')
    T_end2001 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA21_simulation_10.npy')
    T_end3000 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA30_simulation_10.npy')
    T_end3006 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA3_6_simulation_10.npy')
    T_end3001 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA31_simulation_10.npy')


    # array containing surface temperature profiles
#     T_top1000 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA10_surfaceTem_10.npy')
#     T_top1006 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA1_6_surfaceTem_10.npy')
#     T_top1001 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA11_surfaceTem_10.npy')
#     T_top2000 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA20_surfaceTem_10.npy')
#     T_top2006 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA2_6_surfaceTem_10.npy')
#     T_top2001 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA21_surfaceTem_10.npy')
#     T_top3000 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA30_surfaceTem_10.npy')
#     T_top3006 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA3_6_surfaceTem_10.npy')
#     T_top3001 = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_PA31_surfaceTem_10.npy')

    t_grid = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_sample_times_fwd.npy')
    z_grid = np.load('../data_borehole_simulations/GRIP_simulations/GRIP_sample_depths_fwd.npy')

    N_z = z_grid.size
    N_t = t_grid.size

    borehole_simulations=np.zeros((z_grid.size,90))

    borehole_simulations[:,0:10]=T_end1000.T
    borehole_simulations[:,10:20]=T_end1006.T
    borehole_simulations[:,20:30]=T_end1001.T

    borehole_simulations[:,30:40]=T_end2000.T
    borehole_simulations[:,40:50]=T_end2006.T
    borehole_simulations[:,50:60]=T_end2001.T

    borehole_simulations[:,60:70]=T_end3000.T
    borehole_simulations[:,70:80]=T_end3006.T
    borehole_simulations[:,80:90]=T_end3001.T

    used_profiles = []
    idx = 0
    for pa in [1,2,3]:
        for beta in [0,0.6,1]:
            for k in range(10):
                used_profiles.append(idx)
                idx = idx + 1
 

    #Select range in which sensors are placed and number of sensors

    sensor_min_depth = 0
    sensor_max_depth = 200
    sensor_count = sn#16 #12#8
    sensor_std = 0

    error_grid_size = 100000
    error_grid = np.linspace(sensor_min_depth,sensor_max_depth, num=error_grid_size, endpoint=True)
    #error_grid_delta = error_grid[1]-error_grid[0]

    #    Linear sensor placement    
    print('GRIP_res_%dsensor_0derr_linear'%(sn))
    sensor_depths_linear = np.linspace(sensor_min_depth,sensor_max_depth, num=sensor_count, endpoint=True)
    accumulated_error_linear, errors_linear = do_sensor_analysis(borehole_simulations,sensor_depths_linear, sensor_std, z_grid, error_grid, used_profiles)

    #    Exponential sensor placement    
    print('GRIP_res_%dsensor_0derr_logarthmic'%(sn))
    sensor_depths_log = np.geomspace(sensor_min_depth+1,sensor_max_depth+1, num=sensor_count, endpoint=True)
    sensor_depths_logarthmic = sensor_depths_log-1
    accumulated_error_logarthmic, errors_logarthmic = do_sensor_analysis(borehole_simulations,sensor_depths_logarthmic, sensor_std, z_grid, error_grid, used_profiles)

    #    Greedy optimal sensor placement
    if go_status=='okay':
        print('GRIP_res_%dsensor_0derr_go_avg'%(sn))
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

    # sensor_0derr normal sampling techniques with cubicSpline interpolation

    np.save('../output/GRIP_varying_sensors_5_20/device_error_00mK/GRIP_res_%dsensor_0derr_sensor_depths_linear'%(sn),sensor_depths_linear)
    np.save('../output/GRIP_varying_sensors_5_20/device_error_00mK/GRIP_res_%dsensor_0derr_accumulated_error_linear'%(sn),accumulated_error_linear)
    np.save('../output/GRIP_varying_sensors_5_20/device_error_00mK/GRIP_res_%dsensor_0derr_sensor_depths_logarthmic'%(sn),sensor_depths_logarthmic)
    np.save('../output/GRIP_varying_sensors_5_20/device_error_00mK/GRIP_res_%dsensor_0derr_accumulated_error_logarthmic'%(sn),accumulated_error_logarthmic)
    if go_status == 'okay':
        np.save('../output/GRIP_varying_sensors_5_20/device_error_00mK/GRIP_res_%dsensor_0derr_sensor_depths_go_avg'%(sn),sensor_depths_go_avg)
        np.save('../output/GRIP_varying_sensors_5_20/device_error_00mK/GRIP_res_%dsensor_0derr_accumulated_error_go_avg'%(sn),accumulated_error_go_avg)
