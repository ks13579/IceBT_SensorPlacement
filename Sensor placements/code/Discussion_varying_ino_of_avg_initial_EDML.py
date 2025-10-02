import numpy as np
import dill
import random
from scipy.interpolate import CubicSpline
from scipy.stats import qmc


from sampling_error_calculation import do_sensor_analysis
from greedy_optimal_sampling import generate_optimized_sensor_locations_candiate_set



no_avg=[1,20,100,250,500,1000]

for av in no_avg:
    
    random.seed(9001)
    
    # array containing final temperature profiles in borehole
    T_end1000 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA10_simulation_10.npy')
    T_end1006 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA1_6_simulation_10.npy')
    T_end1001 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA11_simulation_10.npy')
    T_end2000 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA20_simulation_10.npy')
    T_end2006 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA2_6_simulation_10.npy')
    T_end2001 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA21_simulation_10.npy')
    T_end3000 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA30_simulation_10.npy')
    T_end3006 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA3_6_simulation_10.npy')
    T_end3001 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA31_simulation_10.npy')


    # array containing surface temperature profiles
    T_top1000 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA10_surfaceTem_10.npy')
    T_top1006 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA1_6_surfaceTem_10.npy')
    T_top1001 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA11_surfaceTem_10.npy')
    T_top2000 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA20_surfaceTem_10.npy')
    T_top2006 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA2_6_surfaceTem_10.npy')
    T_top2001 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA21_surfaceTem_10.npy')
    T_top3000 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA30_surfaceTem_10.npy')
    T_top3006 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA3_6_surfaceTem_10.npy')
    T_top3001 = np.load('../data_borehole_simulations/EDML_simulations/EDML_PA31_surfaceTem_10.npy')

    t_grid = np.load('../data_borehole_simulations/EDML_simulations/EDML_sample_times_fwd.npy')
    z_grid = np.load('../data_borehole_simulations/EDML_simulations/EDML_sample_depths_fwd.npy')

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
    sensor_count = 20
    sensor_std = 0

    error_grid_size = 100000
    error_grid = np.linspace(sensor_min_depth,sensor_max_depth, num=error_grid_size, endpoint=True)
    error_grid_delta = error_grid[1]-error_grid[0]

    
       
    print('EDML_sensor_0derr_go_avges_%d'%(av))
    no_of_cases=av
    sensor_depths_go_1000=[]
    sensor_accumulated_error_go_1000=[] # for min line
    
    for n in range(no_of_cases):
        sensor_depths_go_1000.append(0)
        if no_of_cases == 1000:
            sensor_accumulated_error_go_1000.append(0)

    for rep in range(no_of_cases):        
        sensor_interim, status=generate_optimized_sensor_locations_candiate_set(borehole_simulations, sensor_count, sensor_min_depth, sensor_max_depth, z_grid, used_profiles,sensor_std,error_grid)
        if status != 'okay':
            break
        else:
            sensor_depths_go_1000[rep]=sensor_interim
  
    if status == 'okay':
        if no_of_cases == 1000:  # for min line
            for rep_ac in range(no_of_cases):
                accumulated_error_go_interim, errors_go_interim= do_sensor_analysis(borehole_simulations, sensor_depths_go_1000[rep_ac], sensor_std, z_grid, error_grid, used_profiles)
                sensor_accumulated_error_go_1000[rep_ac]=accumulated_error_go_interim
            sensor_accumulated_error_go_1000=np.array(sensor_accumulated_error_go_1000)
            np.save('../output/Discussion/Numerical_uncertainty/EDML_Model_error/EDML_%d_minLine_sensor_accumulated_error_go_1000'%(av),sensor_accumulated_error_go_1000)
            np.save('../output/Discussion/Numerical_uncertainty/EDML_Model_error/EDML_%d_minLine_sensors_go_1000'%(av),np.array(sensor_depths_go_1000))
              
        sensor_depths_go_1000=np.array(sensor_depths_go_1000)
        sensor_depths_go_avg=np.mean(sensor_depths_go_1000,axis=0)
        accumulated_error_go_avg, errors_go_avg= do_sensor_analysis(borehole_simulations,sensor_depths_go_avg, sensor_std, z_grid, error_grid, used_profiles)
            
        np.save('../output/Discussion/Numerical_uncertainty/EDML_Model_error/EDML_%d_sensor_depth_go_avges'%(av),sensor_depths_go_avg)
        np.save('../output/Discussion/Numerical_uncertainty/EDML_Model_error/EDML_%d_acc_error_go_avges'%(av),accumulated_error_go_avg)
        #np.save('../output/Discussion/Numerical_uncertainty/EDML_Model_error/EDML_%d_errors_go_avges'%(av),errors_go_avg)

    else:
        print(status)

 