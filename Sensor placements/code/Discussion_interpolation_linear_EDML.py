import matplotlib.pyplot as plt
import numpy as np
import dill
import random

from scipy import interpolate

#Function to extract temperatures at sensor locations from the calculated temperature profiles
def get_measured_temperatures(T_end,sensor_depths, sensor_std, z_grid, used_profiles):
    measured_temperatures = []
    randListIndex= sensor_depths.shape[0]
    randList=np.load('RandomLists_device_error/RandomLists_%dsensor_randn_diff_90.npy'%(randListIndex))
    for p in used_profiles:
        
        spline = interpolate.interp1d(z_grid, T_end[:,p],kind='linear')
        evaluated_sensor = spline(sensor_depths)+(sensor_std*randList[p,:])
        measured_temperatures.append(evaluated_sensor)
        
    return measured_temperatures


#Function to reconstruct from the sensor measurements continuous temperature profiles using interpolation.
def build_interpolants(T_end,sensor_depths, measured_temperatures,used_profiles):
    interpolants = []
    for i in range(len(used_profiles)):
        interpolants.append(interpolate.interp1d(sensor_depths, measured_temperatures[i],kind='linear'))
    return interpolants


#Function to calculate error between calculated temperature profiles and reconstructed temperature profiles
def evaluate_measurement_errors(interpolants, z_grid, error_grid, T_end, used_profiles):
    errors = []
    i = 0
    for p in used_profiles:
        T_end_interpolant = interpolate.interp1d(z_grid, T_end[:,p],kind='linear')
        errors.append(np.abs(interpolants[i](error_grid) - T_end_interpolant(error_grid)))
        i = i + 1
    return errors

#Combine all steps in one "analysis" function
def do_sensor_analysis(T_end,sensor_depths, sensor_std, z_grid, error_grid, used_profiles):
    measured_temperatures = get_measured_temperatures(T_end,sensor_depths, sensor_std, z_grid, used_profiles)
    interpolants = build_interpolants(T_end,sensor_depths, measured_temperatures, used_profiles)
    errors = evaluate_measurement_errors(interpolants, z_grid, error_grid, T_end, used_profiles)

    accumulated_error = np.zeros(errors[0].size)
    for i in range(len(used_profiles)):
        accumulated_error = accumulated_error + (errors[i]**2)

    accumulated_error = np.sqrt(accumulated_error / len(used_profiles))


    return accumulated_error,errors

#Function to generate a sensor placement that minimizes the variance in the error by sequentially adding sensors at locations of maximum error variance.
def get_inner_sensors(sensor_min_depth,sensor_max_depth):
    
    import numpy as np
    from scipy.stats import qmc
    
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

def generate_optimized_sensor_locations_candiate_set_1(T_end, sensor_count, sensor_min_depth, sensor_max_depth, z_grid, used_profiles):

    import scipy as sp
    from scipy.stats import rv_continuous,sampling,qmc

    status='okay'
    #    # generate candidates set with candidates being all grid points in z direction (dirty hack, since analysis function
    #    # currently only returns the errors at the grid points)
    #    candidate_set = z_grid

    initial_sensor_count = 4
    print('initial_sensor_count',initial_sensor_count)
    
    # set initial sensor set with just 5 sensors as starting point
    sensors = np.linspace(sensor_min_depth,sensor_max_depth, num=initial_sensor_count, endpoint=True)
    sensor_select=get_inner_sensors(sensor_min_depth,sensor_max_depth)
    sensors[1]=min(sensor_select)
    sensors[2]=max(sensor_select)
    
    print(sensors)
     
    for i in range(sensor_count - initial_sensor_count):
        #    # make sure that initially chosen sensor locations are no longer in candidate set
        #    candidate_set = np.setdiff1d(candidate_set, sensors)
    
        # run analysis with current choice
        sensors = np.sort(sensors)
        accumulated_error, errors = do_sensor_analysis(T_end,sensors, sensor_std, z_grid, error_grid, used_profiles)

        # plt.figure()
        # plt.plot(accumulated_error,error_grid)
        # plt.scatter(np.zeros(sensors.size),sensors,marker='x',c='red')
        # plt.ticklabel_format(scilimits=[-3,6])
        # plt.xlabel("Total error over all calculated temperature profiles")
        # plt.ylabel("Depth")
        # plt.gca().invert_yaxis()
        
        # find indices of grid points with maximum errors
        max_indices = np.argmax(accumulated_error)

        # add new sensor and remove it from candidate set
        sensors = np.append(sensors, [error_grid[max_indices]])
        sensors = np.sort(sensors)
        sn_values, counts = np.unique(sensors, return_counts=True)
        if any(counts>1):
            status='duplicate_%f_for_%d_sensors'%(error_grid[max_indices],sensor_count)
            break
   
    return sensors,status

rm_sensors=[20]

for sn in rm_sensors:
    
    go_status='okay'
    
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
    sensor_count = sn#16 #12#8
    sensor_std = 0

    error_grid_size = 100000
    error_grid = np.linspace(sensor_min_depth,sensor_max_depth, num=error_grid_size, endpoint=True)
    error_grid_delta = error_grid[1]-error_grid[0]

    
    print('EDML_res_%dsensor_0derr_linear'%(sn))
    sensor_depths_linear = np.linspace(sensor_min_depth,sensor_max_depth, num=sensor_count, endpoint=True)
    accumulated_error_linear, errors_linear = do_sensor_analysis(borehole_simulations,sensor_depths_linear, sensor_std, z_grid, error_grid, used_profiles)

    print('EDML_res_%dsensor_0derr_logarthmic'%(sn))
    sensor_depths_log = np.geomspace(sensor_min_depth+1,sensor_max_depth+1, num=sensor_count, endpoint=True)
    sensor_depths_logarthmic = sensor_depths_log-1
    accumulated_error_logarthmic, errors_logarthmic = do_sensor_analysis(borehole_simulations,sensor_depths_logarthmic, sensor_std, z_grid, error_grid, used_profiles)

    if go_status=='okay':
        print('EDML_res_%dsensor_0derr_go_avg'%(sn))
        no_of_cases=1000
        sensor_depths_go_1000=[]
    
        for n in range(no_of_cases):
            sensor_depths_go_1000.append(0)

        for rep in range(no_of_cases):        
            sensor_interim, status=generate_optimized_sensor_locations_candiate_set_1(borehole_simulations, sensor_count, sensor_min_depth, sensor_max_depth, z_grid, used_profiles)
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

    np.save('../output/Discussion/Interpolation_types/EDML_linear_interp_%dsensor_depths_linear'%(sn),sensor_depths_linear)
    np.save('../output/Discussion/Interpolation_types/EDML_linear_interp_%dsensor_0derr_accumulated_error_linear'%(sn),accumulated_error_linear)
    np.save('../output/Discussion/Interpolation_types/EDML_linear_interp_%dsensor_0derr_sensor_depths_logarthmic'%(sn),sensor_depths_logarthmic)
    np.save('../output/Discussion/Interpolation_types/EDML_linear_interp_%dsensor_0derr_accumulated_error_logarthmic'%(sn),accumulated_error_logarthmic)
    
    if go_status == 'okay':
        np.save('../output/Discussion/Interpolation_types/EDML_linear_interp_%dsensor_0derr_sensor_depths_go_avg'%(sn),sensor_depths_go_avg)
        np.save('../output/Discussion/Interpolation_types/EDML_linear_interp_%dsensor_0derr_accumulated_error_go_avg'%(sn),accumulated_error_go_avg)
