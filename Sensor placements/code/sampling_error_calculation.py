import numpy as np
from scipy.interpolate import CubicSpline

#Function to extract temperatures(sensor measurements) at sensor locations from the calculated borehole temperature profiles.
def get_measured_temperatures(T_end,sensor_depths, sensor_std, z_grid, used_profiles):
    measured_temperatures = []
    randListIndex= sensor_depths.shape[0]
    randList=np.load('RandomLists_device_error/RandomLists_%dsensor_randn_diff_90.npy'%(randListIndex))
    for p in used_profiles:
        spline = CubicSpline(z_grid, T_end[:,p])
        evaluated_sensor = spline(sensor_depths)+(sensor_std*randList[p,:])
        measured_temperatures.append(evaluated_sensor)
        
    return measured_temperatures


#Function to reconstruct continuous borehole temperature profiles from the sensor measurements using interpolation.
def build_interpolants(T_end,sensor_depths, measured_temperatures,used_profiles):
    interpolants = []
    for i in range(len(used_profiles)):
        interpolants.append(CubicSpline(sensor_depths, measured_temperatures[i]))
    return interpolants


#Function to calculate error between calculated borehole temperature profiles and reconstructed temperature profiles.
def evaluate_measurement_errors(interpolants, z_grid, error_grid, T_end, used_profiles):
    errors = []
    i = 0
    for p in used_profiles:
        T_end_interpolant = CubicSpline(z_grid, T_end[:,p])
        errors.append(np.abs(interpolants[i](error_grid) - T_end_interpolant(error_grid)))
        i = i + 1
    return errors

#Combine all steps in one "analysis" function. 
#This fuction return the sampling error and the individal errors w.r.t each borehole profile. 
def do_sensor_analysis(T_end,sensor_depths, sensor_std, z_grid, error_grid, used_profiles):
    measured_temperatures = get_measured_temperatures(T_end,sensor_depths, sensor_std, z_grid, used_profiles)
    interpolants = build_interpolants(T_end,sensor_depths, measured_temperatures, used_profiles)
    errors = evaluate_measurement_errors(interpolants, z_grid, error_grid, T_end, used_profiles)

    accumulated_error = np.zeros(errors[0].size)
    for i in range(len(used_profiles)):
        accumulated_error = accumulated_error + (errors[i]**2)

    accumulated_error = np.sqrt(accumulated_error / len(used_profiles))


    return accumulated_error,errors
