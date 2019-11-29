function mod_data= fix_power_level( data, data_Nsamples, maxNsamples)
% this function is used for level normalization, i.e., to fix the power
% level of data to a preset number, and return it to mod_data. 

global Downsample DATAPADDING_MSECS SEARCHBUFFER Fs
global TARGET_AVG_POWER 
TARGET_AVG_POWER= 1e7;

align_filter_dB= [0,-500; 50, -500; 100, -500; 125, -500; 160, -500; 200, -500;
    250, -500; 300, -500; 350,  0; 400,  0; 500,  0; 600,  0; 630,  0;
    800,  0; 1000, 0; 1250, 0; 1600, 0; 2000, 0; 2500, 0; 3000, 0;
    3250, 0; 3500, -500; 4000, -500; 5000, -500; 6300, -500; 8000, -500];    

align_filtered= apply_filter( data, data_Nsamples, align_filter_dB);
power_above_300Hz = pow_of (align_filtered, SEARCHBUFFER* Downsample+ 1, ...
    data_Nsamples- SEARCHBUFFER* Downsample+ DATAPADDING_MSECS* (Fs/ 1000), ...
    maxNsamples- 2* SEARCHBUFFER* Downsample+ DATAPADDING_MSECS* (Fs/ 1000));

global_scale= sqrt( TARGET_AVG_POWER/ power_above_300Hz);
% fprintf( 1, '\tglobal_scale is %f\n', global_scale);
mod_data= data* global_scale;
