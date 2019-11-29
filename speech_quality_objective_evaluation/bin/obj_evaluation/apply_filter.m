function align_filtered= apply_filter( data, data_Nsamples, align_filter_dB)

global Downsample DATAPADDING_MSECS SEARCHBUFFER Fs

align_filtered= data;
n= data_Nsamples- 2* SEARCHBUFFER* Downsample+ DATAPADDING_MSECS* (Fs/ 1000);
% now find the next power of 2 which is greater or equal to n
pow_of_2= 2^ (ceil( log2( n)));

[number_of_points, trivial]= size( align_filter_dB);
overallGainFilter= interp1( align_filter_dB( :, 1), align_filter_dB( :, 2), ...
    1000);

x= zeros( 1, pow_of_2);
x( 1: n)= data( SEARCHBUFFER* Downsample+ 1: SEARCHBUFFER* Downsample+ n);

x_fft= fft( x, pow_of_2);

freq_resolution= Fs/ pow_of_2;

factorDb( 1: pow_of_2/2+ 1)= interp1( align_filter_dB( :, 1), ...
    align_filter_dB( :, 2), (0: pow_of_2/2)* freq_resolution)- ...
    overallGainFilter;
factor= 10.^ (factorDb/ 20);

factor= [factor, fliplr( factor( 2: pow_of_2/2))];
x_fft= x_fft.* factor;

y= ifft( x_fft, pow_of_2);

align_filtered( SEARCHBUFFER* Downsample+ 1: SEARCHBUFFER* Downsample+ n)...
    = y( 1: n);

% fid= fopen( 'log_mat.txt', 'wt');
% fprintf( fid, '%f\n', y( 1: n));
% fclose( fid);




