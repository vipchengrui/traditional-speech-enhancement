function [pesq_mos]= PESQ(ref_wav, deg_wav)

% ----------------------------------------------------------------------
%            PESQ objective speech quality measure
%
%   This function implements the PESQ measure based on the ITU standard
%   P.862 [1].
%
%
%   Usage:  pval=pesq(cleanFile.wav, enhancedFile.wav)
%           
%         cleanFile.wav - clean input file in .wav format
%         enhancedFile  - enhanced output file in .wav format
%         pval          - PESQ value
%
%    Note that the PESQ routine only supports sampling rates of 8 kHz and
%    16 kHz [1]
%
%  Example call:  pval = pesq ('sp04.wav','enhanced.wav')
%
%  
%  References:
%   [1] ITU (2000). Perceptual evaluation of speech quality (PESQ), and 
%       objective method for end-to-end speech quality assessment of 
%       narrowband telephone networks and speech codecs. ITU-T
%       Recommendation P. 862   
%
%   Authors: Yi Hu and Philipos C. Loizou 
%
%
% Copyright (c) 2006 by Philipos C. Loizou
% $Revision: 0.0 $  $Date: 10/09/2006 $
% ----------------------------------------------------------------------

if nargin<2
    fprintf('Usage: [pesq_mos]=pesq(cleanfile.wav,enhanced.wav) \n');
    return;
end;

global Downsample DATAPADDING_MSECS SEARCHBUFFER Fs WHOLE_SIGNAL
global Align_Nfft Window 

[ref_data,sampling_rate]= audioread(ref_wav);
if sampling_rate~=8000 & sampling_rate~=16000
    error('Sampling frequency needs to be either 8000 or 16000 Hz');
end

setup_global( sampling_rate);

% Window= hann( Align_Nfft, 'periodic'); %Hanning window
% Window= Window'; 
TWOPI= 6.28318530717959;
%for count = 0: Align_Nfft- 1
%    Window(1+ count) = 0.5 * (1.0 - cos((TWOPI * count) / Align_Nfft));
%end

count=0:Align_Nfft- 1;
Window= 0.5 * (1.0 - cos((TWOPI * count) / Align_Nfft));
  


ref_data= ref_data';
ref_data= ref_data* 32768;
ref_Nsamples= length( ref_data)+ 2* SEARCHBUFFER* Downsample;
ref_data= [zeros( 1, SEARCHBUFFER* Downsample), ref_data, ...
    zeros( 1, DATAPADDING_MSECS* (Fs/ 1000)+ SEARCHBUFFER* Downsample)];

deg_data= audioread( deg_wav);
deg_data= deg_data';
deg_data= deg_data* 32768;
deg_Nsamples= length( deg_data)+ 2* SEARCHBUFFER* Downsample;
deg_data= [zeros( 1, SEARCHBUFFER* Downsample), deg_data, ...
    zeros( 1, DATAPADDING_MSECS* (Fs/ 1000)+ SEARCHBUFFER* Downsample)];

maxNsamples= max( ref_Nsamples, deg_Nsamples);

ref_data= fix_power_level( ref_data, ref_Nsamples, maxNsamples);
deg_data= fix_power_level( deg_data, deg_Nsamples, maxNsamples);

standard_IRS_filter_dB= [0, -200; 50, -40; 100, -20; 125, -12; 160, -6; 200, 0;...    
    250, 4; 300, 6; 350, 8; 400, 10; 500, 11; 600, 12; 700, 12; 800, 12;...
    1000, 12; 1300, 12; 1600, 12; 2000, 12; 2500, 12; 3000, 12; 3250, 12;...
    3500, 4; 4000, -200; 5000, -200; 6300, -200; 8000, -200]; 

ref_data= apply_filter( ref_data, ref_Nsamples, standard_IRS_filter_dB);
deg_data= apply_filter( deg_data, deg_Nsamples, standard_IRS_filter_dB);
% 



% for later use in psychoacoustical model
model_ref= ref_data;
model_deg= deg_data;

[ref_data, deg_data]= input_filter( ref_data, ref_Nsamples, deg_data, ...
    deg_Nsamples);


[ref_VAD, ref_logVAD]= apply_VAD( ref_data, ref_Nsamples);
[deg_VAD, deg_logVAD]= apply_VAD( deg_data, deg_Nsamples);


crude_align (ref_logVAD, ref_Nsamples, deg_logVAD, deg_Nsamples,...
    WHOLE_SIGNAL);

utterance_locate (ref_data, ref_Nsamples, ref_VAD, ref_logVAD,...
    deg_data, deg_Nsamples, deg_VAD, deg_logVAD);

ref_data= model_ref;
deg_data= model_deg;

% make ref_data and deg_data equal length
if (ref_Nsamples< deg_Nsamples)
    newlen= deg_Nsamples+ DATAPADDING_MSECS* (Fs/ 1000);
    ref_data( newlen)= 0;
elseif (ref_Nsamples> deg_Nsamples)
    newlen= ref_Nsamples+ DATAPADDING_MSECS* (Fs/ 1000);
    deg_data( newlen)= 0;
end


pesq_mos= pesq_psychoacoustic_model (ref_data, ref_Nsamples, deg_data, ...
    deg_Nsamples );




