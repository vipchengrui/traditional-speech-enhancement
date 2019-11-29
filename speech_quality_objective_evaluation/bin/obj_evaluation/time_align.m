function time_align(ref_data, ref_Nsamples, ...
    deg_data, deg_Nsamples, Utt_id)

global Utt_DelayEst Utt_Delay Utt_DelayConf UttSearch_Start UttSearch_End 
global Align_Nfft Downsample Window

estdelay = Utt_DelayEst(Utt_id);

H = zeros( 1, Align_Nfft);
X1= zeros( 1, Align_Nfft);
X2= zeros( 1, Align_Nfft);

startr = (UttSearch_Start(Utt_id)- 1)* Downsample+ 1;
startd = startr + estdelay;
if ( startd < 0 )
    startr = 1 -estdelay;
    startd = 1;
end

while( ((startd + Align_Nfft) <= deg_Nsamples) && ...
        ((startr + Align_Nfft) <= ((UttSearch_End(Utt_id)- 1) * Downsample)) )
    X1= ref_data( startr: startr+ Align_Nfft- 1).* Window;    
    X2= deg_data( startd: startd+ Align_Nfft- 1).* Window;  
    
    % find cross-correlation between X1 and X2
    X1_fft= fft( X1, Align_Nfft );
    X1_fft_conj= conj( X1_fft);
    X2_fft= fft( X2, Align_Nfft );    
    X1= ifft( X1_fft_conj.* X2_fft, Align_Nfft );        

    X1= abs( X1);     
    v_max = max( X1)* 0.99;
    
    X1_greater_vmax= find( X1 > v_max );
    H( X1_greater_vmax )= H( X1_greater_vmax )+ v_max^ 0.125;
    
    startr = startr+ Align_Nfft/ 4;
    startd = startd+ Align_Nfft/ 4;

end

X1= H;
X2= 0;
Hsum = sum( H);

X2(1) = 1.0;
kernel = Align_Nfft / 64;

for count= 2: kernel
    X2( count)= 1- (count- 1)/ kernel;
    X2( Align_Nfft- count+ 2)= 1- (count- 1)/ kernel;
end
    
X1_fft= fft( X1, Align_Nfft );
X2_fft= fft( X2, Align_Nfft );

X1= ifft( X1_fft.* X2_fft, Align_Nfft );

if (Hsum> 0)
    H= abs( X1)/ Hsum;
else
    H= 0;
end

[v_max, I_max] = max( H);
if( I_max- 1 >= (Align_Nfft/2) )
    I_max = I_max- Align_Nfft;
end

Utt_Delay(Utt_id) = estdelay + I_max- 1;
Utt_DelayConf(Utt_id) = v_max; % confidence
    

    
    
    
