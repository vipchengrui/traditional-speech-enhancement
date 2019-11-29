function split_align( ref_data, ref_Nsamples, ref_VAD, ref_logVAD, ...
    deg_data, deg_Nsamples, deg_VAD, deg_logVAD, ...
    Utt_Start_l, Utt_SpeechStart, Utt_SpeechEnd, Utt_End_l, ...
    Utt_DelayEst_l, Utt_DelayConf_l)

global MAXNUTTERANCES Align_Nfft Downsample Window    
global Utt_DelayEst Utt_Delay UttSearch_Start UttSearch_End 
global Best_ED1 Best_D1 Best_DC1 Best_ED2 Best_D2 Best_DC2 Best_BP

Utt_BPs= zeros( 1, 41);
Utt_ED1= zeros( 1, 41);
Utt_ED2= zeros( 1, 41);
Utt_D1= zeros( 1, 41);
Utt_D2= zeros( 1, 41);
Utt_DC1= zeros( 1, 41);
Utt_DC2= zeros( 1, 41);


Utt_Len = Utt_SpeechEnd - Utt_SpeechStart;
Utt_Test = MAXNUTTERANCES;
Best_DC1 = 0.0;
Best_DC2 = 0.0;
kernel = Align_Nfft / 64;
Delta = Align_Nfft / (4 * Downsample);
Step = floor( ((0.801 * Utt_Len + 40 * Delta - 1)/(40 * Delta)));
Step = Step* Delta;
% fprintf( 'Step is %f\n', Step);

Pad = floor( Utt_Len / 10);
if( Pad < 75 ) 
    Pad = 75;
end

Utt_BPs(1) = Utt_SpeechStart + Pad;
N_BPs = 1;
while( 1)
    N_BPs= N_BPs+ 1;
    Utt_BPs(N_BPs)= Utt_BPs(N_BPs- 1)+ Step;
    if (~((Utt_BPs(N_BPs) <= (Utt_SpeechEnd- Pad)) && (N_BPs <= 40) ))
        break;
    end
end

if( N_BPs <= 1 ) 
    return;
end

% fprintf( 'Utt_DelayEst_l, Utt_Start_l, N_BPs is %d,%d,%d\n', ...
%     Utt_DelayEst_l, Utt_Start_l, N_BPs);
for bp = 1: N_BPs- 1
    Utt_DelayEst(Utt_Test) = Utt_DelayEst_l;
    UttSearch_Start(Utt_Test) = Utt_Start_l;
    UttSearch_End(Utt_Test) = Utt_BPs(bp);
%     fprintf( 'bp,Utt_BPs(%d) is %d,%d\n', bp,bp,Utt_BPs(bp)); 
    
    crude_align( ref_logVAD, ref_Nsamples, deg_logVAD, ...
        deg_Nsamples, MAXNUTTERANCES);
    Utt_ED1(bp) = Utt_Delay(Utt_Test);

    Utt_DelayEst(Utt_Test) = Utt_DelayEst_l;
    UttSearch_Start(Utt_Test) = Utt_BPs(bp);
    UttSearch_End(Utt_Test) = Utt_End_l;
    
    crude_align( ref_logVAD, ref_Nsamples, deg_logVAD, ...
        deg_Nsamples, MAXNUTTERANCES);
    Utt_ED2(bp) = Utt_Delay(Utt_Test);
end

% stream = fopen( 'matmat.txt', 'wt' );	
% for count= 1: N_BPs- 1 
%     fprintf( stream, '%d\n', Utt_ED2(count));
% end
% fclose( stream );


Utt_DC1(1: N_BPs-1) = -2.0;
% stream= fopen( 'what_mmm.txt', 'at');
while( 1 )
    bp = 1;
    while( (bp <= N_BPs- 1) && (Utt_DC1(bp) > -2.0) )
        bp = bp+ 1;
    end
    if( bp >= N_BPs )
        break;
    end
    
    estdelay = Utt_ED1(bp);
%     fprintf( 'bp,estdelay is %d,%d\n', bp, estdelay);
    H(1: Align_Nfft)= 0;
    Hsum = 0.0;
    
    startr = (Utt_Start_l- 1) * Downsample+ 1;
    startd = startr + estdelay;
%     fprintf( 'startr/startd is %d/%d\n', startr, startd);
    
    if ( startd < 0 )
        startr = -estdelay+ 1;
        startd = 1;
    end

    while( ((startd + Align_Nfft) <= 1+ deg_Nsamples) &&...
            ((startr + Align_Nfft) <= (1+ (Utt_BPs(bp)- 1) * Downsample)) )
        X1= ref_data(startr: startr+ Align_Nfft- 1).* Window;
        X2= deg_data(startd: startd+ Align_Nfft- 1).* Window;
        
        X1_fft= fft( X1, Align_Nfft );
        X1_fft_conj= conj( X1_fft);
        X2_fft= fft( X2, Align_Nfft );
        X1= ifft( X1_fft_conj.* X2_fft, Align_Nfft);
        
        X1= abs( X1);
        v_max= max( X1)* 0.99;        
        n_max = (v_max^ 0.125 )/ kernel;
%         fprintf( stream, '%f %f\n', v_max, n_max);
        
        for count = 0: Align_Nfft- 1
            if( X1(count+ 1) > v_max )
                Hsum = Hsum+ n_max * kernel;
                for k = 1-kernel: kernel- 1
                    H(1+ rem( count+ k+ Align_Nfft, Align_Nfft))= ...
                        H(1+ rem(count+ k+ Align_Nfft, Align_Nfft))+ ...
                        n_max* (kernel- abs(k));
                end
            end
        end

        startr = startr+ (Align_Nfft / 4);
        startd = startd+ (Align_Nfft / 4);
    end

    [v_max, I_max] = max( H);
    if( I_max- 1 >= (Align_Nfft/2) )
        I_max = I_max- Align_Nfft;
    end

    Utt_D1(bp) = estdelay + I_max- 1;
    if( Hsum > 0.0 )
%         if (Utt_Len== 236)
%             fprintf( 'v_max, Hsum is %f, %f\n', v_max, Hsum);
%         end
        Utt_DC1(bp) = v_max / Hsum;
    else
        Utt_DC1(bp) = 0.0;
    end

%     fprintf( 'bp/startr/startd is %d/%d/%d\n', bp, startr, startd);
    while( bp < (N_BPs - 1) )
        bp = bp + 1;
        
        if( (Utt_ED1(bp) == estdelay) && (Utt_DC1(bp) <= -2.0) )
%             loopno= 0;
            while(((startd+ Align_Nfft)<= 1+ deg_Nsamples) && ...
                    ((startr+ Align_Nfft)<= ...
                    ((Utt_BPs(bp)- 1)* Downsample+ 1) ))
                X1= ref_data( startr: startr+ Align_Nfft- 1).* ...
                    Window;
% %                 if (Utt_Len== 321)
%                     fid= fopen( 'what_mat.txt', 'at');
%                     fprintf( fid, '%f\n', Window);
%                     fclose( fid);
% %                     fprintf( '\n');
% %                 end
                X2= deg_data( startd: startd+ Align_Nfft- 1).* ...
                    Window;
                X1_fft= fft( X1, Align_Nfft );
                X1_fft_conj= conj( X1_fft);
                X2_fft= fft( X2, Align_Nfft );
                X1= ifft( X1_fft_conj.* X2_fft, Align_Nfft);
                
                X1= abs( X1);
                v_max = 0.99* max( X1);
                n_max = (v_max^ 0.125)/ kernel;
%                 fprintf( 'v_max n_max is %f %f\n', v_max, n_max);
                
                for count = 0: Align_Nfft- 1
                    if( X1(count+ 1) > v_max )
                        Hsum = Hsum+ n_max * kernel;
                        for k = 1-kernel: kernel-1
                            H(1+ rem( count+ k+ Align_Nfft, Align_Nfft))= ...
                                H(1+ rem(count+ k+ Align_Nfft, Align_Nfft))+ ...
                                n_max* (kernel- abs(k));
                        end
                    end
                end

                startr = startr+ (Align_Nfft / 4);
                startd = startd+ (Align_Nfft / 4);
                
%                 loopno= loopno+ 1;
            end
%             fprintf( 'loopno is %d\n', loopno);

            [v_max, I_max] = max( H);
%             fprintf( 'I_max is %d ', I_max);
            if( I_max- 1 >= (Align_Nfft/2) )
                I_max = I_max- Align_Nfft;
            end
            

            Utt_D1(bp) = estdelay + I_max- 1;
            if( Hsum > 0.0 )
%                 fprintf( 'v_max Hsum is %f %f\n', v_max, Hsum);
                Utt_DC1(bp) = v_max / Hsum;
            else
                Utt_DC1(bp) = 0.0;
            end
        end
    end
end
% fclose( stream);

for bp= 1: N_BPs- 1
    if( Utt_DC1(bp) > Utt_DelayConf_l )
        Utt_DC2(bp) = -2.0;
    else
        Utt_DC2(bp) = 0.0;
    end
end

while( 1 )
    bp = N_BPs- 1;
    while( (bp >= 1) && (Utt_DC2(bp) > -2.0) )
        bp = bp- 1; 
    end
    if( bp < 1 )
        break;
    end 

    estdelay = Utt_ED2(bp);
    H( 1: Align_Nfft)= 0;
    Hsum = 0.0;
    
    startr = (Utt_End_l- 1)* Downsample+ 1- Align_Nfft;
    startd = startr + estdelay;
    
%     fprintf( '***NEW startr is %d\n', startr);
    
%     fprintf( 'startr/d, deg_Nsamples is %d/%d, %d\n', startr,startd, ...
%         deg_Nsamples);
%     fprintf( 'deg_data has %d elements\n', numel( deg_data));
    
    if ( (startd + Align_Nfft) > deg_Nsamples+ 1 )
        startd = deg_Nsamples - Align_Nfft+ 1;
        startr = startd - estdelay;
    end

    while( (startd>= 1) && (startr>= (Utt_BPs(bp)- 1)* Downsample+ 1) )
        X1= ref_data( startr: startr+ Align_Nfft- 1).* Window;
        X2= deg_data( startd: startd+ Align_Nfft- 1).* Window;
        
        X1_fft= fft( X1, Align_Nfft);
        X1_fft_conj= conj( X1_fft);
        X2_fft= fft( X2, Align_Nfft);
        
        X1= ifft( X1_fft_conj.* X2_fft, Align_Nfft );
        X1= abs( X1);
        
        v_max = max( X1)* 0.99;
        n_max = ( v_max^ 0.125 )/ kernel;
        
        for count = 0: Align_Nfft- 1
            if( X1(count+ 1) > v_max )
                Hsum = Hsum+ n_max * kernel;
                for k = 1-kernel: kernel- 1
                    H(1+ rem(count+ k+ Align_Nfft, Align_Nfft))= ...
                        H(1+ rem(count+ k+ Align_Nfft, Align_Nfft))+ ...
                        n_max* (kernel- abs(k));
                end
            end
        end

        startr = startr- (Align_Nfft / 4);
        startd = startd- (Align_Nfft / 4);
    end

    [v_max, I_max] = max( H);
    if( I_max- 1 >= (Align_Nfft/2) )
        I_max = I_max- Align_Nfft;
    end

    Utt_D2(bp) = estdelay + I_max- 1;
    if( Hsum > 0.0 )
        Utt_DC2(bp) = v_max / Hsum;
    else
        Utt_DC2(bp) = 0.0;
    end

    while( bp > 1 )
        bp = bp - 1;
        if( (Utt_ED2(bp) == estdelay) && (Utt_DC2(bp) <= -2.0) )
            while( (startd >= 1) && (startr >= (Utt_BPs(bp)- 1) * Downsample+ 1)) 
                 X1= ref_data( startr: startr+ Align_Nfft- 1).* Window;
                 X2= deg_data( startd: startd+ Align_Nfft- 1).* Window;
                 X1_fft_conj= conj( fft( X1, Align_Nfft));
                 X2_fft= fft( X2, Align_Nfft);
                 X1= ifft( X1_fft_conj.* X2_fft, Align_Nfft);
                 
                 X1= abs( X1);
                 v_max = max( X1)* 0.99;
                 n_max = (v_max^ 0.125)/ kernel;
                 
                 for count = 0: Align_Nfft- 1
                     if( X1(count+ 1) > v_max )
                         Hsum = Hsum+ n_max * kernel;
                         for k = 1-kernel: kernel- 1
                             H(1+ rem( count+ k+ Align_Nfft, Align_Nfft))= ...
                                 H(1+ rem(count+ k+ Align_Nfft, Align_Nfft))+ ...
                                 n_max* (kernel- abs(k));
                         end
                     end
                 end

                 startr = startr- (Align_Nfft / 4);
                 startd = startd- (Align_Nfft / 4);
            end

            [v_max, I_max] = max( H);
            if( I_max- 1 >= (Align_Nfft/2) )
                I_max = I_max- Align_Nfft;
            end
            

            Utt_D2(bp) = estdelay + I_max- 1;
            if( Hsum > 0.0 )
                Utt_DC2(bp) = v_max / Hsum;
            else
                Utt_DC2(bp) = 0.0;
            end
        end
    end
end

% fid= fopen( 'uttinfo_mat.txt', 'wt');
% fprintf( fid, '%f\n', Utt_D2);
% fprintf( fid, '\n');
% fprintf( fid, '%f\n', Utt_DC2);
% fclose( fid);

% fprintf( 'Utt_Len, N_BPs is %d, %d\n', Utt_Len, N_BPs);
for bp = 1: N_BPs- 1
    if( (abs(Utt_D2(bp) - Utt_D1(bp)) >= Downsample) && ...
            ((Utt_DC1(bp)+ Utt_DC2(bp))> (Best_DC1 + Best_DC2)) &&...
            (Utt_DC1(bp) > Utt_DelayConf_l) && ...
            (Utt_DC2(bp) > Utt_DelayConf_l) )
        Best_ED1 = Utt_ED1(bp);
        Best_D1 = Utt_D1(bp);
        Best_DC1 = Utt_DC1(bp);
        Best_ED2 = Utt_ED2(bp);
        Best_D2 = Utt_D2(bp);
        Best_DC2 = Utt_DC2(bp);
        Best_BP = Utt_BPs(bp);
%         fprintf( 'in loop...');
    end
end

% if (Utt_Len== 236)
%     fid= fopen( 'matmat.txt', 'wt');
%     fprintf( fid, 'N_BPs is %d\n', N_BPs);
%     fprintf( fid, 'Utt_DelayConf is %f\n', Utt_DelayConf_l);
%     fprintf( fid, 'ED2\t ED1\t D2\t D1\t DC2\t DC1\t BPs\n');
%     for bp= 1: N_BPs- 1
%         fprintf( fid, '%d\t %d\t %d\t %d\t %f\t %f\t %d\n', Utt_ED2( bp), ...
%             Utt_ED1( bp), Utt_D2(bp), Utt_D1(bp), Utt_DC2(bp),...
%             Utt_DC1( bp), Utt_BPs( bp));
%     end
%     fclose( fid);
% end























