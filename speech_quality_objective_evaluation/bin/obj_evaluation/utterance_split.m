function utterance_split( ref_data, ref_Nsamples, ref_VAD, ref_logVAD, ...
    deg_data, deg_Nsamples, deg_VAD, deg_logVAD)

global Nutterances MAXNUTTERANCES Downsample SEARCHBUFFER
global Utt_DelayEst Utt_Delay Utt_DelayConf UttSearch_Start
global Utt_Start Utt_End Largest_uttsize UttSearch_End
global Best_ED1 Best_D1 Best_DC1 Best_ED2 Best_D2 Best_DC2 Best_BP

Utt_id = 1;
while( (Utt_id <= Nutterances) && (Nutterances <= MAXNUTTERANCES) )
    Utt_DelayEst_l = Utt_DelayEst(Utt_id);
    Utt_Delay_l = Utt_Delay(Utt_id);
    Utt_DelayConf_l = Utt_DelayConf(Utt_id);
    Utt_Start_l = Utt_Start(Utt_id);
    Utt_End_l = Utt_End(Utt_id);
    
    Utt_SpeechStart = Utt_Start_l;
%     fprintf( 'SpeechStart is %d\n', Utt_SpeechStart);
    while( (Utt_SpeechStart < Utt_End_l) && ...
            (ref_VAD(Utt_SpeechStart)<= 0.0) )
        Utt_SpeechStart = Utt_SpeechStart + 1;
    end %find the SpeechStart for each utterance
    Utt_SpeechEnd = Utt_End_l;
%     fprintf( 'SpeechEnd is %d\n', Utt_SpeechEnd);
    while( (Utt_SpeechEnd > Utt_Start_l) && ...
            (ref_VAD(Utt_SpeechEnd) <= 0))
        Utt_SpeechEnd = Utt_SpeechEnd- 1;
    end
    Utt_SpeechEnd = Utt_SpeechEnd+ 1;    
    %find SpeechEnd for each utterance
    Utt_Len = Utt_SpeechEnd - Utt_SpeechStart;
    
%     fprintf( 'Utt_Len is %d\n', Utt_Len);
    
    if( Utt_Len >= 200 )
        split_align( ref_data, ref_Nsamples, ref_VAD, ref_logVAD, ...
            deg_data, deg_Nsamples, deg_VAD, deg_logVAD, ...
            Utt_Start_l, Utt_SpeechStart, Utt_SpeechEnd, Utt_End_l, ...
            Utt_DelayEst_l, Utt_DelayConf_l);
%         fprintf( '\nBest_ED1, Best_D1, Best_DC1 is %d, %d, %f\n',...
% 				Best_ED1, Best_D1, Best_DC1);
%         fprintf( 'Best_ED2, Best_D2, Best_DC2 is %d, %d, %f\n',...
% 				Best_ED2, Best_D2, Best_DC2);
%         fprintf( 'Best_BP is %d\n', Best_BP);
                
        if( (Best_DC1 > Utt_DelayConf_l) && (Best_DC2 > Utt_DelayConf_l) )
            for step = Nutterances: -1: Utt_id+ 1
                Utt_DelayEst(step+ 1) = Utt_DelayEst(step);
                Utt_Delay(step+ 1) = Utt_Delay(step);
                Utt_DelayConf(step+ 1) = Utt_DelayConf(step);
                Utt_Start(step+ 1) = Utt_Start(step);
                Utt_End(step+ 1) = Utt_End(step);
                UttSearch_Start(step+ 1) = Utt_Start( step);
                UttSearch_End(step+ 1) = Utt_End( step);
            end

            Nutterances = Nutterances+ 1;
            
            Utt_DelayEst(Utt_id) = Best_ED1;
            Utt_Delay(Utt_id) = Best_D1;
            Utt_DelayConf(Utt_id) = Best_DC1;
            
            Utt_DelayEst(Utt_id +1) = Best_ED2;
            Utt_Delay(Utt_id +1) = Best_D2;
            Utt_DelayConf(Utt_id +1) = Best_DC2;
            
            UttSearch_Start(Utt_id +1) = UttSearch_Start(Utt_id);
            UttSearch_End(Utt_id +1) = UttSearch_End( Utt_id);
            if( Best_D2 < Best_D1 )
                Utt_Start(Utt_id) = Utt_Start_l;
                Utt_End(Utt_id) = Best_BP;
                Utt_Start(Utt_id +1) = Best_BP;
                Utt_End(Utt_id +1) = Utt_End_l;
            else
                Utt_Start( Utt_id) = Utt_Start_l;
                Utt_End( Utt_id) = Best_BP + ...
                    floor( (Best_D2- Best_D1)/ (2 * Downsample));
                Utt_Start( Utt_id +1) = Best_BP - ...
                    floor( (Best_D2- Best_D1)/ (2 * Downsample));
                Utt_End( Utt_id +1) = Utt_End_l;
            end

            if( (Utt_Start(Utt_id)- SEARCHBUFFER- 1)* Downsample+ 1+ ...
                    Best_D1 < 0 )
                Utt_Start(Utt_id) = SEARCHBUFFER+ 1+  ...
                    floor( (Downsample - 1 - Best_D1) / Downsample);
            end

            if( ((Utt_End( Utt_id +1)- 1)* Downsample+ 1 + Best_D2) >...
                    (deg_Nsamples - SEARCHBUFFER * Downsample) )
                Utt_End( Utt_id +1) = floor( (deg_Nsamples - Best_D2)...
                    / Downsample)- SEARCHBUFFER+ 1;
            end
        else
            Utt_id= Utt_id+ 1;
        end
    else
        Utt_id = Utt_id+ 1;
    end
end

Largest_uttsize = max( Utt_End- Utt_Start);

% fid= fopen( 'uttinfo_mat.txt', 'wt');
% fprintf( fid, 'Number of Utterances is:\n');
% fprintf( fid, '%d\n', Nutterances);
% fprintf( fid, 'Utterance Delay Estimation:\n');
% fprintf( fid, '%d\n', Utt_DelayEst( 1: Nutterances) );
% fprintf( fid, 'Utterance Delay:\n');
% fprintf( fid, '%d\n', Utt_Delay( 1: Nutterances));
% fprintf( fid, 'Utterance Delay Confidence:\n');
% fprintf( fid, '%f\n', Utt_DelayConf( 1: Nutterances));
% fprintf( fid, 'Utterance Start:\n');
% fprintf( fid, '%d\n', Utt_Start( 1: Nutterances));
% fprintf( fid, 'Utterance End:\n');
% fprintf( fid, '%d\n', Utt_End( 1: Nutterances));
% fprintf( fid, 'Largest utterance length:\n');
% fprintf( fid, '%d\n', Largest_uttsize);
% fclose( fid);



