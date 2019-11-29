function id_searchwindows( ref_VAD, ref_Nsamples, deg_VAD, deg_Nsamples);

global MINUTTLENGTH Downsample MINUTTLENGTH SEARCHBUFFER
global Crude_DelayEst Nutterances UttSearch_Start UttSearch_End

Utt_num = 1;
speech_flag = 0;

VAD_length= floor( ref_Nsamples/ Downsample);
del_deg_start= MINUTTLENGTH- Crude_DelayEst/ Downsample;
del_deg_end= floor((deg_Nsamples- Crude_DelayEst)/ Downsample)-...
    MINUTTLENGTH;

for count= 1: VAD_length
    VAD_value= ref_VAD(count);
    if( (VAD_value> 0) && (speech_flag== 0) ) 
        speech_flag= 1;
        this_start= count;
        UttSearch_Start(Utt_num)= count- SEARCHBUFFER;
        if( UttSearch_Start(Utt_num)< 0 )
            UttSearch_Start(Utt_num)= 0;
        end
    end

    if( ((VAD_value== 0) || (count == (VAD_length-1))) && ...
            (speech_flag == 1) ) 
        speech_flag = 0;
        UttSearch_End(Utt_num) = count + SEARCHBUFFER;
        if( UttSearch_End(Utt_num) > VAD_length - 1 )
            UttSearch_End(Utt_num) = VAD_length -1;
        end

        if( ((count - this_start) >= MINUTTLENGTH) &&...
                (this_start < del_deg_end) &&...
                (count > del_deg_start) )
            Utt_num= Utt_num + 1;            
        end
    end
end
Utt_num= Utt_num- 1;
Nutterances = Utt_num;
    
% fprintf( 1, 'Nutterances is %d\n', Nutterances);

% fid= fopen( 'mat_utt.txt', 'wt');
% fprintf( fid, '%d\n', UttSearch_Start( 1: Nutterances));
% fprintf( fid, '\n');
% fprintf( fid, '%d\n', UttSearch_End( 1: Nutterances));
% fclose(fid);



















