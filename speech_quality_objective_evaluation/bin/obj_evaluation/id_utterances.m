function id_utterances( ref_Nsamples, ref_VAD, deg_Nsamples)

global Largest_uttsize MINUTTLENGTH MINUTTLENGTH Crude_DelayEst
global Downsample SEARCHBUFFER Nutterances Utt_Start
global Utt_End Utt_Delay

Utt_num = 1;
speech_flag = 0;
VAD_length = floor( ref_Nsamples / Downsample);
% fprintf( 1, 'VAD_length is %d\n', VAD_length);

del_deg_start = MINUTTLENGTH - Crude_DelayEst / Downsample;
del_deg_end = floor((deg_Nsamples- Crude_DelayEst)/ Downsample) ...
    - MINUTTLENGTH;

for count = 1: VAD_length 
    VAD_value = ref_VAD(count);
    if( (VAD_value > 0.0) && (speech_flag == 0) ) 
        speech_flag = 1;
        this_start = count;
        Utt_Start (Utt_num) = count;
    end

    if( ((VAD_value == 0) || (count == VAD_length)) && ...
            (speech_flag == 1) ) 
        speech_flag = 0;
        Utt_End (Utt_num) = count;
        
        if( ((count - this_start) >= MINUTTLENGTH) && ...
                (this_start < del_deg_end) && ... 
                (count > del_deg_start) )
            Utt_num = Utt_num + 1;   
        end
    end
end

Utt_Start(1) = SEARCHBUFFER+ 1;
Utt_End(Nutterances) = VAD_length - SEARCHBUFFER+ 1;

for Utt_num = 2: Nutterances
    this_start = Utt_Start(Utt_num)- 1;
    last_end = Utt_End(Utt_num - 1)- 1;
    count = floor( (this_start + last_end) / 2);
    Utt_Start(Utt_num) = count+ 1;
    Utt_End(Utt_num - 1) = count+ 1;
end

this_start = (Utt_Start(1)- 1) * Downsample + Utt_Delay(1);
if( this_start < (SEARCHBUFFER * Downsample) )
    count = SEARCHBUFFER + floor( ...
        (Downsample - 1 - Utt_Delay(1)) / Downsample);
    Utt_Start(1) = count+ 1;
end

last_end = (Utt_End(Nutterances)- 1) * Downsample + 1 + ...
    Utt_Delay(Nutterances);
% fprintf( 'Utt_End(%d) is %d\n', Nutterances, Utt_End(Nutterances));
% fprintf( 'last_end is %d\n', last_end);
% fprintf( 'Utt_Delay(%d) is %d\n', Nutterances, Utt_Delay(Nutterances));
if( last_end > (deg_Nsamples - SEARCHBUFFER * Downsample+ 1) )
    count = floor( (deg_Nsamples - Utt_Delay(Nutterances)) / Downsample) ...
        - SEARCHBUFFER;
    Utt_End(Nutterances) = count+ 1;
end

for Utt_num = 2: Nutterances
    this_start = (Utt_Start(Utt_num)- 1) * Downsample + Utt_Delay(Utt_num);
    last_end = (Utt_End(Utt_num - 1)- 1) * Downsample + Utt_Delay(Utt_num - 1);
    if( this_start < last_end )
        count = floor( (this_start + last_end) / 2);
        this_start = floor( (Downsample- 1+ count- Utt_Delay(Utt_num))...
            / Downsample);
        last_end = floor( (count - Utt_Delay(Utt_num - 1))...
            / Downsample);
        Utt_Start(Utt_num) = this_start+ 1;
        Utt_End(Utt_num- 1) = last_end+ 1;
    end
end

Largest_uttsize= max( Utt_End- Utt_Start);    
    
    
    
    
  