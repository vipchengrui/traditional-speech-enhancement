function pesq_mos= pesq_psychoacoustic_model (ref_data, ref_Nsamples, deg_data, ...
    deg_Nsamples )

global CALIBRATE Nfmax Nb Sl Sp
global nr_of_hz_bands_per_bark_band centre_of_band_bark
global width_of_band_hz centre_of_band_hz width_of_band_bark
global pow_dens_correction_factor abs_thresh_power
global Downsample SEARCHBUFFER DATAPADDING_MSECS Fs Nutterances
global Utt_Start Utt_End Utt_Delay NUMBER_OF_PSQM_FRAMES_PER_SYLLABE 
global Fs Plot_Frame

% Plot_Frame= 75; % this is the frame whose spectrum will be plotted

FALSE= 0;
TRUE= 1;
NUMBER_OF_PSQM_FRAMES_PER_SYLLABE= 20;

maxNsamples = max (ref_Nsamples, deg_Nsamples);
Nf = Downsample * 8;
MAX_NUMBER_OF_BAD_INTERVALS = 1000;

start_frame_of_bad_interval= zeros( 1, MAX_NUMBER_OF_BAD_INTERVALS);
stop_frame_of_bad_interval= zeros( 1, MAX_NUMBER_OF_BAD_INTERVALS);
start_sample_of_bad_interval= zeros( 1, MAX_NUMBER_OF_BAD_INTERVALS);
stop_sample_of_bad_interval= zeros( 1, MAX_NUMBER_OF_BAD_INTERVALS);
number_of_samples_in_bad_interval= zeros( 1, MAX_NUMBER_OF_BAD_INTERVALS);
delay_in_samples_in_bad_interval= zeros( 1, MAX_NUMBER_OF_BAD_INTERVALS);
number_of_bad_intervals= 0;
there_is_a_bad_frame= FALSE;

Whanning= hann( Nf, 'periodic');
Whanning= Whanning';

D_POW_F = 2;
D_POW_S = 6;
D_POW_T = 2;
A_POW_F = 1;
A_POW_S = 6;
A_POW_T = 2;
D_WEIGHT= 0.1;
A_WEIGHT= 0.0309;

CRITERIUM_FOR_SILENCE_OF_5_SAMPLES = 500;
samples_to_skip_at_start = 0;
sum_of_5_samples= 0;
while ((sum_of_5_samples< CRITERIUM_FOR_SILENCE_OF_5_SAMPLES) ...
        && (samples_to_skip_at_start < maxNsamples / 2))
    sum_of_5_samples= sum( abs( ref_data( samples_to_skip_at_start...
        + SEARCHBUFFER * Downsample + 1: samples_to_skip_at_start...
        + SEARCHBUFFER * Downsample + 5)));

    if (sum_of_5_samples< CRITERIUM_FOR_SILENCE_OF_5_SAMPLES)
        samples_to_skip_at_start = samples_to_skip_at_start+ 1;
    end
end
% fprintf( 'samples_to_skip_at_start is %d\n', samples_to_skip_at_start);

samples_to_skip_at_end = 0;
sum_of_5_samples= 0;
while ((sum_of_5_samples< CRITERIUM_FOR_SILENCE_OF_5_SAMPLES) ...
        && (samples_to_skip_at_end < maxNsamples / 2))
    sum_of_5_samples= sum( abs( ref_data( maxNsamples - ...
        SEARCHBUFFER* Downsample + DATAPADDING_MSECS* (Fs/ 1000) ...
        - samples_to_skip_at_end - 4: maxNsamples - ...
        SEARCHBUFFER* Downsample + DATAPADDING_MSECS* (Fs/ 1000) ...
        - samples_to_skip_at_end)));
    if (sum_of_5_samples< CRITERIUM_FOR_SILENCE_OF_5_SAMPLES)
        samples_to_skip_at_end = samples_to_skip_at_end+ 1;
    end
end
% fprintf( 'samples_to_skip_at_end is %d\n', samples_to_skip_at_end);

start_frame = floor( samples_to_skip_at_start/ (Nf/ 2));
stop_frame = floor( (maxNsamples- 2* SEARCHBUFFER* Downsample ...
    + DATAPADDING_MSECS* (Fs/ 1000)- samples_to_skip_at_end) ...
    / (Nf/ 2))- 1;
% number of frames in speech data plus DATAPADDING_MSECS
% fprintf( 'start/end frame is %d/%d\n', start_frame, stop_frame);

D_disturbance= zeros( stop_frame+ 1, Nb);
DA_disturbance= zeros( stop_frame+ 1, Nb);

power_ref = pow_of (ref_data, SEARCHBUFFER* Downsample, ...
    maxNsamples- SEARCHBUFFER* Downsample+ DATAPADDING_MSECS* (Fs/ 1000),...
    maxNsamples- 2* SEARCHBUFFER* Downsample+ DATAPADDING_MSECS* (Fs/ 1000));
power_deg = pow_of (deg_data, SEARCHBUFFER * Downsample, ...
    maxNsamples- SEARCHBUFFER* Downsample+ DATAPADDING_MSECS* (Fs/ 1000),...
    maxNsamples- 2* SEARCHBUFFER* Downsample+ DATAPADDING_MSECS* (Fs/ 1000));
% fprintf( 'ref/deg power is %f/%f\n', power_ref, power_deg);

hz_spectrum_ref             = zeros( 1, Nf/ 2);
hz_spectrum_deg             = zeros( 1, Nf/ 2);
frame_is_bad                = zeros( 1, stop_frame + 1);
smeared_frame_is_bad        = zeros( 1, stop_frame + 1);
silent                      = zeros( 1, stop_frame + 1);

pitch_pow_dens_ref          = zeros( stop_frame + 1, Nb);
pitch_pow_dens_deg          = zeros( stop_frame + 1, Nb);

frame_was_skipped           = zeros( 1, stop_frame + 1);
frame_disturbance           = zeros( 1, stop_frame + 1);
frame_disturbance_asym_add  = zeros( 1, stop_frame + 1);

avg_pitch_pow_dens_ref      = zeros( 1, Nb);
avg_pitch_pow_dens_deg      = zeros( 1, Nb);
loudness_dens_ref           = zeros( 1, Nb);
loudness_dens_deg           = zeros( 1, Nb);
deadzone                    = zeros( 1, Nb);
disturbance_dens            = zeros( 1, Nb);
disturbance_dens_asym_add   = zeros( 1, Nb);

time_weight                 = zeros( 1, stop_frame + 1);
total_power_ref             = zeros( 1, stop_frame + 1);

% fid= fopen( 'tmp_mat.txt', 'wt');

for frame = 0: stop_frame
    start_sample_ref = 1+ SEARCHBUFFER * Downsample + frame* (Nf/ 2);
    hz_spectrum_ref= short_term_fft (Nf, ref_data, Whanning, ...
        start_sample_ref);

    utt = Nutterances;
    while ((utt >= 1) && ((Utt_Start(utt)- 1)* Downsample+ 1 ...
            > start_sample_ref))
        utt= utt - 1;
    end

    if (utt >= 1)
        delay = Utt_Delay(utt);
    else
        delay = Utt_Delay(1);
    end

    start_sample_deg = start_sample_ref + delay;

    if ((start_sample_deg > 0) && (start_sample_deg + Nf- 1 < ...
            maxNsamples+ DATAPADDING_MSECS* (Fs/ 1000)))
        hz_spectrum_deg= short_term_fft (Nf, deg_data, Whanning, ...
            start_sample_deg);
    else
        hz_spectrum_deg( 1: Nf/ 2)= 0;
    end

    pitch_pow_dens_ref( frame+ 1, :)= freq_warping (...
        hz_spectrum_ref, Nb, frame);
    %peak = maximum_of (pitch_pow_dens_ref, 0, Nb);
    pitch_pow_dens_deg( frame+ 1, :)= freq_warping (...
        hz_spectrum_deg, Nb, frame);

    total_audible_pow_ref = total_audible (frame, pitch_pow_dens_ref, 1E2);
    total_audible_pow_deg = total_audible (frame, pitch_pow_dens_deg, 1E2);
    silent(frame+ 1) = (total_audible_pow_ref < 1E7);
    

end
% fclose( fid);

avg_pitch_pow_dens_ref= time_avg_audible_of (stop_frame + 1, ...
    silent, pitch_pow_dens_ref, floor((maxNsamples- 2* SEARCHBUFFER* ...
    Downsample+ DATAPADDING_MSECS* (Fs/ 1000))/ (Nf / 2))- 1);
avg_pitch_pow_dens_deg= time_avg_audible_of (stop_frame + 1, ...
    silent, pitch_pow_dens_deg, floor((maxNsamples- 2* SEARCHBUFFER* ...
    Downsample+ DATAPADDING_MSECS* (Fs/ 1000))/ (Nf/ 2))- 1);

% fid= fopen( 'tmp_mat.txt', 'wt');
% fprintf( fid, '%f\n', avg_pitch_pow_dens_deg);
% fclose( fid);

if (CALIBRATE== 0)
    pitch_pow_dens_ref= freq_resp_compensation (stop_frame + 1, ...
        pitch_pow_dens_ref, avg_pitch_pow_dens_ref, ...
        avg_pitch_pow_dens_deg, 1000);
    if (Plot_Frame>= 0) % plot pitch_pow_dens_ref
        figure;
        subplot( 1, 2, 1);
        plot( centre_of_band_hz, 10* log10( eps+ ...
            pitch_pow_dens_ref( Plot_Frame+ 1, :)));
        axis( [0 Fs/2 0 95]); %xlabel( 'Hz'); ylabel( 'Db');   
        title( 'reference signal bark spectrum with frequency compensation');
        subplot( 1, 2, 2);
        plot( centre_of_band_hz, 10* log10( eps+ ...
            pitch_pow_dens_deg( Plot_Frame+ 1, :)));
        axis( [0 Fs/2 0 95]); %xlabel( 'Hz'); ylabel( 'Db');
        title( 'degraded signal bark spectrum');
    end
        
end
% tmp1= pitch_pow_dens_ref';


MAX_SCALE = 5.0;
MIN_SCALE = 3e-4;
oldScale = 1;
THRESHOLD_BAD_FRAMES = 30;
for frame = 0: stop_frame
    
    total_audible_pow_ref = total_audible (frame, pitch_pow_dens_ref, 1);
    total_audible_pow_deg = total_audible (frame, pitch_pow_dens_deg, 1);        
    total_power_ref (1+ frame) = total_audible_pow_ref;
    
    scale = (total_audible_pow_ref + 5e3)/ (total_audible_pow_deg + 5e3);    
    if (frame > 0) 
        scale = 0.2 * oldScale + 0.8 * scale;
    end
    oldScale = scale;
    
    if (scale > MAX_SCALE) 
        scale = MAX_SCALE;
    elseif (scale < MIN_SCALE) 
        scale = MIN_SCALE;            
    end

    pitch_pow_dens_deg( 1+ frame, :) = ...
        pitch_pow_dens_deg( 1+ frame, :) * scale;
    
    if (frame== Plot_Frame)
        figure;
        subplot( 1, 2, 1);
        plot( centre_of_band_hz, 10* log10( eps+ ...
            pitch_pow_dens_ref( Plot_Frame+ 1, :)));
        axis( [0 Fs/2 0 95]); %xlabel( 'Hz'); ylabel( 'Db');        
        subplot( 1, 2, 2);
        plot( centre_of_band_hz, 10* log10( eps+ ...
            pitch_pow_dens_deg( Plot_Frame+ 1, :)));
        axis( [0 Fs/2 0 95]); %xlabel( 'Hz'); ylabel( 'Db');
    end

    loudness_dens_ref = intensity_warping_of (frame, pitch_pow_dens_ref);
    loudness_dens_deg = intensity_warping_of (frame, pitch_pow_dens_deg);         
    disturbance_dens = loudness_dens_deg - loudness_dens_ref;
    
    if (frame== Plot_Frame)
        figure;
        subplot( 1, 2, 1);
        plot( centre_of_band_hz, 10* log10( eps+ ...
            loudness_dens_ref));
        axis( [0 Fs/2 0 15]); %xlabel( 'Hz'); ylabel( 'Db'); 
        title( 'reference signal loudness density');
        subplot( 1, 2, 2);
        plot( centre_of_band_hz, 10* log10( eps+ ...
            loudness_dens_deg));
        axis( [0 Fs/2 0 15]); %xlabel( 'Hz'); ylabel( 'Db');
        title( 'degraded signal loudness density');        
    end
    
    for band =1: Nb
        deadzone (band) = 0.25* min (loudness_dens_deg (band), ...
            loudness_dens_ref (band));    
    end

    for band = 1: Nb
        d = disturbance_dens (band);
        m = deadzone (band);
        
        if (d > m) 
            disturbance_dens (band) = disturbance_dens (band)- m;
%             disturbance_dens (band) = d- m;
        else
            if (d < -m) 
                disturbance_dens (band) = disturbance_dens (band)+ m;
%                 disturbance_dens (band) = d+ m;
            else
                disturbance_dens (band) = 0;
            end
        end
    end
    
    if (frame== Plot_Frame)
        figure;
        subplot( 1, 2, 1);
        plot( centre_of_band_hz, disturbance_dens);
        axis( [0 Fs/2 -1 50]); %xlabel( 'Hz'); ylabel( 'Db');                
        title( 'disturbance');        
    end
    D_disturbance( frame+ 1, :)= disturbance_dens;

    frame_disturbance (1+ frame) = pseudo_Lp (disturbance_dens, D_POW_F);    
    if (frame_disturbance (1+ frame) > THRESHOLD_BAD_FRAMES) 
        there_is_a_bad_frame = TRUE;
    end
    
    disturbance_dens= multiply_with_asymmetry_factor (...
        disturbance_dens, frame, pitch_pow_dens_ref, pitch_pow_dens_deg);
    
    if (frame== Plot_Frame)        
        subplot( 1, 2, 2);
        plot( centre_of_band_hz, disturbance_dens);
        axis( [0 Fs/2 -1 50]); %xlabel( 'Hz'); ylabel( 'Db');
        title( 'disturbance after asymmetry processing');
    end
    DA_disturbance( frame+ 1, :)= disturbance_dens;


    frame_disturbance_asym_add (1+ frame) = ...
        pseudo_Lp (disturbance_dens, A_POW_F);    
end
% fid= fopen( 'tmp_mat.txt', 'wt');
% fprintf( fid, '%f\n', frame_disturbance);
% fclose( fid);

frame_was_skipped (1: 1+ stop_frame) = FALSE;

for utt = 2: Nutterances
    frame1 = floor (((Utt_Start(utt)- 1- SEARCHBUFFER )* Downsample+ 1+ ...
        Utt_Delay(utt))/ (Nf/ 2));
    j = floor( floor(((Utt_End(utt-1)- 1- SEARCHBUFFER)* Downsample+ 1+ ...
        Utt_Delay(utt-1)))/(Nf/ 2));
    delay_jump = Utt_Delay(utt) - Utt_Delay(utt-1);
    if (frame1 > j) 
        frame1 = j;    
    elseif (frame1 < 0) 
        frame1 = 0;
    end
%     fprintf( 'frame1, j, delay_jump is %d, %d, %d\n', frame1, ...
%         j, delay_jump);

    if (delay_jump < -(Nf/ 2)) 
        frame2 = floor (((Utt_Start(utt)- 1- SEARCHBUFFER)* Downsample+ 1 ...
            + max (0, abs (delay_jump)))/ (Nf/ 2)) + 1; 
        
        for frame = frame1: frame2
            if (frame < stop_frame) 
                frame_was_skipped (1+ frame) = TRUE;
                frame_disturbance (1+ frame) = 0;
                frame_disturbance_asym_add (1+ frame) = 0;
            end
        end
    end
end

nn = DATAPADDING_MSECS* (Fs/ 1000) + maxNsamples;
tweaked_deg = zeros( 1, nn);
% fprintf( 'nn is %d\n', nn);

for i= SEARCHBUFFER* Downsample+ 1: nn- SEARCHBUFFER* Downsample
    utt = Nutterances;
    
    while ((utt >= 1) && ((Utt_Start (utt)- 1)* Downsample> i)) 
        utt = utt- 1;
    end
    if (utt >= 1) 
        delay = Utt_Delay (utt);        
    else
        delay = Utt_Delay (1);
    end

    j = i + delay;
    if (j < SEARCHBUFFER * Downsample+ 1) 
        j = SEARCHBUFFER * Downsample+ 1;
    end
    if (j > nn - SEARCHBUFFER * Downsample) 
        j = nn - SEARCHBUFFER * Downsample;
    end
    tweaked_deg (i) = deg_data (j);
end

if (there_is_a_bad_frame) 
    
    for frame = 0: stop_frame
        frame_is_bad (1+ frame) = (frame_disturbance (1+ frame)...
            > THRESHOLD_BAD_FRAMES);       
        smeared_frame_is_bad (1+ frame) = FALSE;
    end
    frame_is_bad (1) = FALSE;
    SMEAR_RANGE = 2;
    
    for frame = SMEAR_RANGE: stop_frame- 1- SMEAR_RANGE
        max_itself_and_left = frame_is_bad (1+ frame);
        max_itself_and_right = frame_is_bad (1+ frame);
        
        for i = -SMEAR_RANGE: 0
            if (max_itself_and_left < frame_is_bad (1+ frame+ i)) 
                max_itself_and_left = frame_is_bad (1+ frame+ i);
            end
        end

        for i = 0: SMEAR_RANGE
            if (max_itself_and_right < frame_is_bad (1+ frame + i)) 
                max_itself_and_right = frame_is_bad (1+ frame + i);
            end
        end

        mini = max_itself_and_left;
        if (mini > max_itself_and_right) 
            mini = max_itself_and_right;
        end

        smeared_frame_is_bad (1+ frame) = mini;
    end
    
    MINIMUM_NUMBER_OF_BAD_FRAMES_IN_BAD_INTERVAL = 5;
    number_of_bad_intervals = 0;    
    frame = 0; 
    while (frame <= stop_frame) 
        while ((frame <= stop_frame) && (~smeared_frame_is_bad (1+ frame)))
            frame= frame+ 1;
        end

        if (frame <= stop_frame) 
            start_frame_of_bad_interval(1+ number_of_bad_intervals)= ...
                1+ frame;
            
            while ((frame <= stop_frame) && (...
                    smeared_frame_is_bad (1+ frame))) 
                frame= frame+ 1; 
            end

            if (frame <= stop_frame)
                stop_frame_of_bad_interval(1+ number_of_bad_intervals)= ...
                    1+ frame; 
                if (stop_frame_of_bad_interval(1+ number_of_bad_intervals)- ...
                        start_frame_of_bad_interval(1+ number_of_bad_intervals)...
                        >= MINIMUM_NUMBER_OF_BAD_FRAMES_IN_BAD_INTERVAL) 
                    number_of_bad_intervals= number_of_bad_intervals+ 1;
                end
            end
        end
    end

    for bad_interval = 0: number_of_bad_intervals - 1
        start_sample_of_bad_interval(1+ bad_interval) = ...
            (start_frame_of_bad_interval(1+ bad_interval)- 1) * (Nf/ 2) ...
            + SEARCHBUFFER * Downsample+ 1;
        stop_sample_of_bad_interval(1+ bad_interval) = ...
            (stop_frame_of_bad_interval(1+ bad_interval)- 1) * (Nf/ 2) ...
            + Nf + SEARCHBUFFER* Downsample;
        if (stop_frame_of_bad_interval(1+ bad_interval) > stop_frame+ 1) 
            stop_frame_of_bad_interval(1+ bad_interval) = stop_frame+ 1; 
        end

        number_of_samples_in_bad_interval(1+ bad_interval) = ...
            stop_sample_of_bad_interval(1+ bad_interval) - ...
            start_sample_of_bad_interval(1+ bad_interval)+ 1;
    end        
%     fprintf( 'number of bad intervals %d\n', number_of_bad_intervals);
%     fprintf( '%d %d\n', number_of_samples_in_bad_interval(1), ...
%         number_of_samples_in_bad_interval(2));
%     fprintf( '%d %d\n', start_sample_of_bad_interval(1), ...
%         start_sample_of_bad_interval(2));

    SEARCH_RANGE_IN_TRANSFORM_LENGTH = 4;    
    search_range_in_samples= SEARCH_RANGE_IN_TRANSFORM_LENGTH * Nf;
    
    for bad_interval= 0: number_of_bad_intervals- 1
        ref = zeros (1, 2 * search_range_in_samples + ...
            number_of_samples_in_bad_interval (1+ bad_interval));
        deg = zeros (1, 2 * search_range_in_samples + ...
            number_of_samples_in_bad_interval (1+ bad_interval));
        
        ref(1: search_range_in_samples) = 0;

        ref (search_range_in_samples+ 1: search_range_in_samples+ ...
                number_of_samples_in_bad_interval (1+ bad_interval)) = ...
                ref_data (start_sample_of_bad_interval( 1+ bad_interval) + 1: ...
                start_sample_of_bad_interval( 1+ bad_interval) + ...
                number_of_samples_in_bad_interval (1+ bad_interval));
        
        ref (search_range_in_samples + ...
                number_of_samples_in_bad_interval (1+ bad_interval) + 1: ...
                search_range_in_samples + ...
                number_of_samples_in_bad_interval (1+ bad_interval) + ...
                search_range_in_samples) = 0;
        
        for i = 0: 2 * search_range_in_samples + ...
                number_of_samples_in_bad_interval (1+ bad_interval) - 1
            j = start_sample_of_bad_interval (1+ bad_interval) - ...
                search_range_in_samples + i;
            nn = maxNsamples - SEARCHBUFFER * Downsample + ...
                DATAPADDING_MSECS  * (Fs / 1000);
            if (j <= SEARCHBUFFER * Downsample) 
                j = SEARCHBUFFER * Downsample+ 1;
            end
            if (j > nn) 
                j = nn;
            end
            deg (1+ i) = tweaked_deg (j);
        end

        [delay_in_samples, best_correlation]= compute_delay ...
            (1, 2 * search_range_in_samples + ...
            number_of_samples_in_bad_interval (1+ bad_interval), ...
            search_range_in_samples, ref, deg);
        delay_in_samples_in_bad_interval (1+ bad_interval) =  ...
            delay_in_samples;
%         fprintf( 'delay_in_samples, best_correlation is \n\t%d, %f\n', ...
%             delay_in_samples, best_correlation);
%         
        if (best_correlation < 0.5) 
            delay_in_samples_in_bad_interval  (1+ bad_interval) = 0;
        end
    end

    if (number_of_bad_intervals > 0) 
        doubly_tweaked_deg = tweaked_deg( 1: maxNsamples + ...
            DATAPADDING_MSECS  * (Fs / 1000));
        for bad_interval= 0: number_of_bad_intervals- 1
            delay = delay_in_samples_in_bad_interval (1+ bad_interval);
        
            for i = start_sample_of_bad_interval (1+ bad_interval): ...
                    stop_sample_of_bad_interval (1+ bad_interval)
                j = i + delay;
                if (j < 1) 
                    j = 1;
                end
                if (j > maxNsamples) 
                    j = maxNsamples;
                end
                h = tweaked_deg (j);
                doubly_tweaked_deg (i) = h;
            end
        end

        untweaked_deg = deg_data;
        deg_data = doubly_tweaked_deg;
        
        for bad_interval= 0: number_of_bad_intervals- 1
            for frame = start_frame_of_bad_interval (1+ bad_interval): ...
                    stop_frame_of_bad_interval (1+ bad_interval)- 1
                frame= frame- 1;
                start_sample_ref = SEARCHBUFFER * Downsample + ...
                    frame * Nf / 2+ 1;
                start_sample_deg = start_sample_ref;
                hz_spectrum_deg= short_term_fft (Nf, deg_data, ...
                    Whanning, start_sample_deg);    
                pitch_pow_dens_deg( 1+ frame, :)= freq_warping (...
                    hz_spectrum_deg, Nb, frame);
            end

            oldScale = 1;
            for frame = start_frame_of_bad_interval (1+ bad_interval): ...
                    stop_frame_of_bad_interval (1+ bad_interval)- 1
                frame= frame- 1;    
                % see implementation for detail why 1 needed to be
                % subtracted
                total_audible_pow_ref = total_audible (frame, ...
                    pitch_pow_dens_ref, 1);
                total_audible_pow_deg = total_audible (frame, ...
                    pitch_pow_dens_deg, 1);        
                scale = (total_audible_pow_ref + 5e3) / ...
                    (total_audible_pow_deg + 5e3);
                if (frame > 0) 
                    scale = 0.2 * oldScale + 0.8*scale;
                end
                oldScale = scale;
                if (scale > MAX_SCALE) 
                    scale = MAX_SCALE;
                end
                if (scale < MIN_SCALE) 
                    scale = MIN_SCALE;   
                end

                pitch_pow_dens_deg (1+ frame, :) = ...
                    pitch_pow_dens_deg (1+ frame, :)* scale;
                loudness_dens_ref= intensity_warping_of (frame, ...
                    pitch_pow_dens_ref); 
                loudness_dens_deg= intensity_warping_of (frame, ...
                    pitch_pow_dens_deg); 
                disturbance_dens = loudness_dens_deg - loudness_dens_ref;
                
                for band = 1: Nb
                    deadzone(band) = min (loudness_dens_deg(band), ...
                        loudness_dens_ref(band));    
                    deadzone(band) = deadzone(band)* 0.25;
                end

                for band = 1: Nb
                    d = disturbance_dens (band);
                    m = deadzone (band);
                    
                    if (d > m) 
                        disturbance_dens (band) = ...
                            disturbance_dens (band)- m;
                    else
                        if (d < -m) 
                            disturbance_dens (band) = ...
                                disturbance_dens (band)+ m;
                        else
                            disturbance_dens (band) = 0;
                        end
                    end
                end

                frame_disturbance( 1+ frame) = min (...
                    frame_disturbance( 1+ frame), pseudo_Lp(...
                    disturbance_dens, D_POW_F));
                disturbance_dens= multiply_with_asymmetry_factor ...
                    (disturbance_dens, frame, pitch_pow_dens_ref, ...
                    pitch_pow_dens_deg);
                frame_disturbance_asym_add(1+ frame) = min (...
                    frame_disturbance_asym_add(1+ frame), ...
                    pseudo_Lp (disturbance_dens, A_POW_F));    
            end
        end
        deg_data = untweaked_deg;
    end
end     

for frame = 0: stop_frame
    h = 1;
    if (stop_frame + 1 > 1000) 
        n = floor( (maxNsamples - 2 * SEARCHBUFFER * Downsample)...
            / (Nf / 2)) - 1;
        timeWeightFactor = (n - 1000) / 5500;
        if (timeWeightFactor > 0.5) 
            timeWeightFactor = 0.5;
        end
        h = (1.0 - timeWeightFactor) + timeWeightFactor * frame / n;
    end

    time_weight (1 +frame) = h;
end

% fid= fopen( 'tmp_mat1.txt', 'at');
% fprintf( '\n');
for frame = 0: stop_frame
    h = ((total_power_ref (1+ frame) + 1e5) / 1e7)^ 0.04; 
%     if (frame== 118)
%         fprintf( '%f\n', h);    
%         fprintf( '%f\n', frame_disturbance( 1+ frame));
%     end
    frame_disturbance( 1+ frame) = frame_disturbance( 1+ frame)/ h;
    
%     if (frame== 118)
%         fprintf( '%f\n', frame_disturbance( 1+ frame));
%     end
%         
    frame_disturbance_asym_add( 1+ frame) = ...
        frame_disturbance_asym_add( 1+ frame)/ h;
    if (frame_disturbance( 1+ frame) > 45) 
        frame_disturbance( 1+ frame) = 45;  
    end
    if (frame_disturbance_asym_add( 1+ frame)> 45) 
        frame_disturbance_asym_add( 1+ frame) = 45;
    end
end
% fclose ( fid);

d_indicator = Lpq_weight (start_frame, stop_frame, ...
    D_POW_S, D_POW_T, frame_disturbance, time_weight);
a_indicator = Lpq_weight (start_frame, stop_frame, ...
    A_POW_S, A_POW_T, frame_disturbance_asym_add, time_weight);       

pesq_mos = 4.5 - D_WEIGHT * d_indicator - A_WEIGHT * a_indicator; 

if (Plot_Frame> 0)
    figure;
    subplot( 1, 2, 1);
    mesh( 0: stop_frame, centre_of_band_hz, D_disturbance');
    title( 'disturbance');
    subplot( 1, 2, 2);
    mesh( 0: stop_frame, centre_of_band_hz, DA_disturbance');
    title( 'disturbance after asymmetry processing');
end

% fid= fopen( 'tmp_mat.txt', 'wt');
% fprintf( fid, 'time_weight\n');
% fprintf( fid, '%f\n', time_weight);
% fprintf( fid, 'frame_disturbance:\n');
% fprintf( fid, '%f\n', frame_disturbance);
% fprintf( fid, 'frame_disturbance_asym_add\n');
% fprintf( fid, '%f\n', frame_disturbance_asym_add);
% fclose( fid);
    
function result_time= Lpq_weight(start_frame, stop_frame, ...
        power_syllable, power_time, frame_disturbance, time_weight)

global NUMBER_OF_PSQM_FRAMES_PER_SYLLABE

% fid= fopen( 'tmp_mat1.txt', 'at');
% fprintf( 'result_time:\n');

result_time= 0;
total_time_weight_time = 0;
% fprintf( 'start/end frame: %d/%d\n', start_frame, stop_frame);
for start_frame_of_syllable = start_frame: ...
        NUMBER_OF_PSQM_FRAMES_PER_SYLLABE/2: stop_frame
    result_syllable = 0;
    count_syllable = 0;
    
    for frame = start_frame_of_syllable: ...
            start_frame_of_syllable + NUMBER_OF_PSQM_FRAMES_PER_SYLLABE- 1
        if (frame <= stop_frame) 
            h = frame_disturbance(1+ frame);
%             if (start_frame_of_syllable== 101)
%                 fprintf( fid, '%f\n', h);
%             end
            result_syllable = result_syllable+ (h^ power_syllable);
        end
        count_syllable = count_syllable+ 1;
    end

    result_syllable = result_syllable/ count_syllable;
    result_syllable = result_syllable^ (1/power_syllable);     
    
    result_time= result_time+ (time_weight (...
        1+ start_frame_of_syllable - start_frame) * ...
        result_syllable)^ power_time; 
    total_time_weight_time = total_time_weight_time+ ...
        time_weight (1+ start_frame_of_syllable - start_frame)^ power_time;
    
%     fprintf( fid, '%f\n', result_time);
end
% fclose (fid);

% fprintf( 'total_time_weight_time is %f\n', total_time_weight_time);
result_time = result_time/ total_time_weight_time;
result_time= result_time^ (1/ power_time);
% fprintf( 'result_time is %f\n\n', result_time);

    
function [best_delay, max_correlation] = compute_delay (...
    start_sample, stop_sample, search_range, ...
    time_series1, time_series2) 

n = stop_sample - start_sample+ 1;   
power_of_2 = 2^ (ceil( log2( 2 * n)));

power1 = pow_of (time_series1, start_sample, stop_sample, n)* ...
    n/ power_of_2;
power2 = pow_of (time_series2, start_sample, stop_sample, n)* ...
    n/ power_of_2;
normalization = sqrt (power1 * power2);
% fprintf( 'normalization is %f\n', normalization);

if ((power1 <= 1e-6) || (power2 <= 1e-6)) 
    max_correlation = 0;
    best_delay= 0;
end

x1( 1: power_of_2)= 0;
x2( 1: power_of_2)= 0;
y( 1: power_of_2)= 0;

x1( 1: n)= abs( time_series1( start_sample: ...
    stop_sample));
x2( 1: n)= abs( time_series2( start_sample: ...
    stop_sample));

x1_fft= fft( x1, power_of_2)/ power_of_2;
x2_fft= fft( x2, power_of_2);
x1_fft_conj= conj( x1_fft);
y= ifft( x1_fft_conj.* x2_fft, power_of_2);

best_delay = 0;
max_correlation = 0;

% these loop can be rewritten
for i = -search_range: -1
    h = abs (y (1+ i + power_of_2)) / normalization;
    if (h > max_correlation) 
        max_correlation = h;
        best_delay= i;
    end
end
for i = 0: search_range- 1
    h = abs (y (1+i)) / normalization;
    if (h > max_correlation) 
        max_correlation = h;
        best_delay= i;
    end
end
best_delay= best_delay- 1;
    
function mod_disturbance_dens= multiply_with_asymmetry_factor (...
    disturbance_dens, frame, pitch_pow_dens_ref, pitch_pow_dens_deg) 

global Nb
for i = 1: Nb
    ratio = (pitch_pow_dens_deg(1+ frame, i) + 50)...
        / (pitch_pow_dens_ref (1+ frame, i) + 50);
    h = ratio^ 1.2;    
    if (h > 12) 
        h = 12;
    elseif (h < 3) 
        h = 0.0;
    end
    mod_disturbance_dens (i) = disturbance_dens (i) * h;
end


function loudness_dens = intensity_warping_of (...
    frame, pitch_pow_dens)

global abs_thresh_power Sl Nb centre_of_band_bark
ZWICKER_POWER= 0.23;
for band = 1: Nb
    threshold = abs_thresh_power (band);
    input = pitch_pow_dens (1+ frame, band);
    
    if (centre_of_band_bark (band) < 4) 
        h =  6 / (centre_of_band_bark (band) + 2);
    else
        h = 1;
    end

    if (h > 2) 
        h = 2;
    end
    h = h^ 0.15;
    modified_zwicker_power = ZWICKER_POWER * h;
    if (input > threshold) 
        loudness_dens (band) = ((threshold / 0.5)^ modified_zwicker_power)...
            * ((0.5 + 0.5 * input / threshold)^ modified_zwicker_power- 1);
    else
        loudness_dens (band) = 0;
    end

    loudness_dens (band) = loudness_dens (band)* Sl;
end
    
function result= pseudo_Lp (x, p)

global Nb width_of_band_bark
totalWeight = 0;
result = 0;
for band = 2: Nb
    h = abs (x (band));
    w = width_of_band_bark (band);
    prod = h * w;
    
    result = result+ prod^ p;
    totalWeight = totalWeight+ w;
end
result = (result/ totalWeight)^ (1/p);
result = result* totalWeight;

    
function mod_pitch_pow_dens_ref= freq_resp_compensation (number_of_frames, ...
    pitch_pow_dens_ref, avg_pitch_pow_dens_ref, ...
    avg_pitch_pow_dens_deg, constant)

global Nb

for band = 1: Nb
    x = (avg_pitch_pow_dens_deg (band) + constant) / ...
        (avg_pitch_pow_dens_ref (band) + constant);
    if (x > 100.0) 
        x = 100.0;
    elseif (x < 0.01) 
        x = 0.01;
    end

    for frame = 1: number_of_frames
        mod_pitch_pow_dens_ref(frame, band) = ...
            pitch_pow_dens_ref(frame, band) * x;
    end
end



function avg_pitch_pow_dens= time_avg_audible_of(number_of_frames, ...
    silent, pitch_pow_dens, total_number_of_frames) 

global Nb abs_thresh_power

for band = 1: Nb
    result = 0;
    for frame = 1: number_of_frames
        if (~silent (frame)) 
            h = pitch_pow_dens (frame, band);
            if (h > 100 * abs_thresh_power (band)) 
                result = result + h;
            end
        end

        avg_pitch_pow_dens (band) = result/ total_number_of_frames;
    end
end  



function hz_spectrum= short_term_fft (Nf, data, Whanning, start_sample)

x1= data( start_sample: start_sample+ Nf-1).* Whanning;
x1_fft= fft( x1);
hz_spectrum= abs( x1_fft( 1: Nf/ 2)).^ 2;
hz_spectrum( 1)= 0;


function pitch_pow_dens= freq_warping( hz_spectrum, Nb, frame)

global nr_of_hz_bands_per_bark_band pow_dens_correction_factor
global Sp

hz_band = 1;
for bark_band = 1: Nb
    n = nr_of_hz_bands_per_bark_band (bark_band);    
    sum = 0;
    for i = 1: n
        sum = sum+ hz_spectrum( hz_band);
        hz_band= hz_band+ 1;
    end
    sum = sum* pow_dens_correction_factor (bark_band);
    sum = sum* Sp;
    pitch_pow_dens (bark_band) = sum;
    
end


function total_audible_pow = total_audible (frame, ...
    pitch_pow_dens, factor)

global Nb abs_thresh_power

total_audible_pow = 0;
for band= 2: Nb
    h = pitch_pow_dens (frame+ 1,band);
    threshold = factor * abs_thresh_power (band);
    if (h > threshold) 
        total_audible_pow = total_audible_pow+ h;
    end
end








