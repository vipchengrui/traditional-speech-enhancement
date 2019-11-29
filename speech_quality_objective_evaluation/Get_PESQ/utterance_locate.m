function utterance_locate (ref_data, ref_Nsamples, ref_VAD, ref_logVAD,...
    deg_data, deg_Nsamples, deg_VAD, deg_logVAD);

global Nutterances Utt_Delay Utt_DelayConf Utt_Start Utt_End Utt_DelayEst

id_searchwindows( ref_VAD, ref_Nsamples, deg_VAD, deg_Nsamples);

for Utt_id= 1: Nutterances
    %fprintf( 1, 'Utt_id is %d\n', Utt_id);
    crude_align( ref_logVAD, ref_Nsamples, deg_logVAD, deg_Nsamples, Utt_id);
    time_align(ref_data, ref_Nsamples, ...
        deg_data, deg_Nsamples, Utt_id);
end

id_utterances( ref_Nsamples, ref_VAD, deg_Nsamples);


utterance_split( ref_data, ref_Nsamples, ref_VAD, ref_logVAD, ...
    deg_data, deg_Nsamples, deg_VAD, deg_logVAD); 







