function mod_data= apply_filters( data, Nsamples)
%IIRFilt( InIIR_Hsos, InIIR_Nsos, data, data_Nsamples);

global InIIR_Hsos InIIR_Nsos DATAPADDING_MSECS Fs
% data_Nsamples= Nsamples+ DATAPADDING_MSECS* (Fs/ 1000);

% now we construct the second order section matrix
sosMatrix= zeros( InIIR_Nsos, 6);
sosMatrix( :, 4)= 1; %set a(1) to 1
% each row of sosMatrix holds [b(1*3) a(1*3)] for each section
sosMatrix( :, 1: 3)= InIIR_Hsos( :, 1: 3);
sosMatrix( :, 5: 6)= InIIR_Hsos( :, 4: 5);
%sosMatrix

% now we construct second order section direct form II filter
iirdf2= dfilt.df2sos( sosMatrix);

mod_data= filter( iirdf2, data);








