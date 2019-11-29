clear;
clc;
addpath('./bin');
%LogSD_FFT.exe CH01M001.dat CH01M001-white+10dB.dat > test.txt 
[enhancedFiles]=find_wav('.\EnhancedSpeech' ); 
numEnhancedFile=size(enhancedFiles,1);
fid=fopen('LSD.txt','wt');

cleanFilePath='.\cleanSpeech\TIMIT_1_TEST.wav';
cleanSig = audioread(cleanFilePath);
cleanSig = cleanSig.*32768;
fid_ = fopen('temp_clean.dat','w');
fwrite(fid_,cleanSig,'short');

for i=1:numEnhancedFile
    enhancedFilePath=enhancedFiles(i,:);
    enhancedSig = audioread(enhancedFilePath);
    enhancedSig = enhancedSig.*32768;
    fid_1 = fopen('temp_enhanced.dat','w');
    fwrite(fid_1,enhancedSig,'short');
    
    
    system(['.\bin\LogSD_AR ','temp_clean.dat',' ','temp_enhanced.dat',' > LSD_result_singleFile.txt']);
    LSD=importdata('LSD_result_singleFile.txt');
    delete('LSD_result_singleFile.txt');
    fclose(fid_1);
    delete('temp_enhanced.dat');
    fprintf(fid,'%s\t\t',enhancedFilePath);
    fprintf(fid,'%s\n',num2str(LSD));
end
fclose(fid_);
delete('temp_clean.dat');
fclose(fid);