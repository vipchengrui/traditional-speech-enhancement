clear;
clc;
addpath('bin')
%SSNR.exe CH01F051.dat vts.dat > old.txt
path = '.\EnhancedSpeech';
[enhancedFiles]=find_wav(path); 
numEnhancedFile=size(enhancedFiles,1);

fid=fopen('SSNR.txt','wt');

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
    
    system(['.\bin\SSNR ','temp_clean.dat',' ','temp_enhanced.dat',' > SSNR_result_singleFile.txt']);
    SSNR=importdata('SSNR_result_singleFile.txt');
    delete('SSNR_result_singleFile.txt')
    
    fprintf(fid,'%s\t\t',enhancedFilePath(length(path)+2:end));
    fprintf(fid,'%s\n',num2str(SSNR));
    
    fclose(fid_1);
    delete('temp_enhanced.dat');
end
fclose(fid_);
delete('temp_clean.dat');
fclose(fid);