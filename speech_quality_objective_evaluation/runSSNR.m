clear;
clc;
addpath('bin')
%SSNR.exe CH01F051.dat vts.dat > old.txt  
[enhancedFiles]=find_wav('.\EnhancedSpeech\1' ); 
[noisyFiles] = find_wav('.\Noisy_irm');
numEnhancedFile=size(enhancedFiles,1);
numNoisyFile=size(noisyFiles,1);
assert(numEnhancedFile == numNoisyFile, ['error!']);

fid=fopen('SSNR_result1.txt','wt');

cleanFilePath='.\cleanSpeech\TIMIT_1_TEST.wav';
cleanSig = audioread(cleanFilePath);
cleanSig = cleanSig.*32768./max(cleanSig);
fid_ = fopen('temp_clean.dat','w');
fwrite(fid_,cleanSig,'short');

for i=1:numEnhancedFile

    enhancedFilePath=enhancedFiles(i,:);
    noisyFilePath = noisyFiles(i,:);
    
    enhancedSig = audioread(enhancedFilePath);
    noisySig = audioread(noisyFilePath);
    enhancedSig = enhancedSig.*32768./max(enhancedSig);
    noisySig = noisySig.*32768./max(noisySig);
    fid_1 = fopen('temp_enhanced.dat','w');
    fwrite(fid_1,enhancedSig,'short');
    fid_2 = fopen('temp_noisy.dat','w');
    fwrite(fid_2,noisySig,'short');
    
    system(['.\bin\SSNR ','temp_clean.dat',' ','temp_noisy.dat',' > SSNR_result_singleFile.txt']);
    SSNR1=importdata('SSNR_result_singleFile.txt');
    delete('SSNR_result_singleFile.txt')
    
    system(['.\bin\SSNR ','temp_clean.dat',' ','temp_enhanced.dat',' > SSNR_result_singleFile.txt']);
    SSNR2=importdata('SSNR_result_singleFile.txt');
    delete('SSNR_result_singleFile.txt')
    
    SSNRI=SSNR2-SSNR1;
    fprintf(fid,'%s\t\t',enhancedFilePath);
    fprintf(fid,'%s\n',num2str(SSNRI));
    
    fclose(fid_1);
    fclose(fid_2);
    delete('temp_enhanced.dat');
    delete('temp_noisy.dat');
end
fclose(fid_);
delete('temp_clean.dat');
fclose(fid);