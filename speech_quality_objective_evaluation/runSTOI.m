clear;
clc;
%PESQ_LQ.exe +8000 CH01M001.dat CH01M001-white+10dB.dat > test.txt
addpath('bin');
[enhancedFiles]=find_wav('.\EnhancedSpeech' ); 
numEnhancedFile=size(enhancedFiles,1);
fid=fopen('STOI.txt','wt');

cleanFilePath='.\cleanSpeech\TIMIT_1_TEST.wav';
[cleanSig,fs_signal] = audioread(cleanFilePath);
for i=1:numEnhancedFile
    enhancedFilePath=enhancedFiles(i,:);
        
    enhancedSig = audioread(enhancedFilePath);
    short = min(length(enhancedSig),length(cleanSig));
    enhancedSig = enhancedSig(1:1:short);
    cleanSig = cleanSig(1:1:short);
   
    STOI = stoi(cleanSig, enhancedSig, fs_signal);   %¿É¸ÄžéœySTOI
    
    fprintf(fid,'%s\t\t',enhancedFilePath);
    fprintf(fid,'%s\n',num2str(STOI));
end
fclose(fid);