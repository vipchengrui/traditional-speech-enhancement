clear;
clc;
%PESQ_LQ.exe +8000 CH01M001.dat CH01M001-white+10dB.dat > test.txt
addpath('bin');
addpath('bin\obj_evaluation')
path = '.\EnhancedSpeech\';
[enhancedFiles]=find_wav(path); 
numEnhancedFile=size(enhancedFiles,1);
fid=fopen('PESQ.txt','wt');

cleanFilePath='.\cleanSpeech\TIMIT_1_TEST.wav';
cleanSig = audioread(cleanFilePath);
for i=1:numEnhancedFile
    enhancedFilePath=enhancedFiles(i,:);
    
    k=0;
    enhancedSig = audioread(enhancedFilePath);
    short = min(length(enhancedSig),length(cleanSig));
    
    cleanSig = cleanSig(1:1:short);
    enhancedSig = enhancedSig(1:1:short);
   
    audiowrite('temp_clean.wav',cleanSig,8000,'BitsPerSample',16);
    audiowrite('temp_enhanced.wav',enhancedSig,8000,'BitsPerSample',16);
    PESQ_O=pesq('temp_clean.wav','temp_enhanced.wav'); 
    delete('temp_clean.wav');
    delete('temp_enhanced.wav');
    
    fprintf(fid,'%s\t\t',enhancedFilePath(length(path)+2:end));
    fprintf(fid,'%s\n',num2str(PESQ_O));
end
fclose(fid);