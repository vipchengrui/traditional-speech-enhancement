function [ wav_files ] = find_dat( path )  
%FIND_WAV, find all wav file recursively  
wav_files = [];  
if(isdir(path) == 0)  
    
    return;  
end  
path_files = dir(path);  
fileNum = length(path_files);  
for k= 3:fileNum  
    file = [path,'\', path_files(k).name];  
    if (path_files(k).isdir == 1)  
        ret = find_wav(file);  
        if(isempty(ret) ~= 1)  
            if(isempty(wav_files))  
                wav_files = char(ret);  
            else  
                wav_files = char(wav_files, ret);  
            end  
        end  
    elseif strfind(path_files(k).name, '.dat')  
        if(isempty(wav_files))  
            wav_files = char(file);  
        else  
            wav_files = char(wav_files, file);  
        end  
    end  
end  
end  

