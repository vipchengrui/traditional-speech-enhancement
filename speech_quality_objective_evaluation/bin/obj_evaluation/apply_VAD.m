function [VAD, logVAD]= apply_VAD( data, Nsamples)

global Downsample MINSPEECHLGTH JOINSPEECHLGTH

Nwindows= floor( Nsamples/ Downsample);
%number of 4ms window

VAD= zeros( 1, Nwindows);
for count= 1: Nwindows
    VAD( count)= sum( data( (count-1)* Downsample+ 1: ...
        count* Downsample).^ 2)/ Downsample;   
end
%VAD is the power of each 4ms window 

LevelThresh = sum( VAD)/ Nwindows;
%LevelThresh is set to mean value of VAD

LevelMin= max( VAD);
if( LevelMin > 0 )
    LevelMin= LevelMin* 1.0e-4;
else
    LevelMin = 1.0;
end
%fprintf( 1, 'LevelMin is %f\n', LevelMin);

VAD( find( VAD< LevelMin))= LevelMin;

for iteration= 1: 12    
    LevelNoise= 0;
    len= 0;
    StDNoise= 0;    
    
    VAD_lessthan_LevelThresh= VAD( find( VAD<= LevelThresh));
    len= length( VAD_lessthan_LevelThresh);
    LevelNoise= sum( VAD_lessthan_LevelThresh);
    if (len> 0)
        LevelNoise= LevelNoise/ len;
        StDNoise= sqrt( sum( ...
        (VAD_lessthan_LevelThresh- LevelNoise).^ 2)/ len);
    end
    LevelThresh= 1.001* (LevelNoise+ 2* StDNoise);  
end
%fprintf( 1, 'LevelThresh is %f\n', LevelThresh);

LevelNoise= 0;
LevelSig= 0;
len= 0;
VAD_greaterthan_LevelThresh= VAD( find( VAD> LevelThresh));
len= length( VAD_greaterthan_LevelThresh);
LevelSig= sum( VAD_greaterthan_LevelThresh);

VAD_lessorequal_LevelThresh= VAD( find( VAD<= LevelThresh));
LevelNoise= sum( VAD_lessorequal_LevelThresh);

if (len> 0)
    LevelSig= LevelSig/ len;
else
    LevelThresh= -1;
end
%fprintf( 1, 'LevelSig is %f\n', LevelSig);

if (len< Nwindows)
    LevelNoise= LevelNoise/( Nwindows- len);
else
    LevelNoise= 1;
end
%fprintf( 1, 'LevelNoise is %f\n', LevelNoise);

VAD( find( VAD<= LevelThresh))= -VAD( find( VAD<= LevelThresh));
VAD(1)= -LevelMin;
VAD(Nwindows)= -LevelMin;


start= 0;
finish= 0;
for count= 2: Nwindows
    if( (VAD(count) > 0.0) && (VAD(count-1) <= 0.0) )
        start = count;
    end
    if( (VAD(count) <= 0.0) && (VAD(count-1) > 0.0) )
        finish = count;
        if( (finish - start)<= MINSPEECHLGTH )
            VAD( start: finish- 1)= -VAD( start: finish- 1);
        end
    end
end
%to make sure finish- start is more than 4

if( LevelSig >= (LevelNoise* 1000) )
    for count= 2: Nwindows
        if( (VAD(count)> 0) && (VAD(count-1)<= 0) )
            start= count;
        end
        if( (VAD(count)<= 0) && (VAD(count-1)> 0) )
            finish = count;
            g = sum( VAD( start: finish- 1));
            if( g< 3.0* LevelThresh* (finish - start) )
                VAD( start: finish- 1)= -VAD( start: finish- 1);
            end
        end
    end
end

start = 0;
finish = 0;
for count= 2: Nwindows
    if( (VAD(count) > 0.0) && (VAD(count-1) <= 0.0) )
        start = count;
        if( (finish > 0) && ((start - finish) <= JOINSPEECHLGTH) )
            VAD( finish: start- 1)= LevelMin;
        end        
    end
    if( (VAD(count) <= 0.0) && (VAD(count-1) > 0.0) )
        finish = count;
    end
end

start= 0;
for count= 2: Nwindows
    if( (VAD(count)> 0) && (VAD(count-1)<= 0) )
        start= count;
    end
end
if( start== 0 )
    VAD= abs(VAD);
    VAD(1) = -LevelMin;
    VAD(Nwindows) = -LevelMin;
end

count = 4;
while( count< (Nwindows-1) )
    if( (VAD(count)> 0) && (VAD(count-2) <= 0) )
        VAD(count-2)= VAD(count)* 0.1;
        VAD(count-1)= VAD(count)* 0.3;
        count= count+ 1;
    end
    if( (VAD(count)<= 0) && (VAD(count-1)> 0) )
        VAD(count)= VAD(count-1)* 0.3;
        VAD(count+ 1)= VAD(count-1)* 0.1;
        count= count+ 3;
    end
    count= count+ 1;
end

VAD( find( VAD< 0))= 0;

% fid= fopen( 'mat_vad.txt', 'wt');
% fprintf( fid, '%f\n', VAD);
% fclose( fid);

if( LevelThresh<= 0 )
    LevelThresh= LevelMin;
end

logVAD( find( VAD<= LevelThresh))= 0;
VAD_greaterthan_LevelThresh= find( VAD> LevelThresh);
logVAD( VAD_greaterthan_LevelThresh)= log( VAD( ...
    VAD_greaterthan_LevelThresh)/ LevelThresh);




