function [SNRI,SNRO,segSNRI,segSNRO,LSDI,LSDO,PESQ_I,PESQ_O]=res_eva(spe,noi,eha,fs,frame,fstep)
% res_eva result evaluation 结果评价，包括信噪比，分段信噪比，LSD（log-spectral
% distortion），噪声估计错误segERR
%
%输入：spe 纯净语音信号 noi 纯净噪声信号 eha 增强后语音信号 （noisy）
%      frame 窗口长度，fstep 移动步长
%输出：   SNRI，SNRO        输入输出信噪比
%      segSNRI，segSNRI     输入输出分段信噪比
%         LSDI，LSDO        输入输出LSD
%         LSD2I，LSD2O      算法2求出的LSD， Frames are 25ms long with 60 percent
%                           (15ms) overlap，第四个参数表示是否去掉静音段，0不去，1去
%         PESQ_I PESQ_O      输入与输出PESQ
%版本与版权：
%        res_eva v1.0  2010-11-07   原型
%        res_eva v2.0  2011-01-07  （1）添加了PESQ；（2）添加了LSD2，作为参考
%               Copyright (c) 2010. Infocarrier.
%               All rights reserved. 
%% 将三段语音对齐
le=length(eha);
spe=spe(1:le);
noi=noi(1:le);
noisy=noi+spe;
%% 信噪比,SNRI 输入信噪比，SNRO 输出信噪比
SNR1=mean(spe.^2)/mean(noi.^2);
SNRI=10*log10(SNR1);

SNR1=mean(spe.^2)/mean((spe-eha).^2);
SNRO=10*log10(SNR1);

%% 分段信噪比 segSNR，segSNRI 输入分段信噪比，SNRO 输出分段信噪比
N=frame;
M=fstep;
le=length(eha);
ss=fix(le/fstep); 
ts=ss-3;
SNRt=zeros(1,ss);
for t=1:ts
    SNRt(t)=10*log10(sum(spe((t-1)*M+1:(t-1)*M+N).^2)/sum(noi((t-1)*M+1:(t-1)*M+N).^2));
end
SNRtp=min(max(SNRt,-10),35);
segSNRI=sum(SNRtp)/ts;

for t=1:ts
    SNRt(t)=10*log10(sum(spe((t-1)*M+1:(t-1)*M+N).^2)/...
sum((spe((t-1)*M+1:(t-1)*M+N)-eha((t-1)*M+1:(t-1)*M+N)).^2));% 公式44.95,96
end
SNRtp=min(max(SNRt,-10),35);
segSNRO=sum(SNRtp)/ts;

%% LSD

%%%纯净语音
speX=spectrogram(spe,frame,frame-fstep,frame);
speX=20*log10(abs(speX));%speX(K,T),列是时间列
theta=max(max(speX))-50;
LspeX=max(speX,theta);
%%%带噪语音
noisY=spectrogram(noisy,frame,frame-fstep,frame);
noisY=20*log10(abs(noisY));
thetaY=max(max(noisY))-50;
LnoisY=max(noisY,thetaY);
%%%增强后语音
speX1=spectrogram(eha,frame,frame-fstep,frame);
speX1=20*log10(abs(speX1));
theta1=max(max(speX1))-50;
LspeX1=max(speX1,theta1);
%%% 输入LSD
LnoisY(1,:)=0;
LspeX(1,:)=0;
LspeX1(1,:)=0;%傅里叶系数的第一个DC数置零，下式没有用到
LSDI=(1/ts)*sum((2/N*sum((LnoisY-LspeX).^2)).^(1/2));
%%%计算输出LSD
LSDO=(1/ts)*sum((2/N*sum((LspeX-LspeX1).^2)).^(1/2));
%% LSD2
% LSD2I=LogSpectralDistance(spe,noisy,fs,0);%第四个参数，0不去掉静音段，1去掉静音段
% LSD2O=LogSpectralDistance(spe,eha,fs,0);%第四个参数，0不去掉静音段，1去掉静音段
%% PESQ
%%% [pesq_mos]=pesq(sfreq,cleanfile.wav,enhanced.wav) 
addpath('eva_composite');
fname_noisy='temp_noisy.wav';
fname_spe='temp_spe.wav';
fname_eha='temp_eha.wav';
wavwrite(noisy,fs,fname_noisy);
wavwrite(spe,fs,fname_spe);
wavwrite(eha,fs,fname_eha);
PESQ_I=pesq(fs,fname_spe,fname_noisy); 
PESQ_O=pesq(fs,fname_spe,fname_eha); 
delete(fname_noisy);
delete(fname_spe);
delete(fname_eha);


