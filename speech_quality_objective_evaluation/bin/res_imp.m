function [iSNR,isegSNR,iLSD,iPESQ]=res_imp(spe,noi,eha,fs,frame,fstep)
% 计算各指标的提高量
[SNRI,SNRO,segSNRI,segSNRO,LSDI,LSDO,PESQ_I,PESQ_O]=res_eva(spe,noi,eha,fs,frame,fstep);
iSNR=SNRO-SNRI;
isegSNR=segSNRO-segSNRI;
iLSD=LSDI-LSDO;
iPESQ=PESQ_O-PESQ_I;