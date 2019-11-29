%omlsa+wavelete 自动化处理+评估
%% 加噪声，调整输入信噪比，得到带噪语音
noi2=noi*1.1;                   %后面的整数代表噪声的含量
noisy=clean+noi2;               % 带噪语音生成
[SNRI,SNRO,segSNRI,segSNRO,LSDI,LSDO,PESQ_I,PESQ_O]=res_eva(clean,noi2,noisy,fs,frame,fstep);
clear SNRO segSNRO LSDO IPESQ_O;
fout=['d:/noisy_' int2str(SNRI) 'db'];
wavwrite(noisy,fs,fout);
segSNRI
%% omlsa增强处理
fin_om=fout;
fout_om=omlsa_ncre(fin_om);
%% 小波增强处理
fin_wl=fout_om;
fout_wl=wl_thre(fin_wl);
%% 评估增强效果
noisy_omlsa=wavread(fout_om);
noisy_omlsa_wavelet=wavread(fout_wl);
[omsla_iSNR,omsla_isegSNR,omsla_iLSD,omsla_iPESQ]=res_imp(clean,noi2,noisy_omlsa,fs,frame,fstep);
[wavelet_iSNR,wavelet_isegSNR,wavelet_iLSD,wavelet_iPESQ]=res_imp(clean,noi2,noisy_omlsa_wavelet,fs,frame,fstep);

%% 收尾结束
del_file=1;% 是否删除输出的文件
if del_file
    delete([fout '.wav']);
    delete([fout_om '.wav']);
    delete([fout_wl '.wav']);
end
msgbox('程序运行完毕','结束标志','warn');
