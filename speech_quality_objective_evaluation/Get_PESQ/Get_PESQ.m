clc
clear all

speechfile = 'female';
SNR = '12dB';
eva = zeros(10,5);
% %-------------------------------clean-----------------
% savedir = ['C:\Users\74114\Desktop\Multi-Channel_Dataset\时延估计验证\re_clean\single\'];
% speechdir = ['C:\Users\74114\Desktop\Multi-Channel_Dataset\时延估计验证\clean\'];
% for ch = 1: 10
% % 
%  sp=[speechdir ['micT_clean',num2str(ch),'.wav']];
% resp = [savedir ['re-micT_clean','_',speechfile,'-',num2str(ch),'_codec','.wav']];
% % sp='C:\Users\74114\Desktop\Multi-Channel_Dataset\时延估计验证\6dB\mic\micT_babble_female-3-6dB.wav';
% % resp='C:\Users\74114\Desktop\Multi-Channel_Dataset\时延估计验证\re_6dB\single\re-micT_babble_female-3-6dB_codec-test.wav';
%    PESQ_01 = PESQ(sp, resp);
%    eva(ch,1) = PESQ_01;
% end


%%--------------------------------------------------noise-----------
savedir = ['D:\Python_Project\Multi_Codec\Experiment\EVS_ch1_ch10\re_',SNR,'\'];
speechdir = ['D:\Python_Project\Multi_Codec\Experiment\EVS_ch1\',SNR,'\mic\'];
for  Pnum =1:5
 switch Pnum
                case 1
                    noisetype = 'babble';
                case 2
                    noisetype = 'Office';
                case 3
                    noisetype = 'street';
                case 4
                    noisetype = 'VOLVO';
                case 5
                    noisetype = 'White';
 end
 
% [ch1, fs] = audioread([bitdir ['micT_',noisetype,'_',speechfile,'-1-',SNR,'_codec','.wav']]);
%  


for ch = 1: 10

 sp=[speechdir ['micT_',noisetype,'_',speechfile,'-',num2str(ch),'-',SNR,'.wav']];
 resp = [savedir ['re-micT_',noisetype,'_',speechfile,'-',num2str(ch),'-',SNR,'_codec','.wav']];
  PESQ_01 = PESQ(sp, resp);
  eva(ch,Pnum) = PESQ_01;
%  fprintf(fid,'%s\t\t',['re-micT_',noisetype,'_',speechfile,'-',num2str(ch),'-',SNR,'_codec']);

end

 
 
end






