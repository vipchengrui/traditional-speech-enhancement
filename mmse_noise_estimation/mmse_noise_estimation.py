import numpy as np
import wave
import matplotlib.pyplot as plt
import math 
import scipy.integrate as inte
import Init_Noises as IN
import Est_Noises as EN

# input wave file 
f1 = wave.open('sp01.wav')

# read format information
# (nchannels, sampwidth, framerate, nframes, comptype, compname)
params1 = f1.getparams()
nchannels1, sampwidth1, framerate1, nframes1 = params1[:4]
fs1 = framerate1
# read wave data
str_data1 = f1.readframes(nframes1)
# close .wav file
f1.close()

# convert waveform data to an array
x1 = np.fromstring(str_data1, dtype=np.short)

# input wave file 
f = wave.open('in_SNR15_sp01.wav')

# read format information
# (nchannels, sampwidth, framerate, nframes, comptype, compname)
params = f.getparams()
nchannels, sampwidth, framerate, nframes = params[:4]
fs = framerate
# read wave data
str_data = f.readframes(nframes)
# close .wav file
f.close()

# convert waveform data to an array
x = np.fromstring(str_data, dtype=np.short)

# calculation parameters
len_ = 20 * fs // 1000      # frame size in samples
PERC = 50                   # window overlop in percent of frame
len1 = len_ * PERC // 100    # overlop'length
len2 = len_ - len1          # window'length - overlop'length

# setting default parameters
aa = 0.98
eta = 0.15
Thres = 3
mu = 0.98
c = np.sqrt(np.pi) / 2
qk = 0.3
ksi_min = 10 ** (-25 / 10)   #-25dB

# hamming window
win = np.hamming(len_)
# normalization gain for overlap+add with 50% overlap
winGain = len2 / sum(win)

# setting inital noise
nFFT = 2 * 2 ** 8

# initialize various variables
k = 1
img = 1j
x_old = np.zeros(len2)
Nframes = len(x) // len2 - 1
xfinal = np.zeros(Nframes * len2)

# === Start Processing ==== #
for n in range(0, Nframes):

    # Windowing
    insign = win * x[k - 1 : k + len_ - 1]

    # Take fourier transform of frame
    spec = np.fft.fft(insign , nFFT)
    sig = abs(spec)
    sig2 = sig ** 2
    # save the noisy phase information
    theta = np.angle(spec)  

    # Noise Estimation
    #Init_Weight、ConMinTrack、MCRA、MCRA2
    if n == 0:
        para = IN.Init_MCRA2(sig2,fs).info()    
    else:
        para = EN.Est_MCRA2(sig2,para).est()

    noise_ps = para['noise_ps']
    noise_mu = np.sqrt(noise_ps)
    noise_mu2 = noise_mu ** 2

    SNRpos = 10 * np.log10(np.linalg.norm(sig, 2) ** 2 / np.linalg.norm(noise_mu, 2) ** 2)

    # posteriori SNR
    gammak = np.minimum(sig2 / noise_mu2 , 40) 
    
    # decision-direct estimate of a priori SNR   P231 [7.75]
      
    if n == 0:
        ksi = aa + (1 - aa) * np.maximum(gammak - 1 , 0)
    else:
        ksi = aa * Xk_prev / noise_mu2 + (1 - aa) * np.maximum(gammak - 1 , 0)
        # limit ksi to -25 dB 
        ksi = np.maximum(ksi_min , ksi)  

    #Log-MMSE estimator[7.84]
    def integrand(t):
        return np.exp(-t) / t
    A = ksi / (1 + ksi) 
    vk = A * gammak
    ei_vk = np.zeros(len(vk))
    for i in range(len(vk)): 
        ei_vk[i] = 0.5 * inte.quad(integrand,vk[i],np.inf)[0]
    Glsa = A * np.exp(ei_vk)
 
    # get X(w)  [7.205]
    pSAP = (1 - qk) / (1 - qk + qk * (1 + ksi) * np.exp(-vk))
    Gmin = 10 ** (-20 / 10)   #-20dB
    hw = (Glsa ** pSAP) * (Gmin ** (1 - pSAP))
    mmse_speech = hw * sig

    # save for estimation of a priori SNR in next frame
    Xk_prev = mmse_speech ** 2  

    # IFFT
    x_phase = mmse_speech * np.exp(img * theta)
    xi_w = np.fft.ifft(x_phase , nFFT).real

    # overlop add
    xfinal[k - 1 : k + len2 - 1] = x_old + xi_w[0 : len1]
    x_old = xi_w[len1 + 0 : len_]

    k = k + len2
    
# save wave file
wf = wave.open('out_SNR15_sp01.wav', 'wb')

# setting parameters
wf.setparams(params)
# set waveform file .tostring()Convert array to data
wave_data = (winGain * xfinal).astype(np.short)
wf.writeframes(wave_data.tostring())
# close wave file
wf.close()

# plot wave
plt.figure(1)
plt.suptitle('Log-MMSE_SPU based on Noise Estimation(MCRA2) (SNR=15dB)')
plt.subplot(311)
plt.plot(x1)
plt.title('Clean Speech')
plt.subplot(312)
plt.plot(x)
plt.title('Noisy Speech')
plt.subplot(313)
plt.plot(winGain * xfinal)
plt.title('Enhanced Speech')
plt.show()
