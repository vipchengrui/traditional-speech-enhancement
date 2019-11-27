import numpy as np
import wave
import matplotlib.pyplot as plt
import math
#import nextpow2

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

# noisy speech FFT
#x1_FFT = abs(np.fft.fft(x1))

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
print('11111111111')
print(x)

# noisy speech FFT
x_FFT = abs(np.fft.fft(x))

# calculation parameters
len_ = 20 * fs // 1000      # frame size in samples
PERC = 50                   # window overlop in percent of frame
len1 = len_ * PERC // 100   # overlop'length
len2 = len_ - len1          # window'length - overlop'length

# setting default parameters
Thres = 3       # VAD threshold in dB SNRseg
Expnt = 2.0
beta = 0.002    
G = 0.9

# hamming window
#win = np.hamming(len_)

# sine window
i = np.linspace(0,len_ - 1,len_)
win = np.sqrt(2/(len_ + 1)) * np.sin(np.pi * (i + 1) / (len_ + 1))

# normalization gain for overlap+add with 50% overlap
winGain = len2 / sum(win)

# nFFT = 2 * 2 ** (nextpow2.nextpow2(len_))
nFFT = 2 * 2 ** 8
noise_mean = np.zeros(nFFT)
j = 1
for k in range(1, 6):
    noise_mean = noise_mean + abs(np.fft.fft(win * x[j : j + len_] , nFFT))
    j = j + len_ 
noise_mu = abs(noise_mean / 5)

# initialize various variables
k = 1
img = 1j
x_old = np.zeros(len1)
Nframes = len(x) // len2 - 1
xfinal = np.zeros(Nframes * len2)

# === Start Processing ==== #
for n in range(0, Nframes):

    # Windowing
    insign = win * x[k - 1 : k + len_ - 1]    
    # compute fourier transform of a frame
    spec = np.fft.fft(insign, nFFT)    
    # compute the magnitude
    sig = abs(spec)     
    # save the noisy phase information
    theta = np.angle(spec)  
    
    # 1 over subtraction
    sub_speech = sig ** Expnt - noise_mu ** Expnt;
    # the pure signal is less than the noise signal power
    SNRpos = 10 * np.log10(np.linalg.norm(sig, 2) ** 2 / np.linalg.norm(noise_mu, 2) ** 2)

    
    # Priori SNR
    SNRpri = 10 * np.log10(np.linalg.norm(sub_speech ** (1 / Expnt), 2) ** 2 / np.linalg.norm(noise_mu, 2) ** 2)
    # parameter to deal mel
    mel_max = 10
    mel_0 = (1 + 4 * mel_max) / 5
    s = 25 / (mel_max - 1)
    # deal mel

    mel = 0
    if -5.0 <= SNRpri <= 20.0:
        mel = mel_0 - SNRpri / s
    else:
        if SNRpri < -5.0:
            mel = mel_max
        if SNRpri > 20:
            mel = 1

    

    # 2 gain function Gk
    G_k = sub_speech / (sub_speech + mel * noise_mu ** Expnt)
    wf_speech = G_k * sig
    
    # --- implement a simple VAD detector --- #
    if SNRpos < Thres:  # Update noise spectrum
        noise_temp = G * noise_mu ** Expnt + (1 - G) * sig ** Expnt  # Smoothing processing noise power spectrum
        noise_mu = noise_temp ** (1 / Expnt)  # New noise amplitude spectrum
    
    # add phase    
    #wf_speech[nFFT // 2 + 1:nFFT] = np.flipud(wf_speech[1:nFFT // 2])
    x_phase = wf_speech * np.exp(img * theta)

    # take the IFFT
    xi = np.fft.ifft(x_phase).real
    
    # --- Overlap and add --- #
    xfinal[k - 1 : k + len2 - 1] = x_old + xi[0 : len1]
    x_old = xi[0 + len1 : len_]

    k = k + len2

# save wave file
wf = wave.open('out_SNR15_sp01000.wav', 'wb')

# setting parameters
wf.setparams(params)
# set waveform file .tostring()Convert array to data
wave_data = (winGain * xfinal).astype(np.short)
wf.writeframes(wave_data.tostring())
# close wave file
wf.close()

# enchanced speech FFT
es_FFT = abs(np.fft.fft(winGain * xfinal))

# plot wave
plt.figure(1)
plt.suptitle('Wiener Filtering(OS) (SNR=15dB)')
plt.subplot(311)
plt.plot(x1)
plt.title('Clean Speech')
plt.subplot(312)
plt.title('Noisy Speech')
plt.plot(x)
plt.subplot(313)
plt.title('Enhanced Speech')
plt.plot(winGain * xfinal)
plt.show()
