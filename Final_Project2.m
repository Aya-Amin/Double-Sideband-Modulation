%Saves the actual samples of the file in modulatingSignal & its sampling
%rate in Fs
%modulatingSignal is an array of samples. Fs, a number,equals 48000
[modulatingSignal, Fs] = audioread('eric.wav');

%To visualize information in modulatingSignal ...
plot(modulatingSignal);
title('Input modulating signal');
xlabel('sample number');
ylabel('amplitude');
figure;

%Get modulatingSignal's waveform (plot of amp vs time)...
%Fs is in samples/sec. Sampling period, Ts, in seconds/sample: 
Ts = 1/Fs;
%Array for t values: samples * seconds/sample
t = [0:length(modulatingSignal)-1]*Ts;
plot(t,modulatingSignal);
title('Input modulating signal waveform');
xlabel('time / s');
ylabel('amplitude');
figure;

%f-domain representation before fftShift:
modulatingSignalF = fft(modulatingSignal);
plot(abs(modulatingSignalF));
title('Input modulating signal spectrum before fftshift');
xlabel('frequency / bins');                                                %change to hz
ylabel('amplitude');
figure;

%f-domain representation after fftshift:
modulatingSignalF =  fftshift(fft(modulatingSignal), 2);
f = (-Fs/2:Fs/length(modulatingSignal):Fs/2-Fs/length(modulatingSignal));   %by default, it plots frequency in bins, that's to change it to Hz 
n = length(modulatingSignal);
fshift = (-n/2:n/2-1)*(Fs/n); % zero-centered frequency range
plot(fshift, abs(modulatingSignalF));
title('Input modulating signal spectrum after fftshift');
xlabel('frequency / Hz');                                                
ylabel('amplitude');
figure;

%Ideal low-pass filter
Wcutoff = 4000 / (Fs/2);
[b, a] = butter(20, Wcutoff, 'low');

%Filtered signal in time domain
modulatingSignalLowFreq = filter(b, a, modulatingSignal);
plot(t, modulatingSignalLowFreq);
title('Modulating signal after lowpass filter');
xlabel('time / s');
ylabel('amplitude');                                                       
figure;

%Filtered signal in frequency domain
modulatingSignalFLowFreq = fftshift(fft(modulatingSignalLowFreq));
plot(f, abs(modulatingSignalFLowFreq) );
title('Modulating signal after lowpass filter');
xlabel('frequency / bins');
ylabel('amplitude');
figure;

%The pdf asked to ensure that there is only a small error in the filtered
%signal
Error=immse(modulatingSignal, modulatingSignalLowFreq);
display(Error);

sound( real(double(modulatingSignalLowFreq)), Fs);                         %with/out Fs?? default=8192 hz


%DSB-SC modulation...
Fc = 100000;

%pdf asked to increase the sampling frequency of filtered signal to be Fs=5*Fc before modulation process
%y = resample(x,p,q) resamples the input sequence, x, at p/q times the original sample rate. 
%resample applies an antialiasing FIR lowpass filter to x. 
%p has to be an int. Error is generated otherwise. Thus, use rat.
%R = rat(X) returns the rational fraction approximation of X to within the default tolerance,
%[N,D] = rat(___) returns two arrays, N and D, such that N./D approximates X, using any of the above syntaxes
[N,D] = rat((5*Fc)/Fs);
modulatingSignalLowFreqResampled = resample(modulatingSignalLowFreq, N, D);

t1=(0:1/500000:8.56766);   %to keep the same signal width on the time scale despite the higher sampling rate
carrierSignal = 0.5*cos(2*pi*Fc*t1) + 0.5;   %0.5's to limit the signal between 0 & 1 
% plot(t1,carrierSignal);
% title('Carrier signal')
% xlabel('time');
% ylabel('amplitude');
% figure;

%carrierSignal' = Transpose(carrierSignal). To match matrix sizes for .*
DSBSC = modulatingSignalLowFreqResampled(1:length(carrierSignal)).*carrierSignal';
plot(t1, DSBSC);
title('DSB-SC signal (resampled)')
xlabel('time');
ylabel('amplitude');
figure;

%Demodulation...
% awgn(in,snr) adds white Gaussian noise to the vector signal in. awgn(in,snr) adds white Gaussian noise to the vector signal in.
%This syntax assumes that the power of in is 0 dBW. Thus, use 'measured' o
%have the function measure the power of in before adding noise.

%Add noise...
noisyModulatedSignal_1 = awgn(DSBSC,0,'measured'); %should be the same as the modulated signal as noise power = 0
plot(t1,[DSBSC noisyModulatedSignal_1]);
legend('Original Signal','Signal with AWGN');
title('DSB-SC signal - added noise with SNR = 0');
xlabel('time');
figure;
%calculate sampling rate of noisy modulated signal
% n_N1 = length(noisyModulatedSignal_1);   %num of samples
% Fs_N1 = ( max(noisyModulatedSignal_1) - min(noisyModulatedSignal_1) ) / (n_N1 - 1);
% sound( real(double(noisyModulatedSignal_1)), Fs_N1);
%plot in noisy modulated signal in frequency domain
plot( abs( fft(noisyModulatedSignal_1)));
title('DSB-SC signal - added noise with SNR = 0');
xlabel('frequency / bins');
figure;
sound( real(double(noisyModulatedSignal_1)));

noisyModulatedSignal_2 = awgn(DSBSC,10,'measured');
plot(t1,[DSBSC noisyModulatedSignal_2]);
legend('Original Signal','Signal with AWGN');
title('DSB-SC signal - added noise with SNR = 10');
xlabel('time');
figure;
%plot in noisy modulated signal in frequency domain
plot( abs( fft(noisyModulatedSignal_2)));
title('DSB-SC signal - added noise with SNR = 10');
xlabel('frequency / bins');
figure;
sound( real(double(noisyModulatedSignal_2)));

noisyModulatedSignal_3 = awgn(DSBSC,30,'measured');
plot(t1,[DSBSC noisyModulatedSignal_3]);
legend('Original Signal','Signal with AWGN');
title('DSB-SC signal - added noise with SNR = 30');
xlabel('time');
figure;
%plot in noisy modulated signal in frequency domain
plot( abs( fft(noisyModulatedSignal_3)));
title('DSB-SC signal - added noise with SNR = 30');
xlabel('frequency / bins');
figure;
sound( real(double(noisyModulatedSignal_3)));

%Demodulation...
[A, B] = rat(48000/500000);
demodulatedNoisySignal_1 = noisyModulatedSignal_1.*carrierSignal';
% plot(t1,demodulatedNoisySignal_1);
% title('DSB-SC DEMODULATED signal - added noise with SNR = 0');
% xlabel('time');
% figure;
resampledDemodulated_1 = resample(demodulatedNoisySignal_1, A, B);
sound( real(double( resampledDemodulated_1 )));

demodulatedNoisySignal_2 = noisyModulatedSignal_2.*carrierSignal';
% plot(t1,demodulatedNoisySignal_2);
% title('DSB-SC DEMODULATED signal - added noise with SNR = 10');
% xlabel('time');
% figure;
resampledDemodulated_2 = resample(demodulatedNoisySignal_2, A, B);
sound( real(double( resampledDemodulated_2 )));

demodulatedNoisySignal_3 = noisyModulatedSignal_3.*carrierSignal';
% plot(t1,demodulatedNoisySignal_3);
% title('DSB-SC DEMODULATED signal - added noise with SNR = 30');
% xlabel('time');
% figure;
resampledDemodulated_3 = resample(demodulatedNoisySignal_3, A, B);
sound( real(double( resampledDemodulated_3 )));

filter1 = filter(b, a, resampledDemodulated_1);
plot(t, filter1);
title('Demodulated and resampled (SNR = 0)');
xlabel('time');
figure;
filter2 = filter(b, a, resampledDemodulated_2);
plot(t, filter2);
title('Demodulated and resampled (SNR = 10)');
xlabel('time');
figure;
filter3 = filter(b, a, resampledDemodulated_3);
plot(t, filter3);
title('Demodulated and resampled (SNR = 30)')
xlabel('time');
figure;

demodulatedF_1 = fftshift(fft(filter1));
plot(t, demodulatedF_1);
title('Demodulated and resampled (SNR = 0) - Fdomain')
xlabel('F / Hz');
figure;
demodulatedF_2 = fftshift(fft(filter2));
plot(t, demodulatedF_2);
title('Demodulated and resampled (SNR = 10) - Fdomain')
xlabel('F / Hz');
figure;
demodulatedF_3 = fftshift(fft(filter3));
plot(t, demodulatedF_3);
title('Demodulated and resampled (SNR = 30) - Fdomain')
xlabel('F / Hz');
figure;

%coherent detection with frequency error...
FcError = 100100;
carrierSignalError = 0.5*cos(2*pi*FcError*t1) + 0.5;
demodulatedNoisySignal_2Error = noisyModulatedSignal_2.*carrierSignalError';
plot(t1, demodulatedNoisySignal_2Error);
title('DSB-SC DEMODULATED signal using frequency error - added noise with SNR = 10');
xlabel('time');
sound( real(double( resample(demodulatedNoisySignal_2Error, A, B) )));

%coherent detection with phase error...
carrierSignalError = 0.5*cos(2*pi*Fc*t1 + degtorad(20)) + 0.5;
demodulatedNoisySignal_2Error = noisyModulatedSignal_2.*carrierSignalError';
plot(t1, demodulatedNoisySignal_2Error);
title('DSB-SC DEMODULATED signal using frequency error - added noise with SNR = 10');
xlabel('time');
sound( real(double( resample(demodulatedNoisySignal_2Error, A, B) )));


 







