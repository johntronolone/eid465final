%% Part I: Acoustic Measurements
% Helmuth Rosales and John Tronolone

% constants
fs = 44100;

%% SET UP THE DAQ DEVICE

% https://www.mathworks.com/help/daq/examples/acquire-data-using-ni-devices.html
% https://www.mathworks.com/help/daq/examples/generate-signals-on-ni-devices-that-output-voltage.html


% get conneceted devices
devices = daq.getDevices;

% create a session associated with a vendor
s = daq.createSession('ni');

% add analog output channel
addAnalogOutputChannel(s, devices.ID, 'ao0', 'Voltage');

% add analog input channel
addAnalogInputChannel(s, devices.ID, 'ai0', 'IEPE');

% set the session rate (scans/sec)
s.Rate = fs;

% set scan duration
%s.DurationInSeconds = 2;



%% SET UP SIGNAL
% Create your sine sweep and averaging


% chirp parameters
f_0 = 1; % Hz
f_1 = 1e4; % Hz
%T = s.DurationInSeconds;
T_p = 2;


% generate logarithmic sinesweep
k = (f_1/f_0).^(1/T_p);
t = (0:(fs*2-1))/fs;
phase_init = 0;
sinesweep = sin(phase_init + 2*pi*f_0*(((k.^t)-1)/log(k)));

sinesweep = repmat([sinesweep, zeros(1, length(sinesweep))], 1, 5);

% something is up with my spectrogram function but the sweep sounds ok
%soundsc(sinesweep, fs);

% window_length =  256;
% overlap = 128;
% fft_length = 1024;
% 
% spect = spectrogram_plus(sinesweep, fs, fft_length, window_length, overlap);
% time_vector = linspace(0, length(sinesweep)/fs,  floor(length(sinesweep)/(window_length-overlap)));
% freq_vec = linspace(0, fs/2, fft_length/2);
% 
% figure(102)
% image(time_vector, abs(freq_vec-6000), 10*log10(abs(spect)), 'CDataMapping','scaled')
% set(gca,'YDir','normal')
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% title('Chirp')
% colorbar


%% DO THE MEASUREMENT
% Interfacing with the 4431

% use the startForeground function to start the analog output operation and
% block MATLAB execution until all data is generated.
queueOutputData(s, sinesweep');
[output_data,time] = s.startForeground;

%% PLOT (MAKING SURE THAT YOU ARE COLLECTING DATA)
figure(1)
plot(time, output_data);

%% ANALYSIS
% Calculate the Impulse Response (IR)

fft_size = 1024;

X_f = fft(sinesweep', fft_size);
Y_f = fft(output_data, fft_size);

%X_f = fft_plus(sinesweep', fs, fft_size);%, fs, fft_size);
%Y_f = fft_plus(output_data, fs, fft_size);%, fs, fft_size);

%%
impulseResponse_f = Y_f./X_f;
impulseResponse = ifft(impulseResponse_f, fft_size);


%%

% octave band filtering:
%   Center frequencies of filter are:
%   125 Hz, 250 Hz, 500 Hz, 1 kHz, 2 kHz, 4 kHz

% Plot the Octave Band Filtered IR and FRF

% plot the impulse response filtered by octave band filters
w = [31.25, 62.5, 125, 250, 500, 1e3, 2e3, 4e3];
w2 = [w*0.65; w*1.35];

y = zeros(length(w), length(impulseResponse));

figure(2)
for i = 1:length(w)
    y(i, :) = bandpass(impulseResponse, w2(:,i), fs);
    plot(1:length(y(i,:)), y(i,:))
    hold on
end
hold off

%% plot the bandpass filters (frequency response function)
lenFFT =1024;
filter_shape = zeros(length(w), lenFFT);

figure(3)
for i = 1:length(w)
    filter_shape(i,:) = fft_plus(bandpass(sinesweep, w2(:,i), fs), fs, lenFFT);
    plot(1:lenFFT, 20*log10(abs(filter_shape(i,:))))
    hold on
end
hold off

%% Calculation of reverb time

% plot Schroeder curves

for i = 1:5
    dataToAvg(i,:) = output_data((i-1)*length(output_data)/5+1:i*(length(output_data)/5));
end
avgdData = sum(dataToAvg, 1)/5;

%%
figure(4)
plot(1:length(avgdData), 20*log10(avgdData))
hold on
for i = 1:10
    
    Tc = i*0.01;
    dt = 1/fs;

    NUM = [0 1];
    DEN = [Tc/dt 1-Tc/dt];

    %transfer_function = tf(NUM, DEN);

    exponential_filtered = filter(NUM, DEN, avgdData.^2);
    plot(1:length(exponential_filtered), 10*log10(abs(exponential_filtered)))
end
hold off
