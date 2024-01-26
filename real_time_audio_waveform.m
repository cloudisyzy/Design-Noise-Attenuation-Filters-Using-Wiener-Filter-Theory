clear

[signal, Fs] = audioread('EQ2401project1data2024.wav');
% Calculate time interval for each sample
dt = 1 / Fs;
% Calculate total duration of the signal
time = 0:dt:(length(signal) / Fs) - dt;

% Create a time domain plot
figure;
hPlot = plot(time(1), signal(1));
xlim([time(1) time(end)]);
ylim([min(signal) max(signal)]);
xlabel('Time');
ylabel('Amplitude');
title('Real-time Audio Waveform');

% Set batch size for processing
batchSize = 100; % Adjust as needed

% Play the audio
sound(signal, Fs);

% Update the plot
for i = 1:batchSize:length(signal)
    % Determine the end index for the current batch
    endIndex = min(i + batchSize - 1, length(signal));
    
    % Update plot data
    set(hPlot, 'XData', time(1:endIndex), 'YData', signal(1:endIndex));
    drawnow;
    
    % Ensure update rate is synchronized with audio playback
    pause(dt * 0.79 * batchSize);
end
