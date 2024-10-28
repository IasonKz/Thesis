clear;
close all;
clc;
c = 3e8;

freq = 100e6;  % Frequency of the sinusoidal component

ddx = 0.075;
dt = ddx / 6e8;

fsamp = freq * 80; % Sampling frequency
T = 1 / fsamp;     % Time period
L = 451;           % Length of the signal
sigma = 60.0;      % Pulse duration
t = (0:450);       % Time base
t0 = 240.0;

% Generate the Gaussian modulated sinusoidal pulse
E = (exp(-0.5 * (t - t0) .^ 2 / (sigma) ^ 2)) .* sin(2 * pi * freq * dt * (t));

% Plot the time domain signal
subplot(2,1,1);
plot(t, real(E), 'b');
title(['Gaussian Pulse \sigma =', num2str(sigma), ' s']);
xlabel('Time (fs)');
ylabel('Amplitude');
grid on;

% Perform FFT
NFFT = 2 ^ nextpow2(L); % Next power of 2 for zero-padding
X = fft(E, NFFT) / L;   % Compute the FFT and normalize

% Define the frequency vector for the FFT
freq_fft = fsamp * (0:(NFFT/2)) / NFFT;  % Frequency vector for plotting
omega2 = 2 * pi * freq_fft;              % Angular frequency (rad/s)

% Plot the magnitude of the FFT
subplot(2,1,2);
plot(omega2, abs(X(1:NFFT / 2 + 1)), 'g--', 'LineWidth', 1.5);
title('Magnitude of FFT');
xlabel('Frequency (rad/s)');
ylabel('Magnitude |X(f)|');
xlim([0 10e8]);
grid on;

% Mark a specific frequency with a vertical line
line('XData', [2 * pi * freq 2 * pi * freq], 'YData', [-0.1 max(abs(X))], 'LineStyle', '--', ...
    'LineWidth', .5, 'Color', 'k');
