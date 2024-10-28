clear;

ddx = 1.5e-5;
dt = ddx / 6e8;

nsteps = 1500;

fp = 1e12;
omegap = 2 * pi * fp;

v = 0.01e12;

eps_star = zeros;
n_index = zeros;

freq = 2e12;
fs = freq * 20; % Sampling frequency
T = 1 / fs; % Unit time [fs]
spread = 30; % Pulse duration
t0 = 50; % Used to center the pulse
N = 0:nsteps; % Time base

% Generate Gaussian modulated sinusoidal pulse
pulse = (exp(-0.5 * (N - t0) .^ 2 ./ (spread) ^ 2)) .* sin(2 * pi * freq * dt * (N));

% Perform FFT
n = 2 ^ nextpow2(nsteps); % Next power of 2 for zero-padding
Z = fft(pulse, n) / nsteps;

% Frequency vector for FFT
freq_fft = fs * (0:(n/2)) / n;
omega2 = 2 * pi * freq_fft;  % Angular frequency

% Define Drude model parameters
freq1 = linspace(0, 20e12, 20000);
omega = 2 * pi * freq1;
eps_star = 1 - (omegap .^ 2) ./ (omega .^ 2 + v .* 1i .* omega);

% Extract real and imaginary parts of dielectric function
real_eps = real(eps_star);
imag_eps = imag(eps_star);

% Plot the pulse in the time domain
subplot(2,1,1);
plot(N, pulse, 'b', 'LineWidth', 2);
grid on;
ylabel('Pulse (t)');
xlabel('t(ps)');
axis([0 200 -1 1]);
title('Source, f_0 = 2 THz, \sigma = 30, time domain');

% Plot the spectrum in the frequency domain
subplot(2,1,2);
plot(omega2 / 1e12, 100 * abs(Z(1:n/2+1)), 'b', 'LineWidth', 2);
hold on;
plot(omega / 1e12, real_eps, 'g', 'LineWidth', 2);
hold off;
grid on;
axis([0 15 -1 2]);
title('Source spectrum, f_0 = 2 THz, frequency domain');

% Add annotations and vertical lines for specific frequencies
text(2, -0.5, 'Real \epsilon(\omega)', 'Color', 'g');
line('XData', [6.28319 6.28319], 'YData', [-1 2], 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'g');
text(8, -0.5, 'Pulse spectrum', 'Color', 'b');
line('XData', [12.5664 12.5664], 'YData', [-1 2], 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'b');

ylabel('Amplitude');
xlabel('\omega (THz)');
legend('Pulse spectrum', 'Real \epsilon(\omega)');
