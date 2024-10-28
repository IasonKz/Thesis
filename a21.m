% Source + spectrum (time-domain, freq-domain) Metamaterial Debye medium

clear;

ddx = 4.3e-5; % space step
dt = ddx / 6e8; % time step
nsteps = 1500;

fp = 1e12;
omegap = 2 * pi * fp;
v = 0.001e12;

eps_star = zeros;
n_index = zeros;

N = nsteps;

freq = 0.7e12;
freq3 = freq;
omega3 = 2 * pi * freq;
fs = freq * 20; % Sampling frequency
T = 1 / fs; % Unit time [fs]
L = 1000; % Length of signal
spread = 30; % Pulse duration
N = (0:nsteps); % Time base
t0 = 100; % Used to center the pulse

% Generate the sinusoidal pulse
pulse = sin(2 * pi * freq * dt * (N - t0));

% Compute the FFT
n = 2^nextpow2(L);
Z = fft(pulse, n) / L;

% Frequency and omega values
freq1 = linspace(0, 20e12, 20000);
omega = 2 * pi * freq1;

% Calculate permittivity in Debye medium
eps_star = 1 - (omegap.^2) ./ (omega.^2 + v .* 1i .* omega);

% Extract real and imaginary parts of permittivity and refractive index
real_eps = real(eps_star);
imag_eps = imag(eps_star);
n_index = sqrt(eps_star);
real_n = real(n_index);
imag_n = imag(n_index);

% Frequency vector for FFT
freq = 0.5 * fs * linspace(0, 1, n / 2 + 1);
omega2 = 2 * pi * freq;

% Plot results
subplot(2, 1, 1);
plot(N, pulse, 'b', 'LineWidth', 2);
grid on;
ylabel('Pulse (t)');
xlabel('t (ps)');
axis([0 200 -1 1]);
title('Sinusoidal source, f_0 = 0.7 THz, time domain');

subplot(2, 1, 2);
plot(omega2 / 1e12, 5 * abs(Z(1:n / 2 + 1)), 'g-', 'LineWidth', 2);
axis([0 12 -2.5 2.5]);
grid on;
title('Source spectrum, f_0 = 0.7 THz, frequency domain');
Title_1 = {'\omega_0 =', omega3};
text(1, -1, Title_1, 'Color', 'g');
line('XData', [4.4 4.4], 'YData', [-30 25], 'LineStyle', '- -', 'LineWidth', 0.5, 'Color', 'g');

ylabel('Pulse (\omega)');
xlabel('\omega (THz)');

