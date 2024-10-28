% Lorentz Medium

% Permittivity of Lorentz medium
clear;

epsr = 1;
eps1 = 1;
delta = 0.25;

eps_star = zeros;

ddx = 6e-5;
dt = ddx / 6e8;
nsteps = 1500;

N = nsteps * dt;

freq = 0.5e12;
freq3 = freq;
omega3 = 2 * pi * freq;
fs = freq * 20;
T = 1 / fs;
L = 1500;
sigma = 60;
t = (0:nsteps);
t0 = 240;

% Generate the pulse
pulse = (exp(-0.5 * (t - t0).^2 ./ (sigma)^2)) .* sin(2 * pi * freq * dt * t);

n = 2^nextpow2(L);
Z = fft(pulse, n) / L;

freq2 = 0.5e12;
omega0 = 2 * pi * freq2;
freq1 = linspace(0, 20e12, 20000);
omega = 2 * pi * freq1;

% Permittivity of Lorentz medium
eps_star = epsr + eps1 * (1 + 2 * 1i * delta * (omega / omega0) - (omega / omega0).^2).^(-1);

real_eps = real(eps_star);
imag_eps = imag(eps_star);

% Refractive index
n_index = sqrt(eps_star);
real_n = real(n_index);
imag_n = imag(n_index);

% Plot the refractive index
figure;
plot(omega / 1e12, real_n, 'b-', 'LineWidth', 2);
hold on;
plot(omega / 1e12, -imag_n, 'r-', 'LineWidth', 2);
hold off;
axis([0 10 -1 2.5]);

% Add labels and title
grid on;
title('Refractive index in Lorentz medium');
ylabel('n');
xlabel('Frequency (THz)');

% Add legend and annotations
legend('Real part of n', 'Imaginary part of n');
text(2.5, 1.5, '\delta = 0.25', 'FontSize', 12, 'Color', 'k');
text(2.5, 1.2, ['\omega_0 = ', num2str(freq2 / 1e12), ' THz'], 'FontSize', 12, 'Color', 'b');

% Add vertical line at resonance frequency
resonance_freq = omega0 / 1e12;
line('XData', [resonance_freq resonance_freq], 'YData', [-1 2.5], 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'g');

% Add horizontal reference line
line('XData', [0 10], 'YData', [1 1], 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'k');

