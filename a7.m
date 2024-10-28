clear;
ke = 200;
kc = ke / 2;
ks = 100;
pi = 3.141559;
epsz = 8.8e-12;
ddx = 0.01;
dt = ddx / 6e8;

epsr = 2.0;
sigma = 0.01;
chi1 = 2.0;
tau = 0.001e-6;
del_exp = exp(-dt / tau);

nsteps = 600;
t0 = 50;
spread = 10;

ex_low_m2 = 0.0;
ex_low_m1 = 0.0;
ex_high_m1 = 0.0;
ex_high_m2 = 0.0;

ga = ones(ks, 1);
gb = zeros(ks, 1);
gc = zeros(ks, 1);
dx = zeros(ke, 1);
ex = zeros(1, ke);
hy = zeros(1, ke);
ix = zeros(1, ke);
sx = zeros(ke, 1);

T = 0;

% Initialize video writer if needed
fig = figure;
M(nsteps) = struct('cdata', [], 'colormap', []);

fid = fopen('DebyeRC.txt', 'w');

% Initialize material parameters
for k = ks:ke - 1
    ga(k) = 1 / (epsr + (sigma * dt / epsz) + (chi1 * dt / tau));
    gb(k) = sigma * dt / epsz;
    gc(k) = chi1 * dt / tau;
end

for t = 1:nsteps
    T = T + 1;

    % Update electric field
    for k = 2:ke
        dx(k) = dx(k) + 0.5 * (hy(k - 1) - hy(k));
    end

    % Apply Gaussian pulse source
    pulse = exp(-0.5 * ((t0 - T) / spread)^2);
    dx(5) = dx(5) + pulse;

    % Update Ex field using RC method for Debye medium
    for k = 2:ke - 1
        ex(k) = ga(k) * (dx(k) - ix(k) - sx(k));
        ix(k) = ix(k) + gb(k) * ex(k);
        sx(k) = del_exp * sx(k) + gc(k) * ex(k);
    end

    % Update boundary conditions for Ex
    ex(1) = ex_low_m2;
    ex_low_m2 = ex_low_m1;
    ex_low_m1 = ex(2);

    ex(ke) = ex_high_m2;
    ex_high_m2 = ex_high_m1;
    ex_high_m1 = ex(ke - 1);

    % Update Hy field
    for k = 1:ke - 1
        hy(k) = hy(k) + 0.5 * (ex(k) - ex(k + 1));
    end

    % Plot and capture each frame
    if mod(t, 5) == 0 % Capture every 5 steps for a smoother video
        plot(1:ke, ex, 'LineWidth', 2);
        grid on;
        axis([0 ke -1 1.2]);

        titlestring = 'Gaussian Pulse, RC method';
        title(titlestring, 'color', 'k');

        line('XData', [100 100], 'YData', [-1.5 1.5], 'LineStyle', '--', ...
            'LineWidth', 1, 'Color', 'k');
        line('XData', [0 100], 'YData', [0 0], 'LineStyle', '--', ...
            'LineWidth', 0.5, 'Color', 'k');
        line('XData', [100 200], 'YData', [1 1], 'LineStyle', '--', ...
            'LineWidth', 0.5, 'Color', 'k');

        text(10, 0.8, 'Debye Medium at x = 100', 'FontSize', 12);
        text(160, -0.8, ['t = ', num2str(t)], 'Color', 'r');

        xlabel('x( FDTD cells )');
        ylabel('Ex(V/[x])');

        % Capture frame for the movie
        M(t) = getframe(fig);
    end
end

% Close file
fclose(fid);

% Play the movie at 10 frames per second
movie(fig, M, 1, 10);
