clear;
ke = 200;
kc = ke / 2;
ks = 100;
nsteps = 600;

ddx = .01;
dt = ddx / (2 * 3e8);
freq_in = 7 * 1e8;

epsilon = 4;
sigma = 0.04;
epsz = 8.85419e-12;
eaf = dt * sigma / (2 * epsz * epsilon);

ex_low_m2 = 0;
ex_low_m1 = 0;

ca = ones(ks, 1);
cb = zeros(ks, 1);
ex = zeros(1, ke);
hy = zeros(1, ke);
ix = zeros(1, ke);

T = 0;

% Open file to write the data
fid = fopen('arxeio7.txt', 'w');

% Initialize figure for plotting
fig = figure;

% Initialize the structure to store frames for the movie
M(nsteps) = struct('cdata', [], 'colormap', []);

for t = 1:nsteps
    T = T + 1;

    % Update Ex field for free space
    for k = 2:ks - 1
        ex(k) = ex(k) + .5 * (hy(k - 1) - hy(k));
    end

    % Update Ex field for lossy dielectric medium
    for k = ks:ke - 1
        ca(k, 1) = (1 - eaf) / (1 + eaf);
        cb(k, 1) = .5 / (epsilon * (1 + eaf));
        ex(k) = ca(k, 1) * ex(k) + cb(k, 1) * (hy(k - 1) - hy(k));
        ex(k) = ex(k) - .0125 * ex(k); % Loss in the medium
    end

    % Source pulse
    pulse = sin(2 * pi * freq_in * dt * t); 
    ex(5) = ex(5) + pulse;

    % Apply boundary conditions
    if t > 1
        ex(1) = ex_low_m2;
        ex_low_m2 = ex_low_m1;
        ex_low_m1 = ex(2);
    end

    % Update Hy field
    for k = 1:ke - 1
        hy(k) = hy(k) + .5 * (ex(k) - ex(k + 1));
    end

    % Plot the results and capture the frame at every time step
    plot(1:ke, ex, 'LineWidth', 2);
    axis([0 ke -1.5 1.5]);
    title('Sinusoidal wave (f = 700 MHz)', 'color', 'k');

    % Add reference lines and annotations
    line('XData', [100 100], 'YData', [-1.5 1.2], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
    line('XData', [0 200], 'YData', [0 0], 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'k');

    text(140, 1, 'eps = 4', 'FontSize', 14);
    text(120, -0.8, 'conductivity = 0.04', 'FontSize', 14);
    text(20, 1.3, 'Lossy Dielectric medium at x = 100', 'FontSize', 16);

    % Display time step
    Title_ = ['t = ', num2str(t)];
    text(160, -1.2, Title_, 'Color', 'r');

    xlabel('x (FDTD cells)', 'FontSize', 20);
    ylabel('Ex (V/m)', 'FontSize', 20);

    % Capture the frame for the movie
    M(t) = getframe(fig);
end

% Close the file
fclose(fid);

% Play the movie at a frame rate of 10 frames per second
movie(fig, M, 1, 10);
