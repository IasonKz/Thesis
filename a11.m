clear;

ke = 600;
kc = 300;
epsz = 8.8e-12;
ddx = 1.5e-5; % space step
dt = ddx / 6e8; % time step
kstart = 200;
kend = 400;
nsteps = 1300;

freq = 2e12;

% Initialize fields and variables
ex = zeros(1, ke);
dx = zeros(1, ke);
hy = zeros(1, ke);
ix = zeros(1, ke);
sx = zeros(1, ke);
sx1 = zeros(1, ke);
sx2 = zeros(1, ke);

t0 = 50.0;
spread = 30.0;

ex_low_m1 = 0.0;
ex_low_m2 = 0.0;
ex_high_m1 = 0.0;
ex_high_m2 = 0.0;

T = 0;
fig = figure;

% Initialize structure to store frames for the movie
M(nsteps) = struct('cdata', [], 'colormap', []);

fid = fopen('plasmaADE.txt', 'w');

for t = 1:nsteps
    T = T + 1;

    % Update Dx field
    for k = 2:ke
        dx(k) = dx(k) + 0.5 * (hy(k - 1) - hy(k));
    end

    % Source pulse
    x = sin(2 * pi * freq * dt * T);
    y = exp(-0.5 * ((T - t0)^2 / (spread)^2));
    pulse = x * y; % Source
    dx(5) = dx(5) + pulse;

    % Update Ex field using ADE method in the plasma region
    for k = 2:ke
        if (k >= kc && k <= kend)
            omega = 2 * pi * 1e12;
            vc = 0.01e12;
            ex(k) = (dx(k) - ix(k) - ((1 - 0.5 * vc * dt) / (1 + 0.5 * vc * dt)) * sx(k));
            ix(k) = ix(k) + ((omega^2) * dt / vc) * ex(k);
            sx(k) = ((1 - 0.5 * vc * dt) / (1 + 0.5 * vc * dt)) * sx(k) - ((omega^2) * dt / vc) * ex(k);
        else
            ex(k) = dx(k);
        end
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

    % Capture and store frames every 10 steps for smoother video
    if mod(t, 10) == 0
        plot((1:ke), ex, 'LineWidth', 2);
        grid on;
        axis([0 ke -1.2 1.2]);

        rectangle('Position', [300, -1.2, 100, 2.4], 'Curvature', 1, 'EdgeColor', 'g', 'LineWidth', 2.5);

        titlestring = 'Frequency of propagating wave = 2 THz, ADE method';
        title(titlestring, 'color', 'k');

        line('XData', [300 300], 'YData', [-1.5 1.5], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'g');
        line('XData', [400 400], 'YData', [-1.5 1.5], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'g');
        line('XData', [0 100], 'YData', [0 0], 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'k');
        line('XData', [300 400], 'YData', [1 1], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');

        text(310, 1.1, 'Plasma', 'FontSize', 14, 'Color', 'g');

        text(160, -0.8, ['t = ', num2str(t)], 'Color', 'r');

        xlabel('x( FDTD cells )');
        ylabel('Ex(V/[x])');

        % Capture the frame for the movie
        M(t) = getframe(fig);
    end
end

% Play the movie at 10 frames per second
movie(fig, M, 1, 10);

% Close file
fclose(fid);
