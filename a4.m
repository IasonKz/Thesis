clear;
ke = 200;
kc = ke / 2;
ks = 100;
nsteps = 600;

S = 1;
ddx = .01;
dt = ddx / (2 * 3e8);
freq_in = 700 * 1e6;

ex_low_m2 = 0;
ex_low_m1 = 0;

ex = zeros(1, ke);
hy = zeros(1, ke);
t0 = 40.0;

epsilon = 4;

% Open file to write the data
fid = fopen('arxeio6.txt', 'w');

% Initialize figure for plotting
fig = figure;

% Initialize the structure to store frames for the movie
M(nsteps) = struct('cdata', [], 'colormap', []);

for t = 1:nsteps
    % Update Ex field
    for k = 2:ke - 1
        cb(k) = 0.5;
    end
    for k = ks:ke - 1
        cb(k) = 0.5 / epsilon;
    end
    for k = 2:ke - 1
        ex(k) = ex(k) + cb(k) * (hy(k - 1) - hy(k));
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
        hy(k) = hy(k) + 0.5 * (ex(k) - ex(k + 1));
    end

    % Write data to file
    fprintf(fid, '%f\n', ex);

    % Plot the results
    plot(1:ke, ex, 'LineWidth', 2);
    axis([0 ke -1.5 1.5]);
    title('Sinusoidal wave (f_0 = 700 MHz)', 'color', 'k');

    % Add reference lines and annotations
    line('XData', [100 100], 'YData', [-1.5 1.2], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
    line('XData', [0 200], 'YData', [0 0], 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'k');

    text(140, 1, 'eps = 4', 'FontSize', 14);
    text(20, 1.3, 'Dielectric medium at x = 100', 'FontSize', 14);

    % Display time step
    Title_ = ['t = ', num2str(t)];
    text(160, -1.2, Title_, 'Color', 'r');

    xlabel('x (FDTD cells)', 'FontSize', 20);
    ylabel('Ex (V/m)', 'FontSize', 20);

    % Capture frame for the movie
    M(t) = getframe(fig);
end

% Close the file
fclose(fid);

% Play the movie at a frame rate of 10 frames per second
movie(fig, M, 1, 10);
