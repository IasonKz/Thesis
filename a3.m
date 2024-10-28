clear;
ke = 200;
kc = ke / 2;
ks = 100;  % Start of dielectric medium
nsteps = 1000;
t0 = 40.0;
spread = 12;
epsilon = 4;  % Dielectric constant

% Initialize the electric (ex) and magnetic (hy) fields
ex = zeros(1, ke);
hy = zeros(1, ke);

% Initialize absorbing boundary conditions (ABC) for ex
ex_low_m2 = 0.;
ex_low_m1 = 0.;
ex_high_m1 = 0.;
ex_high_m2 = 0.;

% Precompute the coefficients for dielectric and free space
cb = ones(1, ke) * 0.5;
cb(ks:ke) = 0.5 / epsilon;  % Dielectric region

% Open file to write the data
fid = fopen('arxeio6.txt', 'w');

% Initialize figure for plotting
fig = figure;

% Create a structure to hold frames for the movie
M(nsteps) = struct('cdata', [], 'colormap', []);

for t = 1:nsteps
    % Update Ex field in free space and dielectric medium
    for k = 2:ke - 1
        ex(k) = ex(k) + cb(k) * (hy(k - 1) - hy(k));
    end
    
    % Apply absorbing boundary conditions (ABC)
    ex(1) = ex_low_m2;
    ex_low_m2 = ex_low_m1;
    ex_low_m1 = ex(2);

    ex(ke) = ex_high_m2;
    ex_high_m2 = ex_high_m1;
    ex_high_m1 = ex(ke - 1);

    % Source pulse at x = 2
    pulse = exp(-0.5 * ((t - t0) / spread) ^ 2);
    ex(2) = pulse;

    % Update Hy field
    for k = 1:ke - 1
        hy(k) = hy(k) + 0.5 * (ex(k) - ex(k + 1));
    end

    % Write data to file
    fprintf(fid, '%f\n', ex);

    % Plot the current state
    plot(1:ke, ex, 'LineWidth', 2);

    axis([0 ke -1.5 1.5]);
    title('Gaussian Pulse at x=0', 'color', 'k');

    xlabel('x (FDTD cells)', 'FontSize', 20);
    ylabel('Ex (V/m)', 'FontSize', 20);

    % Vertical and horizontal reference lines
    line('XData', [100 100], 'YData', [-1.5 1.2], 'LineStyle', '--', ...
        'LineWidth', 1, 'Color', 'k');
    line('XData', [0 200], 'YData', [0 0], 'LineStyle', '--', ...
        'LineWidth', 0.5, 'Color', 'k');

    % Display time step and dielectric information
    text(160, -1.2, ['t = ', num2str(t)], 'Color', 'r');
    text(140, 1, 'eps = 4', 'FontSize', 14);
    text(20, 1.3, 'Dielectric medium at x = 100', 'FontSize', 16);

    % Capture the frame for the movie
    M(t) = getframe(fig);
end

% Close the file
fclose(fid);

% Play the movie at a frame rate of 10 frames per second
movie(fig, M, 1, 10);
