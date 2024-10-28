clear;
ke = 200;
kc = ke / 2;
nsteps = 350;
t0 = 40.0;
spread = 12;

% Initialize the electric (ex) and magnetic (hy) fields
ex = zeros(1, ke);
hy = zeros(1, ke);

% Initialize absorbing boundary conditions (ABC) for ex
ex_low_m2 = 0.;
ex_low_m1 = 0.;
ex_high_m1 = 0.;
ex_high_m2 = 0.;

% Open file to write the data
fid = fopen('arxeio1.txt', 'w');

% Initialize figure for plotting
fig = figure;

% Create a structure to hold frames for the movie
M(nsteps) = struct('cdata', [], 'colormap', []);

for t = 1:nsteps
    % Update Ex field
    for k = 2:ke - 1
        ex(k) = ex(k) + 0.5 * (hy(k - 1) - hy(k));
    end
    
    % Apply absorbing boundary conditions (ABC)
    ex(1) = ex_low_m2;
    ex_low_m2 = ex_low_m1;
    ex_low_m1 = ex(2);

    ex(ke) = ex_high_m2;
    ex_high_m2 = ex_high_m1;
    ex_high_m1 = ex(ke - 1);
    
    % Source pulse
    pulse = exp(-0.5 * ((t - t0) / spread) ^ 2);
    ex(kc) = ex(kc) + pulse;

    % Update Hy field
    for k = 1:ke - 1
        hy(k) = hy(k) + 0.5 * (ex(k) - ex(k + 1));
    end

    % Write data to file
    fprintf(fid, '%f\n', ex);

    % Plot the current state
    plot(1:ke, ex, 'LineWidth', 2);
    axis([0 ke -0.5 1.5]);
    xlim([0, 200]);
    ylim([-0.5, 1.5]);

    % Add vertical and horizontal reference lines
    line('XData', [100 100], 'YData', [-0.5 1.2], 'LineStyle', '--', ...
        'LineWidth', 1, 'Color', 'k');
    line('XData', [0 200], 'YData', [0 0], 'LineStyle', '--', ...
        'LineWidth', 0.5, 'Color', 'k');

    % Add note and title
    text(20, 1.3, 'Absorbing Boundary Conditions', 'FontSize', 12);
    title(['t = ', num2str(t)], 'Color', 'r');

    % Axis labels
    xlabel('x (FDTD cells)', 'FontSize', 20);
    ylabel('Ex (V/m)', 'FontSize', 20);
    
    % Capture the frame for the movie
    M(t) = getframe(fig);
end

% Close the file
fclose(fid);

% Play the movie at a frame rate of 10 frames per second
movie(fig, M, 1, 10);
