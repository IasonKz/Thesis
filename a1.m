clear;
% Define the parameters of the problem
ke = 200;
kc = ke / 2;
nsteps = 500;
t0 = 40.0;
spread = 12;

% Create the fields
ex = zeros(1, ke);
hy = zeros(1, ke);

% Open file to write the data
fid = fopen('arxeio1.txt', 'w');

% Initialize figure for plotting
fig = figure;

% Initialize structure to store frames for the movie
M(nsteps) = struct('cdata', [], 'colormap', []);

for t = 1:nsteps
    % Update Ex field
    for k = 2:ke - 1
        ex(k) = ex(k) + 0.5 * (hy(k - 1) - hy(k));
    end

    % Source pulse
    pulse = exp(-0.5 * ((t - t0) / spread) ^ 2);
    ex(kc) = ex(kc) + pulse;

    % Update Hy field
    for k = 1:ke - 1
        hy(k) = hy(k) + 0.5 * (ex(k) - ex(k + 1));
    end

    % Write data to file
    fprintf(fid, '%f\n', ex);

    % Plot the results
    plot(1:ke, ex, 'LineWidth', 1.5);
    axis([0 ke -1.5 1.5]);

    % Title and labels
    Title_ = ['t = ', num2str(t)];
    text(170, 1.2, Title_, 'Color', 'k');

    xlabel('y (FDTD cells)');
    ylabel('E_x (V/m)');

    % Capture the frame for the movie
    M(t) = getframe(fig);
end

% Close the file
fclose(fid);

% Play the movie at a frame rate of 10 frames per second
movie(fig, M, 1, 10);

