clear;

ke = 600;
kc = 300;
epsz = 8.8e-12;
ddx = 6e-5;
dt = ddx / 6e8;
kstart = 300;
kend = 600;
nsteps = 1500;
epsr = 1;
eps1 = 1;
delta = 0.25;

freq = 0.5e12;

% Initialize fields and variables
ex = zeros(1, ke);
dx = zeros(1, ke);
hy = zeros(1, ke);
sx = zeros(1, ke);
sx1 = zeros(1, ke);
sx2 = zeros(1, ke);

t0 = 240.0;
spread = 60.0;

ex_low_m1 = 0.0;
ex_low_m2 = 0.0;
ex_high_m1 = 0.0;
ex_high_m2 = 0.0;

T = 0;
fig = figure;

% Initialize structure to store frames for the movie
M(nsteps) = struct('cdata', [], 'colormap', []);

fid = fopen('LorentzRC.txt', 'w');

for t = 1:nsteps
    T = T + 1;
    
    % Update Dx field
    for k = 2:ke
        dx(k) = dx(k) + 0.5 * (hy(k - 1) - hy(k));
    end
    
    % Source Pulse
    x = sin(2 * pi * freq * dt * T);
    y = exp(-0.5 * ((T - t0)^2 / (spread)^2));
    pulse = x * y; % Source
    
    dx(5) = dx(5) + pulse;
    
    % Update Ex field using Recursive Convolution in Lorentz medium
    for k = 2:ke
        if (k >= kc && k <= ke)
            omega = 2 * pi * 0.5e12;
            ex(k) = (1 / epsr) * (dx(k) - sx(k));
            
            sx(k) = ((2 - (dt * omega)^2) * exp(-delta * omega * dt)) * sx1(k) - ...
                    exp(-2 * delta * omega * dt) * sx2(k) + ...
                    ((dt^2 * omega^2 * eps1) * exp(-delta * omega * dt)) * ex(k);
                    
            sx2(k) = sx1(k);
            sx1(k) = sx(k);
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
        
        rectangle('Position', [300, -1.2, 300, 2.4], 'Curvature', 1, 'EdgeColor', 'g', 'LineWidth', 2.5);
        
        titlestring = 'Frequency of propagating wave = 0.5 THz';
        title(titlestring, 'color', 'k');
        
        line('XData', [300 300], 'YData', [-1.5 1.5], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'g');
        
        text(320, 1.05, 'Lorentz Medium', 'FontSize', 14, 'Color', 'g');
        text(320, 0.9, 'freq = 0.5 THz', 'FontSize', 14, 'Color', 'g');
        
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

