clear;

ddx = 0.3e-5;
dt = ddx / 6e8;
epsr = 1;
eps1 = 1;
delta = 0.5;

nsteps = 500;

IE = 500;
JE = 500;

freq = 6e12;

ga = ones(IE, JE);
dz = zeros(IE, JE);
ez = zeros(IE, JE);
hx = zeros(IE, JE);
hy = zeros(IE, JE);

sz = zeros(IE, JE);
sz1 = zeros(IE, JE);
sz2 = zeros(IE, JE);

ic = 100;
jc = 300;

t0 = 100.0;
spread = 30.0;
T = 0;

% Initialize video recording
v = VideoWriter('2d_free_space_medium.avi');
v.FrameRate = 10;
open(v)

for n = 3:nsteps

    % Update Dz field
    for j = 2:JE-1
        for i = 2:IE-1
            dz(i, j) = dz(i, j) + 0.5 * (hy(i, j) - hy(i-1, j) - hx(i, j) + hx(i, j-1));
        end
    end
    
    % Source pulse
    x = sin(2 * pi * freq * dt * (n - t0));
    y = exp(-0.5 * ((n - t0)^2 / spread^2));
    pulse = x .* y;
    dz(250, 250) = pulse;
    
    % Update Ez field in Lorentz medium
    for i = 1:IE
        for j = 1:JE
            if j >= 300 && j <= 350
                omega = 2 * pi * 2.8e12;
                ez(i, j) = dz(i, j) - sz(i, j);
                
                sz(i, j) = ((2 - (dt * omega)^2) / (1 + dt * delta * omega)) * sz1(i, j) ...
                    - ((1 - dt * delta * omega) / (1 + dt * delta * omega)) * sz2(i, j) ...
                    + ((dt^2 * omega^2 * eps1) / (1 + delta * omega * dt)) * ez(i, j);
                    
                sz2(i, j) = sz1(i, j);
                sz1(i, j) = sz(i, j);
            else
                ez(i, j) = ga(i, j) .* dz(i, j);
            end
        end
    end
    
    % Update Hx and Hy fields
    for j = 2:JE-2
        for i = 2:IE-2
            hx(i, j) = hx(i, j) + 0.5 * (ez(i, j) - ez(i, j + 1));
        end
    end
    
    for j = 2:JE-2
        for i = 2:IE-2
            hy(i, j) = hy(i, j) + 0.5 * (ez(i + 1, j) - ez(i, j));
        end
    end
    
    % Visualization and recording at step n = 200 or every 10 steps
    if mod(n, 10) == 0 || n == 200
        % Plot 2D image of Ez field
        subplot(2, 1, 1)
        imagesc(1:1:IE, 1:1:JE, ez)
        shading interp
        colorbar
        colormap(jet)
        axis([0 500 0 500])
        zlabel('Ez');
        xlabel('x');
        ylabel('y');
        Title_ = {'t=', n};
        text(370, 150, 0.8, Title_, 'Color', 'k');
        Title_ = {'f_0 =', freq};
        text(370, 40, 0.3, Title_, 'Color', 'k');
        Title_ = {'f_m_e_d_i_u_m =', omega / (2 * pi)};
        text(420, 350, 0.3, Title_, 'Color', 'k');
        text(420, 420, 'e_1 = 10', 'Color', 'k');
        rectangle('Position', [300, 0, 50, 500], 'Curvature', 1, 'EdgeColor', 'k', 'LineWidth', 1.5)
        title('Propagating wave in Lorentz medium')
        
        % Plot 3D surface of Ez field
        subplot(2, 1, 2)
        surf(1:1:IE, 1:1:JE, ez)
        shading interp
        colorbar
        colormap(jet)
        axis([0 500 0 500 -0.2 0.2]);
        zlabel('Ez');
        xlabel('x');
        ylabel('y');
        Title_ = {'t=', n};
        text(27, 250, 0.1, Title_, 'Color', 'k');
        rectangle('Position', [300, 0, 50, 500], 'Curvature', 1, 'EdgeColor', 'k', 'LineWidth', 1)
        
        M = getframe;
        writeVideo(v, M); % Save the frame
    end
end

% Close the video recording
close(v);

