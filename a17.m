clear;

IE = 100;
JE = 100;

ga = ones(IE, JE);
dz = zeros(IE, JE);
ez = zeros(IE, JE);
hx = zeros(IE, JE);
hy = zeros(IE, JE);

ic = IE / 2;
jc = JE / 2;

t0 = 20.0;
spread = 6.0;
nsteps = 100;

fig = figure;

% Initialize structure to store frames for the movie
M(nsteps) = struct('cdata', [], 'colormap', []);

for n = 1:nsteps

    % Update Dz field
    for j = 2:JE-1
        for i = 2:IE-1
            dz(i, j) = dz(i, j) + 0.5 * (hy(i, j) - hy(i-1, j) - hx(i, j) + hx(i, j-1));
        end
    end
    
    % Source pulse
    pulse = exp(-0.5 * ((t0 - n) / spread)^2);
    dz(ic, jc) = pulse;
    
    % Update Ez field
    for j = 2:JE
        for i = 2:IE
            ez(i, j) = ga(i, j) .* dz(i, j);
        end
    end
    
    % Update Hx field
    for j = 2:JE-1
        for i = 2:IE-1
            hx(i, j) = hx(i, j) + 0.5 * (ez(i, j) - ez(i, j+1));
        end
    end
    
    % Update Hy field
    for j = 2:JE-1
        for i = 2:IE-1
            hy(i, j) = hy(i, j) + 0.5 * (ez(i+1, j) - ez(i, j));
        end
    end
    
    % Capture frames every 10 steps for a smoother movie
    if mod(n, 10) == 0
        % Plot 3D mesh of Ez field
        subplot(2, 1, 1)
        mesh(1:1:IE, 1:1:JE, ez)
        shading interp
        title('Gaussian pulse propagating in free space')
        axis([0 100 0 100 -0 1]);
        zlabel('Ez');
        xlabel('x');
        ylabel('y');
        Title_ = {'t=', n};
        text(37, 7, 1, Title_, 'Color', 'k');
        
        % Plot 2D Ez field as an image
        subplot(2, 1, 2)
        imagesc((1:1:IE), (1:1:JE), ez, [-0 1])
        shading interp
        colorbar
        axis([0 100 0 100]);
        xlabel('x');
        ylabel('y');
        Title_ = {'t=', n};
        text(27, 27, Title_, 'Color', 'k');
        
        % Capture the frame for movie
        M(n) = getframe(fig);
    end
end

% Play the movie at 10 frames per second
movie(fig, M, 1, 10);

