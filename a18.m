clear;

IE = 100;
JE = 100;

ddx = 4.3e-5;
dt = ddx / 6e8;
c0 = 3e8;
ddy = ddx;
dt1 = dt;

freq = 0.7e12;
ga = ones(IE, JE);
dz = zeros(IE, JE);
ez = zeros(IE, JE);
hx = zeros(IE, JE);
hy = zeros(IE, JE);

ez2 = zeros(IE, JE);
ez1 = zeros(IE, JE);

ic = IE / 2;
jc = JE / 2;

t0 = 20.0;
spread = 6.0;
T = 0;
nsteps = 250;

for n = 1:nsteps
    ez2 = ez1; 
    ez1 = ez;

    % Update Dz field
    for j = 2:JE-1
        for i = 2:IE-1
            dz(i, j) = dz(i, j) + 0.5 * (hy(i, j) - hy(i-1, j) - hx(i, j) + hx(i, j-1));
        end
    end
    
    % Source Pulse
    pulse = sin(2 * pi * freq * dt * (n - t0));
    dz(ic, jc) = pulse;
    
    % Update Ez field
    for i = 2:IE-1
        for j = 2:JE-1
            ez(i, j) = ga(i, j) .* dz(i, j);
        end
    end
    
    % Mur's absorbing boundary conditions
    M1 = c0 * dt - ddx;
    M2 = c0 * dt + ddx;
    M = M1 / M2;
    N1 = (2 * ddx) / M2;

    MM1 = c0 * dt - ddy;
    MM2 = c0 * dt + ddy;
    MM = MM1 / MM2;
    N2 = (2 * ddy) / MM2;

    % Apply Mur boundary conditions
    % x = 0:
    for j = 2:JE-1
        ez(2, j) = -ez2(3, j) + M * (ez(3, j) + ez2(2, j)) + N1 * (ez1(2, j) + ez1(3, j));
    end
    
    % x = h:
    for j = 2:JE-1
        ez(IE, j) = -ez2(IE-1, j) + M * (ez(IE-1, j) + ez2(IE, j)) + N1 * (ez1(IE-1, j) + ez1(IE, j));
    end
    
    % y = 0:
    for i = 2:IE-1
        ez(i, 2) = -ez2(i, 3) + MM * (ez(i, 3) + ez2(i, 2)) + N2 * (ez1(i, 3) + ez1(i, 2));
    end
    
    % y = h:
    for i = 2:IE-1
        ez(i, JE) = -ez2(i, JE-1) + MM * (ez(i, JE-1) + ez2(i, JE)) + N2 * (ez1(i, JE-1) + ez1(i, JE));
    end
    
    % Update Hx and Hy fields
    for j = 2:JE-1
        for i = 2:IE-1
            hx(i, j) = hx(i, j) + 0.5 * (ez(i, j) - ez(i, j+1));
        end
    end
    
    for j = 2:JE-1
        for i = 2:IE-1
            hy(i, j) = hy(i, j) + 0.5 * (ez(i+1, j) - ez(i, j));
        end
    end
    
    % Οπτικοποίηση κάθε 10 βήματα
    if mod(n, 10) == 0
        % Plot 3D mesh of Ez field
        subplot(2, 1, 1)
        mesh(1:1:IE, 1:1:JE, ez)
        shading interp
        title('Sinusoidal wave propagating in free space, 1st order Mur ABC')
        axis([0 100 0 100 -0.5 1]);
        zlabel('Ez');
        xlabel('x');
        ylabel('y');
        Title_ = {'t=', n};
        text(37, 7, 1, Title_, 'Color', 'k');
        
        % Plot 2D Ez field as an image
        subplot(2, 1, 2)
        imagesc((1:1:IE), (1:1:JE), ez, [-1 1])
        shading interp
        colorbar
        axis([0 100 0 100]);
        zlabel('Ez');
        xlabel('x');
        ylabel('y');
        Title_ = {'t=', n};
        text(27, 27, Title_, 'Color', 'k');
        
        % Καθυστέρηση για ψευδαίσθηση βίντεο
        pause(0.1);
    end
end
