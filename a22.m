clear;

ddx = 4.3e-5;
dt = ddx / 6e8;
c0 = 3e8;
ddy = ddx;
dt1 = dt;

% Drude medium properties
fp = 1e12;
omegap = 2 * pi * fp;
v = 0.001e12;

nsteps = 800;

IE = 200;
JE = 200;

freq = 0.7e12;

dz = zeros(IE, JE);
ez = zeros(IE, JE);
hx = zeros(IE, JE);
hy = zeros(IE, JE);
iz = zeros(IE, JE);
sz = zeros(IE, JE);

bx = zeros(IE, JE);
by = zeros(IE, JE);
kx = zeros(IE, JE);
jx = zeros(IE, JE);
ky = zeros(IE, JE);
jy = zeros(IE, JE);
ga = ones(IE, JE);

ez2 = zeros(IE, JE);
ez1 = zeros(IE, JE);

t0 = 50.0;
T = 0;

for n = 1:nsteps
    ez2 = ez1; ez1 = ez;
    
    % Update Dz field
    for i = 2:IE-1
        for j = 2:JE-1
            dz(i,j) = dz(i,j) + 0.5 * (hy(i,j) - hy(i-1,j) - hx(i,j) + hx(i,j-1));
        end
    end
    
    % Sinusoidal source pulse
    x = sin(2 * pi * freq * dt * (n - t0));
    pulse = x;
    dz(100, 70) = pulse;
    
    % Update Ez field
    for i = 2:IE-1
        for j = 2:JE-1
            if (j >= 80 && j <= 100)
                ez(i,j) = (dz(i,j) - iz(i,j) - ((1 - 0.5 * v * dt) / (1 + 0.5 * v * dt)) * sz(i,j));
                iz(i,j) = iz(i,j) + ((omegap^2) * dt1 / v) * ez(i,j);
                sz(i,j) = ((1 - 0.5 * v * dt) / (1 + 0.5 * v * dt)) * sz(i,j) - ((omegap^2) * dt1 / v) * ez(i,j);
            else
                ez(i,j) = ga(i,j) .* dz(i,j);
            end
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

    % Mur x = 0:
    for j = 2:JE-1
        ez(2,j) = -ez2(3,j) + M * (ez(3,j) + ez2(2,j)) + N1 * (ez1(2,j) + ez1(3,j));
    end
    
    % Mur x = h:
    for j = 2:JE-1
        ez(IE,j) = -ez2(IE-1,j) + M * (ez(IE-1,j) + ez2(IE,j)) + N1 * (ez1(IE-1,j) + ez1(IE,j));
    end
    
    % Mur y = 0:
    for i = 2:IE-1
        ez(i,2) = -ez2(i,3) + MM * (ez(i,3) + ez2(i,2)) + N2 * (ez1(i,3) + ez1(i,2));
    end
    
    % Mur y = h:
    for i = 2:IE-1
        ez(i,JE) = -ez2(i,JE-1) + MM * (ez(i,JE-1) + ez2(i,JE)) + N2 * (ez1(i,JE-1) + ez1(i,JE));
    end
    
    % Update Hx and Hy fields
    for i = 2:IE-1
        for j = 2:JE-1
            bx(i,j) = bx(i,j) + 0.5 * (ez(i,j) - ez(i,j+1));
            if (j >= 80 && j <= 100)
                hx(i,j) = (bx(i,j) - kx(i,j) - ((1 - 0.5 * v * dt) / (1 + 0.5 * v * dt)) * jx(i,j));
                kx(i,j) = kx(i,j) + ((omegap^2) * dt1 / v) * hx(i,j);
                jx(i,j) = ((1 - 0.5 * v * dt) / (1 + 0.5 * v * dt)) * jx(i,j) - ((omegap^2) * dt1 / v) * hx(i,j);
            else
                hx(i,j) = hx(i,j) + 0.5 * (ez(i,j) - ez(i,j+1));
            end
        end
    end
    
    for i = 2:IE-1
        for j = 2:JE-1
            by(i,j) = by(i,j) + 0.5 * (ez(i+1,j) - ez(i,j));
            if (j >= 80 && j <= 100)
                hy(i,j) = (by(i,j) - ky(i,j) - ((1 - 0.5 * v * dt) / (1 + 0.5 * v * dt)) * jy(i,j));
                ky(i,j) = ky(i,j) + ((omegap^2) * dt1 / v) * hy(i,j);
                jy(i,j) = ((1 - 0.5 * v * dt) / (1 + 0.5 * v * dt)) * jy(i,j) - ((omegap^2) * dt1 / v) * hy(i,j);
            else
                hy(i,j) = hy(i,j) + 0.5 * (ez(i+1,j) - ez(i,j));
            end
        end
    end
    
    % Οπτικοποίηση κάθε 10 βήματα
    if mod(n, 10) == 0
        subplot(2, 1, 1)
        imagesc(1:1:IE, 1:1:JE, ez)
        shading interp
        colorbar
        colormap(jet)
        axis([0 200 0 200])
        caxis([-0.4 0.4])
        zlabel('Ez');
        xlabel('x');
        ylabel('y');
        Title_ = {'t=', n};
        text(180, 30, 0.3, Title_, 'Color', 'k');
        text(10, 50, 'E_z', 'Color', 'k');
        rectangle('Position', [80, -10, 20, 220], 'Curvature', 1, 'EdgeColor', 'k', 'LineWidth', 1.5)
        title('Metamaterial, sinusoidal source f_0 = 0.7 THz')

        subplot(2, 1, 2)
        surf(1:1:IE, 1:1:JE, ez)
        shading interp
        colorbar
        colormap(jet)
        axis([0 200 0 200 -0.4 0.4])
        caxis([-0.4 0.4])
        zlabel('Ez');
        xlabel('x');
        ylabel('y');
        Title_ = {'t=', n};
        text(180, 30, 0.3, Title_, 'Color', 'k');
        rectangle('Position', [80, -10, 20, 220], 'Curvature', 1, 'EdgeColor', 'k', 'LineWidth', 1.5)
        
        % Προσθήκη καθυστέρησης για τη δημιουργία της ψευδαίσθησης βίντεο
        pause(0.1); % καθυστέρηση 0.1 δευτερολέπτων για 10 καρέ/δευτερόλεπτο
    end
end

