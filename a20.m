
clear;

ddx = 2.5e-5;
dt = ddx / 6e8;

c0 = 3e8;
ddy = ddx;

% Drude medium properties
fp = 1e12;
omegap = 2 * pi * fp;
v = 0.01e12;

nsteps = 500;

IE = 500;
JE = 500;

freq = 1.2e12;

dz = zeros(IE, JE);
ez = zeros(IE, JE);
hx = zeros(IE, JE);
hy = zeros(IE, JE);
iz = zeros(IE, JE);
sz = zeros(IE, JE);

ga = ones(IE, JE);

ic = 100;
jc = 300;

t0 = 100.0;
spread = 30.0;
T = 0;

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
    dz(250, 150) = pulse;
    
    % Update Ez field in Drude medium
    for i = 1:IE
        for j = 1:JE
            if (j >= 200 && j <= 250)
                ez(i, j) = (dz(i, j) - iz(i, j) - ((1 - 0.5 * v * dt) / (1 + 0.5 * v * dt)) * sz(i, j));
                iz(i, j) = iz(i, j) + ((omegap^2) * dt / v) * ez(i, j);
                sz(i, j) = ((1 - 0.5 * v * dt) / (1 + 0.5 * v * dt)) * sz(i, j) - ((omegap^2) * dt / v) * ez(i, j);
            else
                ez(i, j) = ga(i, j) .* dz(i, j);
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

    % Apply Mur boundary conditions
    % x = 0:
    for j = 2:JE-1
        ez(2, j) = -ez(3, j) + M * (ez(3, j) + ez(2, j)) + N1 * (ez(2, j) + ez(3, j));
    end
    
    % x = IE:
    for j = 2:JE-1
        ez(IE, j) = -ez(IE-1, j) + M * (ez(IE-1, j) + ez(IE, j)) + N1 * (ez(IE-1, j) + ez(IE, j));
    end
    
    % y = 0:
    for i = 2:IE-1
        ez(i, 2) = -ez(i, 3) + MM * (ez(i, 3) + ez(i, 2)) + N2 * (ez(i, 3) + ez(i, 2));
    end
    
    % y = JE:
    for i = 2:IE-1
        ez(i, JE) = -ez(i, JE-1) + MM * (ez(i, JE-1) + ez(i, JE)) + N2 * (ez(i, JE-1) + ez(i, JE));
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
    
    % Εμφάνιση του πεδίου κάθε 10 βήματα
    if mod(n, 10) == 0
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
        text(350, 150, 0.3, Title_, 'Color', 'k');
        Title_ = {'f_0 =', freq};
        text(350, 350, 0.3, Title_, 'Color', 'k');
        Title_ = {'f_p_l_a_s_m_a =', omegap / (2 * pi)};
        text(67, 100, 0.3, Title_, 'Color', 'k');
        rectangle('Position', [200, 0, 50, 500], 'Curvature', 1, 'EdgeColor', 'k', 'LineWidth', 1)
        title('Drude medium')

        % Plot 3D surface of Ez field
        subplot(2, 1, 2)
        surf(1:1:IE, 1:1:JE, ez)
        shading interp
        colorbar
        colormap(jet)
        axis([0 500 0 500 -0.15 0.15]);
        zlabel('Ez');
        xlabel('x');
        ylabel('y');
        Title_ = {'t=', n};
        text(270, 200, 0.2, Title_, 'Color', 'k');
        rectangle('Position', [200, 0, 50, 500], 'Curvature', 1, 'EdgeColor', 'k', 'LineWidth', 1)
        
        % Προσθήκη καθυστέρησης ώστε να δημιουργηθεί η ψευδαίσθηση του βίντεο
        pause(0.1); % καθυστέρηση 0.1 δευτερολέπτων για 10 καρέ/δευτερόλεπτο
    end
end
