clear;

% Define grid size
IE = 1000;  % Διάσταση x (1000 κελιά)
JE = 1000;  % Διάσταση y (1000 κελιά)

ddx = 2.5e-5;
dt = ddx / 6e8;

c0 = 3e8;
ddy = ddx;

% Dielectric medium properties
n_dielec = 2;  % Δείκτης διάθλασης = 2
epsilon = n_dielec^2;  % Σχετική διηλεκτρική σταθερά = 4

nsteps = 3000;

freq = 1.2e10;

dz = zeros(IE, JE);
ez = zeros(IE, JE);
ez1 = zeros(IE, JE);  % Για την αποθήκευση της προηγούμενης χρονικής στιγμής
hx = zeros(IE, JE);
hy = zeros(IE, JE);

% Parameters for the circles
radius = 40;
num_circles = 10;
spacing = 6;
circle_spacing = 2 * radius + spacing;
start_x = 50 + 6 + radius;
center_y = JE / 2;

% Define the dielectric region with circles
ga = ones(IE, JE);

for k = 0:num_circles-1
    center_x = start_x + k * circle_spacing;
    for i = 1:IE
        for j = 1:JE
            distance = sqrt((i - center_x)^2 + (j - center_y)^2);
            if distance <= radius
                ga(i, j) = 1 / epsilon;
            end
        end
    end
end

% Χειροκίνητη θέση της πηγής: x = 10, y = 500
source_x = 10;
source_y = JE / 2;

t0 = 100.0;
spread = 30.0;

% Υπολογισμός της σταθεράς M για τις συνθήκες Mur
M = (c0 * dt - ddx) / (c0 * dt + ddx);

for n = 3:nsteps
    % Αποθήκευση της προηγούμενης τιμής του πεδίου
    ez_old = ez1;
    ez1 = ez;

    % Update Dz field
    for i = 2:IE-1
        for j = 2:JE-1
            dz(i, j) = dz(i, j) + 0.5 * (hy(i, j) - hy(i-1, j) - hx(i, j) + hx(i, j-1));
        end
    end

    % Source pulse
    x = sin(2 * pi * freq * dt * (n - t0));
    y = exp(-0.5 * ((n - t0)^2 / spread^2));
    pulse = x .* y;
    dz(source_x, source_y) = pulse;

    % Update Ez field
    for i = 2:IE-1
        for j = 2:JE-1
            ez(i, j) = ga(i, j) * dz(i, j);
        end
    end

    % Εφαρμογή των απορροφητικών συνθηκών Mur
    % Αριστερό όριο (x = 1)
    for j = 2:JE-1
        ez(1,j) = ez1(2,j) + M * (ez(2,j) - ez1(1,j));
    end

    % Δεξί όριο (x = IE)
    for j = 2:JE-1
        ez(IE,j) = ez1(IE-1,j) + M * (ez(IE-1,j) - ez1(IE,j));
    end

    % Κάτω όριο (y = 1)
    for i = 2:IE-1
        ez(i,1) = ez1(i,2) + M * (ez(i,2) - ez1(i,1));
    end

    % Πάνω όριο (y = JE)
    for i = 2:IE-1
        ez(i,JE) = ez1(i,JE-1) + M * (ez(i,JE-1) - ez1(i,JE));
    end

    % Update Hx and Hy fields
    for i = 1:IE-1
        for j = 1:JE-1
            hx(i, j) = hx(i, j) + 0.5 * (ez(i, j) - ez(i, j+1));
            hy(i, j) = hy(i, j) + 0.5 * (ez(i+1, j) - ez(i, j));
        end
    end

    % Εμφάνιση του πεδίου κάθε 10 βήματα
    if mod(n, 10) == 0
        % Plot 2D image of Ez field
        subplot(2, 1, 1)
        imagesc(1:IE, 1:JE, ez')
        shading interp
        colorbar
        colormap(jet)
        axis([0 1000 0 1000])
        xlabel('x');
        ylabel('y');
        title(sprintf('Dielectric medium with circles (n = 2), t = %d', n))
        for k = 0:num_circles-1
            center_x = start_x + k * circle_spacing;
            rectangle('Position', [center_x-radius, center_y-radius , 2*radius, 2*radius], ...
                'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 1)
        end

        % Plot 3D surface of Ez field
        subplot(2, 1, 2)
        surf(1:IE, 1:JE, ez')
        shading interp
        colorbar
        colormap(jet)
        axis([0 1000 0 1000 -0.15 0.15]);
        xlabel('x');
        ylabel('y');
        zlabel('Ez');
        title(sprintf('3D View of Ez field at t = %d', n))
        for k = 0:num_circles-1
            center_x = start_x + k * circle_spacing;
            rectangle('Position', [center_x-radius, center_y-radius , 2*radius, 2*radius], ...
                'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 1)
        end

        pause(0.1);
    end
end
