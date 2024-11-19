% Main script: KronigPenney_Main.m

% Clear workspace and figures
clear all;
close all;

% Define global parameters
nb = 10;    % number of barriers
b = 1/6;    % barrier width
V0 = 100;   % barrier height
L = nb;     % total length (a=1 in scaled units)
nx = 10000; % number of spatial points
N = 100;    % number of basis modes

% Part I.A.(b) - Plot potential
x = linspace(0, nb, nx);
V = generate_potential(x, nb, b, V0);

figure(1);
plot(x, V);
title('Kronig-Penney Potential');
xlabel('x'); ylabel('V(x)');
grid on;

% Part I.A.(c) - Generate and verify Hamiltonian matrix
H = generate_hamiltonian(N, x, V);
fprintf('H11 = %.2f, H33 = %.2f, H55 = %.2f\n', H(1,1), H(3,3), H(5,5));

% Part I.B - Calculate and plot band structure
[E, vectors] = eig(H);
E = diag(E);
[E, idx] = sort(E);
vectors = vectors(:, idx);

% Calculate free particle energies
E0 = (pi^2 * (1:N).^2)/nb^2;

% Plot band structure (extended zone)
k = (1:29)/10;
figure(2);
plot(k, E(1:29)/pi^2, 'bo', k, E0(1:29)/pi^2, '--');
title('Band Structure (Extended Zone)');
xlabel('k/\pi'); ylabel('E/\pi^2');
grid on;
legend('With potential', 'Free particle');

% Calculate wavefunctions for selected states
n_plot = [1, 9, 10, 16];
for i = 1:length(n_plot)
    n = n_plot(i);
    psi = calculate_wavefunction(vectors(:,n), x, N, L);
    
    figure(2+i);
    plot(x, psi, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(x, 0.6*sin(n*pi*x/L), 'r:', 'LineWidth', 1);
    title(['Eigenfunction \psi_' num2str(n) '(x)']);
    xlabel('x'); ylabel('\psi(x)');
    legend('Bloch wave', 'Box mode');
    grid on;
end

% Time evolution simulation
dt = 0.004;
t = 0:dt:2;
x0 = 5;
sigma = 1/sqrt(12.5);
Psi0 = sqrt(2)*exp(-(x-x0).^2/(2*sigma^2));

% Calculate expansion coefficients
d = zeros(N,1);
for n = 1:N
    psi = calculate_wavefunction(vectors(:,n), x, N, L);
    d(n) = trapz(x, conj(psi) .* Psi0);
end

% Time evolution
[X, T] = meshgrid(x, t);
Psi = zeros(length(t), length(x));

for i = 1:length(t)
    Psi_t = zeros(1, length(x));
    for n = 1:N
        psi = calculate_wavefunction(vectors(:,n), x, N, L);
        Psi_t = Psi_t + d(n) * psi * exp(-1i*E(n)*t(i));
    end
    Psi(i,:) = abs(Psi_t).^2;
end

figure(7);
imagesc(x, t, Psi);
colormap(jet);
xlabel('x'); ylabel('t');
title('Probability density evolution');
colorbar;

% Function: generate_potential.m
function V = generate_potential(x, nb, b, V0)
    V = zeros(size(x));
    for r = 1:nb
        xr = -0.5 + r;
        V = V + V0 * (abs(x - xr) <= b/2);
    end
end

% Function: generate_hamiltonian.m
function H = generate_hamiltonian(N, x, V)
    L = max(x);
    dx = x(2) - x(1);
    H = zeros(N, N);
    
    for n = 1:N
        for m = 1:N
            % Kinetic energy term
            if n == m
                H(n,m) = (pi^2 * n^2)/(2*L^2);
            end
            
            % Potential energy term
            phi_n = sqrt(2/L) * sin(n*pi*x/L);
            phi_m = sqrt(2/L) * sin(m*pi*x/L);
            H(n,m) = H(n,m) + trapz(x, phi_n .* V .* phi_m);
        end
    end
end

% Function: calculate_wavefunction.m
function psi = calculate_wavefunction(coefficients, x, N, L)
    psi = zeros(size(x));
    for m = 1:N
        phi_m = sqrt(2/L) * sin(m*pi*x/L);
        psi = psi + coefficients(m) * phi_m;
    end
end