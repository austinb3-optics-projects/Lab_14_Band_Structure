clear; close all;


a = 1;                  % periodic potential period
nb = 10;                % number of potential barriers
b = 1/6;                % size scaling of potential barriers
V0 = 100;               % amplitude of potential wells

L = nb*a;                 % total length of propagation

test = 1;

%% Section A: Matrix Method, part b

% plot the periodic potential

nx = 10000;             % number of plotting points

x = linspace(0,L,nx);   % position vector

rect = @(y) double(abs(y) <= 1/2);     % anonymous function to define the rect shape


V = generatePotential(V0,x,nb,b,a,rect);

figure;
plot(x,V,'LineWidth',1.5);
title('Periodic potential');
ylabel('$V(x)$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
set(gca,'FontSize',15);

%% Section A: Matrix Method, part c

% This section is concerned with calculating the Hamiltonian for a given
% potential. This section references Function (2)
N = 100;        % matrix size

H = generateHamiltonian(V,N,L,x);     % Hamiltonian
fprintf('H11 = %.2f, H33 = %.2f, H55 = %.2f\n', H(1,1), H(3,3), H(5,5));

%% Section B: Band Structure, part a: extended zone scheme

% calculate the eigenvalues and eigenvectors of the hamiltonians
kmax = 29;
% free particle eigenvalues
En0 = zeros(1,kmax);
for n = 1:kmax
    En0(n) = pi^2*n^2./nb.^2;
end

[Psi,En] = Eigen(H,kmax);

% this is the normalized k-vector kn/pi
kn = (1:1:kmax)./10;

figure;
plot(kn,En0./pi.^2,'b--','LineWidth',1.2); hold on;
plot(kn,En./pi.^2,'ro'); hold off;
title('Free particle energy and periodic potential well energy');
ylabel('$E_{n}/\pi^{2} , E_{n}^{(0)}/\pi^{2}$','Interpreter','latex');
xlabel('$k_{n}/\pi$','Interpreter','latex');
legend('E_{n}^{(0)}','E_{n}','Location','northwest');
set(gca,'FontSize',15);

%% Section B: Band Structure, part b: reduced zone scheme

kn1 = (1:1:9);
kn2 = (20:-1:10);
kn3 = (21:1:29);


En1 = En(kn1)./pi.^2;
En2 = En(kn2)./pi.^2;
En3 = En(kn3)./pi.^2;

kn1 = (1:1:9)/10;
kn2 = (0:1:10)/10;
kn3 = (1:1:9)/10;

figure;
plot(kn1,En1,'b','LineWidth',1.5); hold on;
plot(kn2,En2,'r','LineWidth',1.5);
plot(kn3,En3,'g','LineWidth',1.5); hold off;
title('Reduced Brillouin zone');
ylabel('$E_{n}/\pi^{2}$','Interpreter','latex');
xlabel('$k_{n}/\pi$','Interpreter','latex');
legend('j=1','j=2','j=3','Location','best');
xlim([0,1]);
set(gca,'FontSize',15);

%% Section C: Bloch standing waves, part a and b

% this section uses Function (4)
% n1 = 16;

[Cn,En] = Eigen(H,N);

Phi = boxModes(x,L,N);
Psi = (Cn.')* (Phi.');
Phi = 0.6.*sqrt(L/2).*boxModes(x,L,16);

figure;     % Psi1
plot(x,-Psi(1,:),'LineWidth',1.5); hold on;
plot(x,Phi(:,1),'--','LineWidth',1.5); hold off;
title('Eigenfunctions and periodic expansion functions');
ylabel('$\psi_{1} , \varphi_{1}$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
legend('$\psi_{1}(x)$','$\varphi_{1}(x)$','Interpreter','latex');
set(gca,'FontSize',15);

figure;     % Psi9
plot(x,Psi(9,:),'LineWidth',1.5); hold on;
plot(x,-Phi(:,9),'--','LineWidth',1.5); hold off;
title('Eigenfunctions and periodic expansion functions');
ylabel('$\psi_{9} , \varphi_{9}$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
legend('$\psi_{9}(x)$','$\varphi_{9}(x)$','Interpreter','latex');
set(gca,'FontSize',15);

figure;     % Psi10
plot(x,-Psi(10,:),'LineWidth',1.5); hold on;
plot(x,Phi(:,10),'--','LineWidth',1.5); hold off;
title('Eigenfunctions and periodic expansion functions');
ylabel('$\psi_{10} , \varphi_{10}$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
legend('$\psi_{10}(x)$','$\varphi_{10}(x)$','Interpreter','latex');
set(gca,'FontSize',15);

figure;     % Psi16
plot(x,-Psi(16,:),'LineWidth',1.5); hold on;
plot(x,-Phi(:,16),'--','LineWidth',1.5); hold off;
title('Eigenfunctions and periodic expansion functions');
ylabel('$\psi_{16} , \varphi_{16}$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
legend('$\psi_{16}(x)$','$\varphi_{16}(x)$','Interpreter','latex');
set(gca,'FontSize',15);

%% Section D: Tight-Binding limit, part a

% plotting the wavefunctions probability
figure;     % Psi1
plot(x,abs(Psi(1,:)).^2,'LineWidth',1.5);
title('Probability in the tight binding limit');
ylabel('$\psi_{1}$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
set(gca,'FontSize',15);

%% Section D: Tight-Binding limit, part b

% this part looks at deconstructing the ground state probability into
% localized peaks that are note independent.

Psi0 = sqrt(2).*exp(-6.25.*(x - 5).^2);     % initial conditions
t = linspace(0,2,500);


% calculate expansion coefficients
d = zeros(N,1);
for p = 1:N
    d(p) = trapz(x,conj(Psi(p,:)).*Psi0);
end


PSI_temp = zeros(length(t),length(x));
PSI = zeros(length(t),length(x));

for j = 1:length(d)
    for k = 1:length(t)
        PSI_temp(k,:) = d(j)*Psi(j,:).*exp(-1i*t(k)*En(j));
    end
    if j == 1 
        PSI = PSI_temp;
    else
        % sum up the weighted eigenvectors
        PSI = PSI + PSI_temp;
    end
end
PSI = PSI.';
figure;
imagesc([0,2],[0,L],abs(PSI).^2);
xlabel('$t$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
title('$|\Psi(x,t)|^2$','Interpreter','latex');
set(gca,'FontSize',15);
colormap turbo;
colorbar;


%% Functions

% Function (1)
function [V] = generatePotential(V0,x,nb,b,a,rect)

    V = zeros(1,length(x)); % initialize potential function
    for r = 1:nb
        xr = -a/2 + r;      % calculate center of rect function
        V = V + V0.*rect((x - xr)./b);
    end

end

% Function (2)
function [H] = generateHamiltonian(V,N,L,x)
    % numeric integration is going to happen with trapz function
    H = zeros(N,N);
    for n = 1:N     % rows
        for m = 1:N         % columns
            if n == m
                H(n,m) = ((pi*n)^2/L^2) + 2./L.*trapz(x,sin(n.*pi.*x./L).*V.*sin(m.*pi.*x./L));
            else
                H(n,m) = 2./L.*trapz(x,sin(n.*pi.*x./L).*V.*sin(m.*pi.*x./L));
            end
        end
    end 
end

% Function (3)
function [Psi,En] = Eigen(H,N)
    En = eigs(H,N,'smallestabs');
    [Psi,~] = eigs(H,N,'smallestabs');
    % I use the eigs function twice, because the formatting of the second
    % expression for the eigenvalues sucks.
end

% Function (4)
function [phi] = boxModes(x,L,N)
    
    phi = zeros(length(x),N);
    for n = 1:N
        kn = n*pi/L;
        phi(:,n) = sqrt(2./L).*sin(kn.*x);
    end
end