%----------------------------------------------------------------
% File:     script.m
%----------------------------------------------------------------
%
% Author:   Marek Rychlik (rychlik@arizona.edu)
% Date:     Mon Oct 21 06:59:33 2024
% Copying:  (C) Marek Rychlik, 2020. All rights reserved.
% 
%----------------------------------------------------------------
% Driver for iterative methods
methods = {'SOR'};
%{'Jacobi','Gauss-Seidel',
N = 29;
omega_guess = [];
tol = 1e-13;
xpos=3./4; ypos=1./3;
%xpos=1./2; ypos=1./2;
imgdir=fullfile(pwd,'images');
mkdir(imgdir);

for j = 1:numel(methods)
    method = methods{j};
    disp('----------------');
    disp(method);
    disp('----------------');
    f = zeros(N,N);
    h = 1./(N+1);
    f(round(xpos.*N),round(ypos.*N)) = 1./h.^2;
    [u,omega,rho,A] = solve_poisson(f, ...
                                    'Method', method, ...
                                    'Omega', omega_guess,...
                                    'Tolerance', tol);
    u = padarray(u,[1,1],0,'both');
    imagesc(u);
    title(sprintf('%s: N=%3d, omega=%-.4g, rho=%-.4g',...
                  method, N, omega, rho));
    drawnow;
    filename=[method,num2str(N),'.png'];
    saveas(gcf, fullfile(imgdir,filename));
    pause(2);
end

% % Quiver plot of the gradient
% [X, Y] = meshgrid(0:(N+1), 0:(N+1));
% [FX, FY] = gradient(u);
% clf;
% quiver(X,Y,FX,FY);
% drawnow;
% pause(2);
% filename=['gradient',num2str(N),'.png'];
% saveas(gcf,fullfile(imgdir,filename));
% 
% % Cuthill-McKee reordering
% d = symrcm(A);
% spy(A(d,d));
% filename=['cuthill_mckee',num2str(N),'.png'];
% saveas(gcf,fullfile(imgdir,filename));
% 
% % Approximate minimum-degree reordering
% r = amd(A);
% spy(A(r,r));
% filename=['amd',num2str(N),'.png'];
% saveas(gcf,fullfile(imgdir,filename));