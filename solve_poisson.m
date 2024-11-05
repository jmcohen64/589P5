function [u, omega, rho, A] = solve_poisson(f, varargin)
% SOLVE_POISSON - solve Poisson PDE on the unit square
%   U = SOLVE_POISSON(F, VARARGIN) solves the Poisson equation
% on the grid of size N. The argument N is an integer >=1.
% The argument F is an N-by-N matrix containing the values
% of the source density function over the uniform grid on [0,1]^2. The
% returned 2D-array U is of the same size as F. The boundary
% values of U must be 0, corresponding to the Dirichlet boundary
% conditions. The function accepts the following keyword parameters.
%
%   'Method'    - the method. One of 'Jacobi', 'Gauss-Seidel' and 'SOR'.
%                 Default: 'Jacobi'
%   'Omega'     - only used if the 'Method' is 'SOR'. The over-relaxation
%                 parameter. If 'Omega' is set to [] then you should optimize
%                 the value of 'Omega' to minimize the relevant spectral radius.
%                 Default: []
%   'Tolerance' - The tolerance in fixed-point iteration.
%                 Default: 1e-6     
%
% [U,OMEGA,RHO] = SOLVE_POISSON(N, F, ..., 'Omega', [], ...) should also
% return the optimized value of 'Omega' and spectral radius 'Rho' for
% that 'Omega'.
    p = inputParser;
    p.addRequired('f');
    p.addParameter('Method','Jacobi');    
    p.addParameter('Omega', []);
    p.addParameter('Tolerance', 1e-6);
    p.parse(f, varargin{:});
        
    [N,Ncol]=size(f);
    assert(N==Ncol);
    A = build_matrix(N);
    switch(p.Results.Method)
      case 'Jacobi',
        [u, omega, rho] = jacobi(A, f(:), p.Results.Tolerance);
      case 'Gauss-Seidel',
        [u, omega, rho] = gauss_seidel(A, f(:), p.Results.Tolerance);
      case 'SOR',
        [u, omega, rho] = sor(A, f(:), ...
                              p.Results.Omega,...
                              p.Results.Tolerance);
      otherwise,
        error(['Invalid method: ',p.Results.Method])
    end
    u = reshape(u,N,[]);
    assert(all(size(u) == size(f)));
end


function A = build_matrix(N)
    M = N.^2;
    L = sparse(M,M);
    % Put code that returns A here
    %calculate remainder when divided by 4
    r4 = mod(M,4);
    %divisor
    m = (M-r4)/4;
    if r4 == 3
        J_r = [4, -1, 0; -1, 4,-1; 0, -1, 4];
        r0 = [0,0,0,0; 0,0,0,0; 0,0,0,0];
        r1 = [1,0,0,0; 0,1,0,0; 0,0,1,0];
    elseif r4 == 2
        J_r = [4, -1; -1, 4];
        r0 = [0,0,0,0; 0,0,0,0];
        r1 = [1,0,0,0; 0,1,0,0];
    else
        J_r = [4];
        r0 = [0,0,0,0];
        r1 = [1,0,0,0];
    end
    c0 = r0.';
    c1 = r1.';
    J_4 = [4, -1, 0, 0; -1, 4,-1, 0; 0, -1, 4,-1; 0, 0, -1, 4];
    I_4 = eye(4);
    dim = m+r4;
    disp(dim);
    if m==0
        A = J_r;
    else
        disp('looping');
        A = cellmat(dim,dim,1,1);
        disp(A)
        for i = 1:dim
            for j = 1:dim
                disp(A)
                disp(j)
                if (j*4 <= M) && (i*4 <= M)
                    if (i == j)
                        A{i,j} = J_4;
                    elseif j==(i+1) || j==(i-1)
                        A{i,j} = I_4;
                    else
                        A{i,j} = sparse(4,4);
                    end
                elseif (i == j)
                    A{i,j} = J_r;
                elseif (i == dim)
                    if (j == (dim-1))
                       A{i,j} = r1;
                    else
                        A{i,j} = r0;
                    end
                elseif (j==dim)
                    if (i == (dim-1))
                       A{i,j} = c1;
                    else
                        A{i,j} = c0;
                    end
                end
            end
        end


    end
    A = sparce(cell2mat(A));
    
end

function [u, omega, rho] = jacobi(A, f, tol)
    D = diag(diag(A));
    L = tril(A,-1);
    U = triu(A,1);    
    Dinv = inv(D);
    B = Dinv*(U+L);
    f = Dinv*f;
    rho = spectral_radius(B);
    assert(rho < 1, 'Jacobi: Spectral radius is not < 1');
    u = sparse(size(D,1),1);
    while true
        u_new = f - B * u;
        if norm(u_new - u) < tol
            break;
        end
        u = u_new;
    end
    omega = NaN;
end

function [u, omega, rho] = gauss_seidel(A, f, tol)
    % Put your Gauss-Seidel code here
    % ... or make it run with mockup code
    D = diag(diag(A));
    u = sparse(size(D,1),1);
    omega = 1;
    rho = NaN;
end

function [u, omega, rho] = sor(A, f, omega, tol)
    % Put your SOR code here.
    % If omega is empty, estimate and return optimized omega
    % ... or make it run with mockup code
    D = diag(diag(A));
    u = sparse(size(D,1),1);
    omega = 1;
    rho = NaN;
end

function rho = spectral_radius(A)
    rho = abs( eigs(A,1,'largestabs','MaxIterations',512) );
end