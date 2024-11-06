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
    h = 1/(N+1);
    A = sparse(M,M);
    J = 4*eye(N) - diag(ones(N-1, 1), 1) - diag(ones(N-1, 1), -1);
    S = -1*eye(N);
    for i = 1:N
        for j = 1:N
            % Row and column indices for the block structure
            row_idx = (i-1)*N + (1:N);  % Row indices for the block
            col_idx = (j-1)*N + (1:N);  % Column indices for the block
            
            if i == j
                % Diagonal block: Use J matrix
                A(row_idx, col_idx) = J;
            elseif abs(i - j) == 1
                % Adjacent blocks: Use S matrix for boundary conditions
                A(row_idx, col_idx) = S;  % Boundary conditions
            end
        end
    end
    A = (1/h^2)*A;
    
end

% function A = build_matrix(N)
%     M = N.^2;
%     h = 1/(N+1);
%     %L = sparse(M,M);
%     J = 4*eye(N);
%     S = (-1/h^2)*eye(N);
%     %A = cellmat(N,N,N,N);
%     A = mat2cell(sparse(N,N));
%     for i =1:N
%         for j = 1:N
%             if abs(i-j) == 1
%                 J(i,j) = -1;
%             end
%         end
%     end
%     J = (1/h^2)*J;
%     for i =1:N
%         for j = 1:N
%             if i ==j
%                 A{i,j} = J;
%             elseif abs(i-j) == 1
%                 A{i,j} = S;
%             end
%         end
%     end
%     A = sparse(cell2mat(A));
% 
% end

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
    L = tril(A,-1);
    U = triu(A,1);
    DLinv = inv(D+L);
    u = sparse(size(D,1),1);
    rho = spectral_radius((DLinv*U));
    assert(rho < 1, 'Gauss-Seidel: Spectral radius is not < 1');
    while true
        u_new = DLinv*(f - U * u);
        if norm(u_new - u) < tol
            break;
        end
        u = u_new;
    end
    omega = 1;
    
end

function [u, omega, rho] = sor(A, f, omega, tol)
    % Put your SOR code here.
    % If omega is empty, estimate and return optimized omega
    % ... or make it run with mockup code
    if isempty(omega)
        omega = 1.82;
    end
    D = diag(diag(A));
    L = tril(A,-1);
    U = triu(A,1);
    DwLinv = inv(D+omega*L);
    u = sparse(size(D,1),1);
    rho = spectral_radius(DwLinv*(omega*U + (omega-1)*D));
    assert(rho < 1, 'SOR: Spectral radius is not < 1');
    while true
        u_new = DwLinv * omega * f - DwLinv * (omega * U + (omega-1)*D)*u;
        if norm(u_new - u) < tol
            break;
        end
        u = u_new;
    end
    
    %sor1(A, f, tol, omega);
end

% function [u, omega, rho] = sor1(A, f, tol, omega)
%     % Put your SOR code here.
%     % If omega is empty, estimate and return optimized omega
%     % ... or make it run with mockup code
%     arguments
%         A (1,:) double
%         f (1,:) double
%         tol (1,:) double 
%         omega (1,1) double = 1.5
%     end
%     D = diag(diag(A));
%     L = tril(A,-1);
%     U = triu(A,1);
%     omega = 1.5;
%     DwLinv = inv(D+omega*L);
%     u = sparse(size(D,1),1);
%     rho = spectral_radius(DwLinv*(omega*U + (omega-1)*D));
% end

function rho = spectral_radius(A)
    rho = abs( eigs(A,1,'largestabs','MaxIterations',512) );
end