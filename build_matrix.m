function A = build_matrix(N)
    M = N.^2;
    h = 1/(N+1);
    %L = sparse(M,M);
    %initialzie the two nxn matrices
    J = 4*eye(N);
    S = (-1/h^2)*eye(N);
    %initialize A as nxn cell array (matrix) with elements that are nxn matrices
    A = cellmat(N,N,N,N);
    % modify J to add -1 to the elements 1 off of the diagonal
    for i =1:N
        for j = 1:N
            if abs(i-j) == 1
                J(i,j) = -1;
            end
        end
    end
    %scale J
    J = (1/h^2)*J;
    %modify our matrix of matrices
    %add J to the diagonal and S to 1 off the diagonal; rest are 0s by
    %construction of A
    for i =1:N
        for j = 1:N
            if i ==j
                A{i,j} = J;
            elseif abs(i-j) == 1
                A{i,j} = S;
            end
        end
    end
    %convert to matrix
    A = full(cell2mat(A));
    
end


% function A = build_matrix(N)
%     M = N.^2;
%     h = 1/(N+1);
%     L = sparse(M,M);
%     % Put code that returns A here
%     %calculate remainder when divided by 4
%     r4 = mod(M,4);
%     %divisor
%     m = (M-r4)/4;
%     if r4 == 3
%         J_r = (1/h^2)*[4, -1, 0; -1, 4,-1; 0, -1, 4];
%         r0 = [0,0,0,0; 0,0,0,0; 0,0,0,0];
%         r1 = (-1/h^2)*[1,0,0,0; 0,1,0,0; 0,0,1,0];
%     elseif r4 == 2
%         J_r = (1/h^2)*[4, -1; -1, 4];
%         r0 = [0,0,0,0; 0,0,0,0];
%         r1 = (-1/h^2)*[1,0,0,0; 0,1,0,0];
%     else
%         J_r = (1/h^2)*[4];
%         r0 = [0,0,0,0];
%         r1 = (-1/h^2)*[1,0,0,0];
%     end
%     c0 = r0.';
%     c1 = r1.';
%     J_4 = (1/h^2)*[4, -1, 0, 0; -1, 4,-1, 0; 0, -1, 4,-1; 0, 0, -1, 4];
%     I_4 = (-1/h^2)*eye(4);
%     dim = m+r4;
%     %disp(dim);
%     if m==0
%         A = J_r;
%     else
%         %disp('looping');
%         A = cellmat(dim,dim,1,1);
%         %disp(A)
%         for i = 1:dim
%             for j = 1:dim
%                 %disp(A)
%                 %disp(j)
%                 if (j*4 <= M) && (i*4 <= M)
%                     if (i == j)
%                         A{i,j} = J_4;
%                     elseif j==(i+1) || j==(i-1)
%                         A{i,j} = I_4;
%                     else
%                         A{i,j} = sparse(4,4);
%                     end
%                 elseif (i == j)
%                     A{i,j} = J_r;
%                 elseif (i == dim)
%                     if (j == (dim-1))
%                        A{i,j} = r1;
%                     else
%                         A{i,j} = r0;
%                     end
%                 elseif (j==dim)
%                     if (i == (dim-1))
%                        A{i,j} = c1;
%                     else
%                         A{i,j} = c0;
%                     end
%                 end
%             end
%         end
% 
% 
%     end
%     A = full(cell2mat(A));
%     %disp(class(A));
% 
% end