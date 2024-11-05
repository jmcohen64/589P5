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