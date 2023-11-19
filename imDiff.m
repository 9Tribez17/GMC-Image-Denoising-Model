function Det = imDiff(A)
    [m n] = size(A);
    A1 = zeros(m/2,n);
    A2 = zeros(m/2,n);
    Det1 = zeros(m/2,n);
    Det2 = zeros(m/2,n);    
    for i = 1:m/2
        for j = 1:n
            A1(i,j) = A(i,j);
        end
    end
    for i = 1:m/2
        for j = 1:n
            A2(i,j) = A(i+m/2,j);
        end
    end
    for i = 1:m/2
        for j = 1:n
            if i == 1
                Det1(i,j) = A1(m/2,j) - A1(i,j);
            end
            if i>1 && i<=(m/2)
                Det1(i,j) = A1(i-1,j) - A1(i,j);
            end
        end
    end
    for i = 1:m/2
        for j =1:n
            if j == 1
                Det2(i,j) = A2(i,n) -A2(i,j);
            end
            if j>1 && j<=n
                Det2(i,j) = A2(i,j-1) - A2(i,j);
            end
        end
    end
    Det = Det1 + Det2;
end