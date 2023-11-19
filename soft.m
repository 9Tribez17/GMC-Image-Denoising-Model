function D = soft(A,lambda)
    [m n] = size(A);
    A1 = zeros(m/2,n);
    A2 = zeros(m/2,n);
  
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
    
    v  = sqrt(A1.^2 + A2.^2);
    d1 = max(v - lambda,0) .* (A1./ v);
    d2 = max(v - lambda,0) .* (A2./ v);
    D  = [d1 ; d2];

