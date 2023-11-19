function Det = Diff(A)
   [m n] = size(A);
   Det = zeros(m*2,n);
   for i = 1:m-1
       for j = 1:n-1
           Det(i,j) = A(i+1,j) - A(i,j);
           Det(i+m,j) = A(i,j+1) - A(i,j);
       end
   end
   for j = 1:n
       Det(m,j) = A(1,j) - A(m,j);
   end
   for i = 1:m
       Det(i+m,n) = A(i,1) - A(i,n);
   end
end
