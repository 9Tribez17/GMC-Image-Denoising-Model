function y = normgrad_a(X)

[i,j,L] = size(X);
normM = zeros(i,j);
for l = 1 : L
    normM = normM+(X(:,:,l)).^2;
end

%y = zeros(i,j/2);
%delta = 0.0001; 
%delta = 0.0001;
% X1 = X(1:i,1:j/2);
% X2 = X(1:i,(j/2+1):j); 

y = sqrt(normM);