function L = pwcp2quality(x,P,N)
 L = 0;
for i=1:N
    for j=1:N
        L = L + (P(i,j)*log(x(i)/(x(i) + x(j))));
    end
end
L=-L;
end

