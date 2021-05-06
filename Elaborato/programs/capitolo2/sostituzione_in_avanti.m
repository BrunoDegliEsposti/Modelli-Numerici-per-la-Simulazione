y(n0:n0+k-1) = initial_conditions;
for n = n0:nmax
	y(n+k) = g(n);
	for i = 1:k
		y(n+k) = y(n+k) - p(i,n)*y(n+k-i);
	end
end
