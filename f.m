function m = f(x,t, epsilon, omega)
	m = a(t,epsilon, omega) .*x.*x + b(t,epsilon, omega) .* x;
end
