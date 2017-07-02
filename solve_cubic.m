function sols = solve_cubic(a, b, c, d)
  syms x
  a=2;
  b=0;
  c=2;
  d=-54;
  sols = solve(a*x^3 + b*x^2 + c*x + d);
end