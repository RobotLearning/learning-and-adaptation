%% Deal function values to an array

fun = @(x) x^2;
A = rand(5,5);

A2 = arrayfun(fun, A); 

%% Get variable output

fun = @(x1,x2) x1*x2;

a = 1:7;
b = 3:9;
[A,B] = ndgrid(a,b);
fun(A,B)