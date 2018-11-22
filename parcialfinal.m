



disp('------------------------------------------------------------------');
disp('EJERCICIO 2:');
disp('------------------------------------------------------------------');

x = 0:0.5:6;
y = [2 5.4375 7.3516 7.5625 8.4453 9.1875 12];
t = 3.5;
plagr(x, y, t);






function suma = plagr(x, y, t)
    n = length(x);
    suma = 0;

    for i = 1:n
        prod = 1;
        for j = 1:n
            if i ~= j
                prod = prod .* (t - x(j))./(x(i) - x(j));
            end
        end
        suma = suma + y(i).*prod
    end

end


function dy = deriv(x, y)
    n = length(x);
    dy(1) = (y(2) - y(1))/(x(2) - x(1));
    
    for i = 2:n - 1
        dy(i) = (y(i + 1) - y(i - 1))/(x(i + 1) - x(i - 1))
    end
    
    dy(n)  =(y(n) - y(n - 1))/(x(n) - x(n - 1))
end

function sum = trapecio(f, a, b, n)
    if a == b
        sum = 0;
        return;
    end 
    h = (b-a)/n
    sum = f(a)
    x = a
    for i =1 : n-1
        x = x+h
        sum = sum + 2*f(x)
    end 
    sum = sum + f(b)
    sum = h*sum/2
    calculoError(sum, f, a, b);
end

function sum = simpson13(f, a, b, n)
    if a == b
        sum = 0;
        return;
    end

    h = (b-a)/n;
    x = a:h:b;
    fx = f(x);
    sum = 0;

    for i=1: n-1
        if mod(i,2) == 0
            sum = sum + 4*fx(i)
        else
            sum = sum + 2*fx(i)
        end
    end

    sum = sum + f(b)
    sum = sum*h/3
    calculoError(sum, f, a, b);
end


function sum = simpson38(f, a, b)
    if a == b
       sum = 0;
       return;
    end
    d = linspace(a, b, 4)
    fd = f(d)
    sum = (fd(1)+3*fd(2)+3*fd(3)+fd(4))*(b-a)/8
    calculoError(sum, f, a, b);
end

function calculoError(Iaprox, f, a, b)
    disp('------------------------------------------------------------------');
    Ireal = integral(f, a, b);
    errorVerdadero = abs(Ireal - Iaprox);
    errorRelativoPorcentual = 100*abs((Ireal - Iaprox)/Ireal);
    fprintf('Ireal = %f\n', Ireal);
    fprintf('Error Verdadero = %f\n', errorVerdadero);
    fprintf('Error Relativo Porcentual = %f\n', errorRelativoPorcentual);
    disp('------------------------------------------------------------------');
end

f = @(x) exp(-0.5*x).*(4-x)-2;
df = @(x) (exp(-x/2)*(x - 4))/2 - exp(-x/2);

x1 = 2.7;
tol = 1e-10;

while abs(f(x)) > tol
    f(x)
    % Newton-Raphson iteration
    x = x1 - f(x1)/df
    input('go')
end