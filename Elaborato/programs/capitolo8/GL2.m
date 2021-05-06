function [A,b,c,order] = GL2()
%GL2 Tableau di Butcher del metodo di collocazione di Gauss-Legendre a due stadi
    r = -sqrt(3)/6 + 1/2;
    A = [1/4,1/4-sqrt(3)/6;
         1/4+sqrt(3)/6,1/4];
    b = [1/2,1/2];
    c = [r, 1-r];
    order = 4;
end
