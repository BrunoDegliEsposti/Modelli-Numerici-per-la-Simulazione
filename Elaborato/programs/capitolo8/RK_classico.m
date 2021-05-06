function [A,b,c,order] = RK_classico()
%RK_CLASSICO Tableau di Butcher del metodo Runge-Kutta classico
    A = [0,  0,  0,  0;
         0.5,0,  0,  0;
         0,  0.5,0,  0;
         0,  0,  1,  0];
    b = [1/6,1/3,1/3,1/6];
    c = [0,0.5,0.5,1];
    order = 4;
end
