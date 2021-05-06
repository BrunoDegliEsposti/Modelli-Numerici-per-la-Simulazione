function [A,b,c,order] = SSPIRK33()
%SSPIRK33 Tableau di Butcher del metodo Runge-Kutta SSPIRK33
    r = -sqrt(2)/4 + 1/2;
    A = [r,0,0;
         sqrt(2)/4,r,0;
         sqrt(2)/4,sqrt(2)/4,r];
    b = [1/3,1/3,1/3];
    c = [r, 1/2, 1-r];
    order = 3;
end
