function [H] = kepler_hamiltonian(y,G,m1,m2)
%KEPLER_HAMILTONIAN Calcola la funzione hamiltoniana del problema di Keplero.
% L'hamiltoniana è associata al moto in R^2 di un corpo con massa m2
% attorno a un corpo inamovibile posizionato nell'origine con massa m1>>m2
% per effetto della forza di gravità.
    if nargin < 2
        G = 1;
    end
    if nargin < 3
        m1 = 1;
    end
    if nargin < 4
        m2 = 1;
    end
    q = y(1:2,:);
    p = y(3:4,:);
    H = 0.5*dot(p,p,1) - G*m1*m2./vecnorm(q,2,1);
end
