function [x1,fx1,n,errcode,T] = secanti(...
    f,x0,x1,tolla,tollr,tollf,nmax,save_table)
% SECANTI Metodo delle secanti per f: R -> R
% Parametri in ingresso:
%   f: funzione di cui si cerca uno zero
%   x0, x1: approssimazioni iniziali distinte di uno zero di f
%   tolla: tolleranza assoluta per il criterio di arresto sulle iterate
%   tollr: tolleranza relativa per il criterio di arresto sulle iterate
%   tollf: tolleranza per il criterio di arresto sulla funzione
%   nmax: il metodo termina dopo aver calcolato x_nmax (evita loop infiniti)
%   save_table: flag booleano; salva in T lo storico delle iterazioni
% Parametri in uscita:
%   x1: approssimazione dello zero calcolata all'ultima iterazione
%   fx1: valore di f in x1 calcolato all'ultima iterazione
%   n: il parametro in uscita x1 corrisponde a x_n
%   errcode: intero che codifica l'esito del metodo:
%     0 -> successo, almeno un criterio di arresto soddisfatto
%     1 -> errore, raggiunto nmax
%     2 -> errore, secante parallela all'asse x
%   T: tabella con informazioni sulla convergenza del metodo

    % Innesco del metodo
    fx0 = f(x0);
    fx1 = f(x1);
    if abs(fx1) <= tollf
        n = 1;
        errcode = 0;
        return
    end
    
    % Inizializza la tabella (facoltativo)
    if save_table
        T = table(0, x0, fx0, NaN, NaN, NaN, ...
                  'VariableNames', {'n','xn','fxn','absdx', ...
                                    'dxratio','dxratiophi'});
        T(2,:) = {1, x1, fx1, abs(x1-x0), NaN, NaN};
    end
    phi = (1+sqrt(5))/2;
    dx_old = x1-x0;
    
    for n = 2:nmax
        % Calcola un passo del metodo delle secanti. Le variabili
        % x0 e x1 vengono riciclate, così da usare meno memoria.
        % La funzione f viene valutata solo una volta per iterazione.
		df = fx1 - fx0;
		if df == 0
			errcode = 2;
			return
		end
		dx = -fx1 * (x1-x0) / df;
        x0 = x1;
        fx0 = fx1;
        x1 = x1 + dx;
        fx1 = f(x1);
        
        % Salva informazioni sull'iterazione corrente (facoltativo)
        if save_table
            T(end+1,:) = {n, x1, fx1, abs(dx), ...
                abs(dx)/abs(dx_old), abs(dx)/abs(dx_old)^phi};
        end
        dx_old = dx;
        
        % Criterio di arresto misto sulle iterate e sulla funzione
        if abs(dx) <= tolla+tollr*abs(x1) || abs(fx1) <= tollf
			errcode = 0;
			return
        end
    end
    % Numero massimo di iterazioni raggiunto
	errcode = 1;
end
