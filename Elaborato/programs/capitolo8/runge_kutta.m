function [yT,tab] = runge_kutta(odefun,t0,T,h,y0,A,b,c,order,toll)
%RUNGE_KUTTA Metodo Runge-Kutta per il problema y'(t)=odefun(t,y(t)), y(t0)=y0.
% L'incognita y(t) è un vettore colonna di dimensione m.
% L'intervallo di tempo [t0,T] viene discretizzato con passo h uniforme.
% Il metodo è caratterizzato dai nodi c, pesi b e matrice di Butcher A.
% L'argomento order serve a calcolare in modo automatico una precisione
% adeguata con cui risolvere i sistemi di equazioni non lineari (utile
% solo nel caso di metodi semi-impliciti o impliciti). L'argomento
% facoltativo toll permette di stabilire manualmente tale precisione.
% L'argomento in uscita yT è l'approssimazione numerica di y(T).
% L'argomento in uscita tab è facoltativo, è di tipo "table" e contiene
% tutte le iterazioni effettuate dal metodo Runge-Kutta.

    % Controlla e normalizza i dati in ingresso
    assert(t0 <= T);
    assert(h > 0);
    N = round((T-t0)/h);
    h = (T-t0)/N;
    m = size(y0,1);
    assert(size(y0,2)==1);
    s = length(c);
    assert(length(b) == s);
    assert(all(size(A) == [s,s]));
    assert(order >= 1);
    if nargin < 10
        toll = 0.1*(h^(order+1));  % un decimo dell'err. locale di troncamento,
        toll = max(toll,eps);      % ma non meno dell'epsilon di macchina
    end
    
    % Tipo di metodo
    if any(any(triu(A,1)~=0))
        algorithm = 2;         % metodo implicito
    elseif any(diag(A)~=0)
        algorithm = 1;         % metodo semi-implicito
    else
        algorithm = 0;         % metodo esplicito
    end
    
    % Riserva memoria per y_n, Y_i, odefun(t,Y_i)
    y = y0;
    Y = zeros(m,s);
    fY = zeros(m,s);
    
    % Inizializza la tabella (facoltativo)
    save_table = (nargout > 1);
    if save_table
        varNames = {'t'};
        varTypes = {'double'};
        for i = 1:m
            varNames{i+1} = ['y',num2str(i)];
            varTypes{i+1} = 'double';
        end
        tab = table('Size',[N+1,1+m],...
            'VariableTypes',varTypes,'VariableNames',varNames);
        tab(1,1) = {t0};
        tab(1,2:end) = array2table(y0');
    end
    
    % Loop principale: a ogni iterazione viene calcolato y_n \approx y(t)
    n = 1;
    t = t0+h;
    while n <= N
        if algorithm == 0
            % metodo esplicito: sostituzioni in avanti
            for i = 1:s
                Y(:,i) = y;
                for j = 1:(i-1)
                    Y(:,i) = Y(:,i) + h*A(i,j)*fY(:,j);
                end
                fY(:,i) = odefun(t+c(i)*h,Y(:,i));
            end
        elseif algorithm == 1
            % metodo semi-implicito: s sistemi con m incognite
            fy = odefun(t,y);
            for i = 1:s
                g = y;
                for j = 1:(i-1)
                    g = g + h*A(i,j)*fY(:,j);
                end
                % previsione basata sul metodo di Eulero esplicito
                Yi_predicted = y + h*c(i)*fy;
                [Y(:,i),errcode] = RK_semiimplicit_solver(...
                    odefun, t+c(i)*h, Yi_predicted, h, A(i,i), g, toll);
                if errcode == 1
                    error('Il sistema nonlineare non è stato risolto.');
                end
                fY(:,i) = odefun(t+c(i)*h,Y(:,i));
            end
        else
            % metodo implicito: un unico sistema con s*m incognite
            fy = odefun(t,y);
            Y_predicted = zeros(m,s);
            for i = 1:s
                % previsione basata sul metodo di Eulero esplicito
                Y_predicted(:,i) = y + h*c(i)*fy;
            end
            [Y(:),errcode] = RK_implicit_solver(...
                odefun, t, y, Y_predicted(:), h, A, c, toll);
            if errcode == 1
                error('Il sistema nonlineare non è stato risolto.');
            end
            for i = 1:s
                fY(:,i) = odefun(t+c(i)*h,Y(:,i));
            end
        end
        
        % aggiornamento di y
        for i = 1:s
            y = y + h*b(i)*fY(:,i);
        end
        if save_table
            tab(n+1,1) = {t};
            tab(n+1,2:end) = array2table(y');
        end
        n = n + 1;
        t = t + h;
    end
    yT = y;
end
