function [yT,tab] = LMF(odefun,t0,T,h,y0,alpha,beta,order,toll)
%LMF Metodo lineare multistep per il problema y'(t)=odefun(t,y(t)), y(t0)=y0.
% L'incognita y(t) è un vettore colonna di dimensione m.
% L'intervallo di tempo [t0,T] viene discretizzato con passo h uniforme.
% Il metodo è caratterizzato dai coefficienti alpha e beta.
% L'argomento order serve a calcolare in modo automatico una precisione
% adeguata con cui risolvere i sistemi di equazioni non lineari (utile
% solo nel caso di metodi impliciti). L'argomento facoltativo toll
% permette di stabilire manualmente tale precisione.
% L'argomento in uscita yT è l'approssimazione numerica di y(T).
% L'argomento in uscita tab è facoltativo, è di tipo "table" e contiene
% tutte le iterazioni effettuate dal metodo lineare multistep.

    % Controlla e normalizza i dati in ingresso
    assert(t0 <= T);
    assert(h > 0);
    N = round((T-t0)/h);
    h = (T-t0)/N;
    m = size(y0,1);
    assert(size(y0,2)==1);
    k = length(alpha)-1;
    assert(length(beta)==k+1);
    alpha_k = alpha(end);
    assert(alpha_k~=0);
    alpha = alpha / alpha_k;
    beta = beta / alpha_k;
    is_explicit = (beta(end)==0);
    assert(order >= 1);
    if nargin < 9
        toll = 0.1*(h^(order+1));  % un decimo dell'err. locale di troncamento,
        toll = max(toll,eps);      % ma non meno dell'epsilon di macchina
    end
    
    % Riserva memoria per i k passi più recenti.
    % La funzione idx serve a rendere circolari gli accessi a y e f
    y = zeros(m,k);
    f = zeros(m,k);
    idx = @(i) 1+mod(i,k);
    
    % Inizializza y e f con i valori iniziali da y_0 a y_{k-1}.
    % Se k > 1, i valori mancanti vengono calcolati con il metodo Runge-Kutta
    % esplicito ode45(), il cui errore di troncamento locale è O(h_rk^5).
    % Dato che l'errore sui valori iniziali si ripercuote su tutte le
    % iterazioni successive, è importante che questi vengano calcolati bene,
    % ossia con un errore dell'ordine di O(h^order)
    if k == 1
        y = y0;
    else
        tstep = [t0, t0+(k-1)*h];
        h_rk = min(h,h^(order/5));
        options = odeset('InitialStep',h_rk,'MaxStep',h_rk);
        sol = ode45(odefun,tstep,y0,options);
        y = deval(sol,t0+(0:k-1)*h);
    end
    for i = 1:k
        f(:,i) = odefun(t0+(i-1)*h,y(:,i));
    end
    
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
        for i = 1:k
            tab(i,1) = {t0 + (i-1)*h};
            tab(i,2:end) = array2table(y(:,i)');
        end
    end
    
    % Loop principale: a ogni iterazione viene calcolato y_n
    n = k;
    t = t0+k*h;
    while n <= N
        g_n = zeros(m,1);
        for i = 1:k
            g_n = g_n - alpha(i)*y(:,idx(n-k+i-1));
            g_n = g_n + h*beta(i)*f(:,idx(n-k+i-1));
        end
        if is_explicit
            y(:,idx(n)) = g_n;
        else
            y_predicted = y(:,idx(n-1)) + h*f(:,idx(n-1));
            [y(:,idx(n)),errcode] = LMF_nonlinear_solver(...
                odefun, t, y_predicted, h, beta(end), g_n, toll);
            if errcode == 1
                error('Il sistema nonlineare non è stato risolto.');
            end
        end
        f(:,idx(n)) = odefun(t,y(:,idx(n)));
        if save_table
            tab(n+1,1) = {t};
            tab(n+1,2:end) = array2table(y(:,idx(n))');
        end
        n = n + 1;
        t = t + h;
    end
    yT = y(:,idx(N));
end
