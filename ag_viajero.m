function ag_viajero(nRuns)
% AG_VIAJERO  Algoritmo genetico para el problema del viajero (13 ubicaciones).
% Estudios: 1 (seleccion), 2 (iteraciones) y 3 (probabilidad de mutacion).
% Uso:  ag_viajero              -> 10 corridas por escenario
%       ag_viajero(2)           -> prueba rapida
%       ag_viajero('verificar') -> comprueba la funcion objetivo
% Genera res_viajero.mat y muestra las tablas en pantalla.
if nargin >= 1 && ischar(nRuns), verificar(); return; end
if nargin < 1, nRuns = 10; end
rng(2026);

cfg.n = 50;                 % tamano de la poblacion
cfg.nRuns = nRuns;          % corridas por escenario
cfg.iter0 = 100;            % iteraciones del estudio 1
cfg.Pm0 = 0.1;              % probabilidad de mutacion base
cfg.iterLista = [10 20 30 50 100 150 200];
cfg.PmLista = [0 0.05 0.1 0.2 0.4 0.8];
metodos = {'elitismo','torneo','ruleta'};
prob = problema_viajero();

%% Estudio 1: metodo de seleccion
S1 = struct();
for k = 1:3
    S1.(metodos{k}) = correr(prob, cfg.n, cfg.iter0, cfg.Pm0, metodos{k}, nRuns);
end
medias = [S1.elitismo.media, S1.torneo.media, S1.ruleta.media];
[~, ib] = min(medias);
mejorSel = metodos{ib};
Rs = S1.elitismo; Rs(2) = S1.torneo; Rs(3) = S1.ruleta;
imprimir_tabla(['Estudio 1 - viajero (mejor seleccion: ' mejorSel ')'], metodos, Rs);

%% Estudio 2: numero de iteraciones
S2 = []; etiq = cell(1, numel(cfg.iterLista));
for k = 1:numel(cfg.iterLista)
    r = correr(prob, cfg.n, cfg.iterLista(k), cfg.Pm0, mejorSel, nRuns);
    if isempty(S2), S2 = r; else, S2(k) = r; end
    etiq{k} = sprintf('iter=%d', cfg.iterLista(k));
end
med = [S2.media];
% criterio: menor numero de iteraciones con media a menos de 5 % de la mejor media
k0 = find(med <= min(med)*1.05, 1, 'first');
iterOpt = cfg.iterLista(k0);
imprimir_tabla(['Estudio 2 - viajero (iteraciones elegidas: ' num2str(iterOpt) ')'], etiq, S2);

%% Estudio 3: probabilidad de mutacion
Rpm = []; etiq = cell(1, numel(cfg.PmLista));
for k = 1:numel(cfg.PmLista)
    r = correr(prob, cfg.n, iterOpt, cfg.PmLista(k), mejorSel, nRuns);
    if isempty(Rpm), Rpm = r; else, Rpm(k) = r; end
    etiq{k} = sprintf('Pm=%.2f', cfg.PmLista(k));
end
[~, ib] = min([Rpm.media]);
pmOpt = cfg.PmLista(ib);
S3 = struct(); S3.pm = Rpm;
imprimir_tabla(['Estudio 3 - viajero (Pm elegida: ' num2str(pmOpt) ')'], etiq, Rpm);

save('res_viajero.mat', 'cfg', 'S1', 'S2', 'S3', 'mejorSel', 'iterOpt', 'pmOpt');
fprintf('\nResultados guardados en res_viajero.mat\n');
end

%% ===================== FUNCION OBJETIVO =====================
function D = matriz_tiempos()
% Tabla I (minutos). 1 UCO, 2 Rionegro, 3 La Ceja, 4 El Carmen, 5 Marinilla,
% 6 Guarne, 7 El Retiro, 8 Santuario, 9 La Union, 10 San Vicente, 11 Cocorna,
% 12 El Penol, 13 Guatape. D(i,j) = tiempo de i hacia j.
D = [  0  10  35  29  16  30  40  27  56  44  67  54  81
       6   0  32  31  16  32  34  28  52  44  63  56  83
      34  32   0  29  45  58  35  57  26  73  91  83 110
      29  29  28   0  32  47  57  31  49  61  67  70  99
      18  18  46  33   0  23  52  19  68  39  54  41  69
      30  30  58  46  25   0  46  36  79  39  71  62  90
      40  37  36  59  52  54   0  64  61  80 104  89 117
      32  32  57  32  18  36  65   0  82  50  41  51  82
      55  52  23  48  67  80  56  73   0  95 116 101 130
      47  47  76  62  42  37  73  53  97   0 102  81 111
      79  79 106  80  66  87 110  59 148 101   0  97 100
      59  59  86  74  44  65  90  55 107  80  93   0  29
      89  89 116 104  74  96 119  86 137 110 103  33   0 ];
end

function T = fobj_viajero(ruta, regreso)
% Tiempo total (min) de un recorrido. ruta debe ser permutacion de 1..13.
if nargin < 2, regreso = false; end
persistent D
if isempty(D), D = matriz_tiempos(); end
ruta = ruta(:)';
if numel(ruta) ~= 13 || ~isequal(sort(ruta), 1:13)
    T = 10000 + 100*(13 - numel(unique(ruta)));   % penalizacion
    return
end
T = sum(D(sub2ind([13 13], ruta(1:end-1), ruta(2:end))));
if regreso, T = T + D(ruta(end), ruta(1)); end
end

function verificar()
fprintf('Ruta 1:13 sin regreso: %d (esperado 641)\n', fobj_viajero(1:13, false));
fprintf('Ruta 1:13 con regreso: %d (esperado 730)\n', fobj_viajero(1:13, true));
fprintf('Ruta invalida: %d (esperado un valor mayor a 10000)\n', fobj_viajero([1 1 3:13], false));
end

%% ===================== POBLACION, CRUCE Y MUTACION =====================
function P = pob_viajero(n)
P = zeros(n, 13);
for i = 1:n, P(i,:) = randperm(13); end
end

function H = cruce_ox(padres, nh)
% Cruce de orden (OX): copia un segmento de un padre y completa con el orden del otro
np = size(padres,1); L = size(padres,2); H = zeros(nh, L);
for k = 1:nh
    ij = randperm(np, 2); p1 = padres(ij(1),:); p2 = padres(ij(2),:);
    c = sort(randperm(L, 2)); a = c(1); b = c(2);
    h = zeros(1, L); h(a:b) = p1(a:b);
    p2r = p2([b+1:L, 1:b]);
    resto = p2r(~ismember(p2r, h(a:b)));
    h([b+1:L, 1:a-1]) = resto;
    H(k,:) = h;
end
end

function x = mut_viajero(x)
% Mutacion por intercambio de dos posiciones
ij = randperm(numel(x), 2); x(ij) = x(ij([2 1]));
end

function prob = problema_viajero()
prob.fobj = @(x) fobj_viajero(x, false);
prob.pob = @(n) pob_viajero(n);
prob.cruce = @(P, nh) cruce_ox(P, nh);
prob.mut = @(x) mut_viajero(x);
end

%% ===================== NUCLEO DEL ALGORITMO GENETICO =====================
function idx = seleccionar(metodo, fit, m)
fit = fit(:); n = numel(fit);
switch metodo
    case 'elitismo'
        [~, o] = sort(fit, 'ascend'); idx = o(1:m);
    case 'torneo'
        k = 4; nGan = 2; pool = 1:n; idx = [];
        while numel(idx) < m
            kk = min(k, numel(pool));
            g = pool(randperm(numel(pool), kk));
            [~, o] = sort(fit(g), 'ascend');
            gan = g(o(1:min(nGan, kk)));
            gan = gan(1:min(numel(gan), m - numel(idx)));
            idx = [idx, gan];
            pool = setdiff(pool, gan);
        end
        idx = idx(:);
    case 'ruleta'
        p = 1./fit; libre = true(n,1); idx = zeros(m,1);
        for j = 1:m
            q = p.*libre; c = cumsum(q)/sum(q); c(end) = 1;
            i = find(rand <= c, 1, 'first');
            idx(j) = i; libre(i) = false;
        end
    otherwise
        error('Metodo de seleccion desconocido');
end
end

function f = evaluar(fobj, P)
f = zeros(size(P,1), 1);
for i = 1:size(P,1), f(i) = fobj(P(i,:)); end
end

function R = ag(prob, n, iter, Pm, metodoSel)
t0 = tic;
P = prob.pob(n); fit = evaluar(prob.fobj, P);
[mejorFit, ib] = min(fit); mejor = P(ib,:);
hist = zeros(iter+1, 1); hist(1) = mejorFit;
m = floor(n/2); nh = n - m;
for g = 1:iter
    idx = seleccionar(metodoSel, fit, m);
    padres = P(idx,:);
    hijos = prob.cruce(padres, nh);
    for k = 1:nh
        if rand < Pm, hijos(k,:) = prob.mut(hijos(k,:)); end
    end
    fh = evaluar(prob.fobj, hijos);
    P = [padres; hijos]; fit = [fit(idx); fh];
    [fmin, ib] = min(fit);
    if fmin < mejorFit
        mejorFit = fmin; mejor = P(ib,:);
    elseif fmin > mejorFit            % el mejor historico se perdio: se reinserta
        [~, iw] = max(fit); P(iw,:) = mejor; fit(iw) = mejorFit;
    end
    hist(g+1) = mejorFit;
end
R.mejor = mejor; R.fit = mejorFit; R.hist = hist; R.t = toc(t0);
end

function res = correr(prob, n, iter, Pm, sel, nRuns)
Y = zeros(nRuns,1); T = Y; H = zeros(iter+1, nRuns); ind = cell(nRuns,1);
for r = 1:nRuns
    R = ag(prob, n, iter, Pm, sel);
    Y(r) = R.fit; T(r) = R.t; H(:,r) = R.hist; ind{r} = R.mejor;
end
[res.mejor, ib] = min(Y);
res.ind = ind{ib};
res.media = mean(Y);
res.std = std(Y);
res.t = mean(T);
res.Y = Y;
res.curva = mean(H, 2);
end

function imprimir_tabla(titulo, etiq, Rs)
fprintf('\n%s\n', titulo);
fprintf('| Configuracion | Tiempo (s) | Mejor f(X) | Media | Desv. est. | Mejor individuo |\n');
fprintf('|---|---|---|---|---|---|\n');
for k = 1:numel(Rs)
    r = Rs(k);
    fprintf('| %s | %.3f | %.4f | %.4f | %.4f | [%s] |\n', etiq{k}, r.t, r.mejor, r.media, r.std, sprintf('%d ', r.ind));
end
end