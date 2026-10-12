function ag_ackley(nRuns)
% AG_ACKLEY  Algoritmo genetico para la funcion de Ackley (d = 3 en los estudios 1 a 3).
% Estudios: 1 (seleccion), 2 (iteraciones), 3 (mutacion: Pm y sigma_m) y 4 (dimensiones).
% Uso:  ag_ackley              -> 10 corridas por escenario
%       ag_ackley(2)           -> prueba rapida
%       ag_ackley('verificar') -> comprueba la funcion objetivo
% Genera res_ackley.mat y muestra las tablas en pantalla.
if nargin >= 1 && ischar(nRuns), verificar(); return; end
if nargin < 1, nRuns = 10; end
rng(2026);

cfg.n = 50;  cfg.nRuns = nRuns;  cfg.iter0 = 100;
cfg.Pm0 = 0.1;  cfg.sig0 = 1;  cfg.d = 3;
cfg.iterLista = [10 20 30 50 100 150 200];
cfg.PmLista = [0 0.05 0.1 0.2 0.4 0.8];
cfg.sigLista = [0.1 0.5 1 2 5];
cfg.dimLista = [2 3 10 20 50 100];
ref = 3 - exp(1);           % optimo teorico de la version usada
metodos = {'elitismo','torneo','ruleta'};
prob = problema_ackley(cfg.d, cfg.sig0);

%% Estudio 1: metodo de seleccion
S1 = struct();
for k = 1:3
    S1.(metodos{k}) = correr(prob, cfg.n, cfg.iter0, cfg.Pm0, metodos{k}, nRuns, ref);
end
medias = [S1.elitismo.media, S1.torneo.media, S1.ruleta.media];
[~, ib] = min(medias);
mejorSel = metodos{ib};
Rs = S1.elitismo; Rs(2) = S1.torneo; Rs(3) = S1.ruleta;
imprimir_tabla(['Estudio 1 - Ackley d=' num2str(cfg.d) ' (mejor seleccion: ' mejorSel ')'], metodos, Rs);

%% Estudio 2: numero de iteraciones
S2 = []; etiq = cell(1, numel(cfg.iterLista));
for k = 1:numel(cfg.iterLista)
    r = correr(prob, cfg.n, cfg.iterLista(k), cfg.Pm0, mejorSel, nRuns, ref);
    if isempty(S2), S2 = r; else, S2(k) = r; end
    etiq{k} = sprintf('iter=%d', cfg.iterLista(k));
end
med = [S2.media];
k0 = find(med <= min(med)*1.05, 1, 'first');   % menor numero de iteraciones a menos de 5 % de la mejor media
iterOpt = cfg.iterLista(k0);
imprimir_tabla(['Estudio 2 - Ackley (iteraciones elegidas: ' num2str(iterOpt) ')'], etiq, S2);

%% Estudio 3: mutacion (primero Pm, luego sigma_m con la mejor Pm)
Rpm = []; etiqPm = cell(1, numel(cfg.PmLista));
for k = 1:numel(cfg.PmLista)
    r = correr(prob, cfg.n, iterOpt, cfg.PmLista(k), mejorSel, nRuns, ref);
    if isempty(Rpm), Rpm = r; else, Rpm(k) = r; end
    etiqPm{k} = sprintf('Pm=%.2f', cfg.PmLista(k));
end
[~, ib] = min([Rpm(2:end).media]);   % se excluye Pm=0 para que el estudio de sigma tenga sentido
pmOpt = cfg.PmLista(ib + 1);
Rsg = []; etiqSg = cell(1, numel(cfg.sigLista));
for k = 1:numel(cfg.sigLista)
    r = correr(problema_ackley(cfg.d, cfg.sigLista(k)), cfg.n, iterOpt, pmOpt, mejorSel, nRuns, ref);
    if isempty(Rsg), Rsg = r; else, Rsg(k) = r; end
    etiqSg{k} = sprintf('sigma=%.2f', cfg.sigLista(k));
end
[~, ib] = min([Rsg.media]);
sigOpt = cfg.sigLista(ib);
S3 = struct(); S3.pm = Rpm; S3.sig = Rsg;
imprimir_tabla(['Estudio 3 - Ackley (Pm=' num2str(pmOpt) ', sigma=' num2str(sigOpt) ')'], [etiqPm, etiqSg], [Rpm, Rsg]);

%% Estudio 4: numero de dimensiones
S4 = []; etiq4 = cell(1, numel(cfg.dimLista));
for k = 1:numel(cfg.dimLista)
    d = cfg.dimLista(k);
    r = correr(problema_ackley(d, sigOpt), cfg.n, iterOpt, pmOpt, mejorSel, nRuns, ref);
    if isempty(S4), S4 = r; else, S4(k) = r; end
    etiq4{k} = sprintf('d=%d', d);
end
imprimir_tabla('Estudio 4 - Ackley vs dimensiones', etiq4, S4);

save('res_ackley.mat', 'cfg', 'S1', 'S2', 'S3', 'S4', 'mejorSel', 'iterOpt', 'pmOpt', 'sigOpt');
fprintf('\nResultados guardados en res_ackley.mat\n');
end

%% ===================== FUNCION OBJETIVO =====================
function f = fobj_ackley(x)
% Ackley con a=20, b=0.2, c=2*pi y constante final 3 en lugar de e.
a = 20; b = 0.2; c = 2*pi; d = numel(x);
f = -a*exp(-b*sqrt(sum(x.^2)/d)) - exp(sum(cos(c*x))/d) + a + 3;
end

function verificar()
fprintf('Ackley en el origen (d=3): %.6f (esperado 3-e = %.6f)\n', fobj_ackley([0 0 0]), 3 - exp(1));
end

%% ===================== POBLACION, CRUCE Y MUTACION =====================
function P = pob_ackley(n, d, lim)
P = lim(1) + (lim(2) - lim(1))*rand(n, d);
end

function H = cruce_blx(padres, nh, lim, alfa)
% Cruce BLX-alfa para vectores reales
np = size(padres,1); L = size(padres,2); H = zeros(nh, L);
for k = 1:nh
    ij = randperm(np, 2); p1 = padres(ij(1),:); p2 = padres(ij(2),:);
    lo = min(p1, p2); hi = max(p1, p2); dd = hi - lo;
    h = (lo - alfa*dd) + rand(1, L).*(dd*(1 + 2*alfa));
    H(k,:) = min(max(h, lim(1)), lim(2));
end
end

function x = mut_ackley(x, sigma, lim)
% Mutacion gaussiana sobre un gen elegido al azar
j = randi(numel(x));
x(j) = min(max(x(j) + sigma*randn, lim(1)), lim(2));
end

function prob = problema_ackley(d, sigma)
lim = [-10 10];
prob.fobj = @(x) fobj_ackley(x);
prob.pob = @(n) pob_ackley(n, d, lim);
prob.cruce = @(P, nh) cruce_blx(P, nh, lim, 0.5);
prob.mut = @(x) mut_ackley(x, sigma, lim);
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

function res = correr(prob, n, iter, Pm, sel, nRuns, refVal)
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
res.rmse = sqrt(mean((Y - refVal).^2));      % RMSE respecto a 3-e
res.Y = Y;
res.curva = mean(H, 2);
end

function s = fmt_ind(x)
if numel(x) > 5
    s = ['[' sprintf('%.4f ', x(1:5)) '... (d=' num2str(numel(x)) ')]'];
else
    s = ['[' sprintf('%.4f ', x) ']'];
end
end

function imprimir_tabla(titulo, etiq, Rs)
fprintf('\n%s\n', titulo);
fprintf('| Configuracion | Tiempo (s) | Mejor f(X) | Media | Desv. est. | RMSE | Mejor individuo |\n');
fprintf('|---|---|---|---|---|---|---|\n');
for k = 1:numel(Rs)
    r = Rs(k);
    fprintf('| %s | %.3f | %.4f | %.4f | %.4f | %.4f | %s |\n', etiq{k}, r.t, r.mejor, r.media, r.std, r.rmse, fmt_ind(r.ind));
end
end