function ag_inventario(nRuns)
% AG_INVENTARIO  Algoritmo genetico para la optimizacion de inventarios (5 productos).
% Estudios: 1 (seleccion), 2 (iteraciones) y 3 (mutacion: Pm y sigma_m).
% Uso:  ag_inventario              -> 10 corridas por escenario
%       ag_inventario(2)           -> prueba rapida
%       ag_inventario('verificar') -> compara la forma cerrada con integral
% Individuo: 14 variables [Q1..Q4 S1..S5 K1..K5]; Q5 se despeja de la Ecuacion 9.
% Genera res_inventario.mat y muestra las tablas en pantalla.
if nargin >= 1 && ischar(nRuns), verificar(); return; end
if nargin < 1, nRuns = 10; end
rng(2026);

cfg.n = 50;  cfg.nRuns = nRuns;  cfg.iter0 = 100;
cfg.Pm0 = 0.1;  cfg.frac0 = 0.1;      % sigma_m como fraccion del rango de la variable
cfg.iterLista = [10 20 30 50 100 150 200];
cfg.PmLista = [0 0.05 0.1 0.2 0.4 0.8];
cfg.fracLista = [0.01 0.05 0.1 0.25 0.5];
metodos = {'elitismo','torneo','ruleta'};
prob = problema_inventario(cfg.frac0);

%% Estudio 1: metodo de seleccion
S1 = struct();
for k = 1:3
    S1.(metodos{k}) = correr(prob, cfg.n, cfg.iter0, cfg.Pm0, metodos{k}, nRuns);
end
medias = [S1.elitismo.media, S1.torneo.media, S1.ruleta.media];
[~, ib] = min(medias);
mejorSel = metodos{ib};
Rs = S1.elitismo; Rs(2) = S1.torneo; Rs(3) = S1.ruleta;
imprimir_tabla(['Estudio 1 - inventarios (mejor seleccion: ' mejorSel ')'], metodos, Rs);

%% Estudio 2: numero de iteraciones
S2 = []; etiq = cell(1, numel(cfg.iterLista));
for k = 1:numel(cfg.iterLista)
    r = correr(prob, cfg.n, cfg.iterLista(k), cfg.Pm0, mejorSel, nRuns);
    if isempty(S2), S2 = r; else, S2(k) = r; end
    etiq{k} = sprintf('iter=%d', cfg.iterLista(k));
end
med = [S2.media];
k0 = find(med <= min(med)*1.05, 1, 'first');   % menor numero de iteraciones a menos de 5 % de la mejor media
iterOpt = cfg.iterLista(k0);
imprimir_tabla(['Estudio 2 - inventarios (iteraciones elegidas: ' num2str(iterOpt) ')'], etiq, S2);

%% Estudio 3: mutacion (primero Pm, luego sigma_m con la mejor Pm)
Rpm = []; etiqPm = cell(1, numel(cfg.PmLista));
for k = 1:numel(cfg.PmLista)
    r = correr(prob, cfg.n, iterOpt, cfg.PmLista(k), mejorSel, nRuns);
    if isempty(Rpm), Rpm = r; else, Rpm(k) = r; end
    etiqPm{k} = sprintf('Pm=%.2f', cfg.PmLista(k));
end
[~, ib] = min([Rpm.media]);
pmOpt = cfg.PmLista(ib);
Rsg = []; etiqSg = cell(1, numel(cfg.fracLista));
for k = 1:numel(cfg.fracLista)
    r = correr(problema_inventario(cfg.fracLista(k)), cfg.n, iterOpt, pmOpt, mejorSel, nRuns);
    if isempty(Rsg), Rsg = r; else, Rsg(k) = r; end
    etiqSg{k} = sprintf('sigma=%.2f', cfg.fracLista(k));
end
[~, ib] = min([Rsg.media]);
fracOpt = cfg.fracLista(ib);
S3 = struct(); S3.pm = Rpm; S3.sig = Rsg;
imprimir_tabla(['Estudio 3 - inventarios (Pm=' num2str(pmOpt) ', sigma=' num2str(fracOpt) ')'], [etiqPm, etiqSg], [Rpm, Rsg]);

save('res_inventario.mat', 'cfg', 'S1', 'S2', 'S3', 'mejorSel', 'iterOpt', 'pmOpt', 'fracOpt');
fprintf('\nResultados guardados en res_inventario.mat\n');
end

%% ===================== FUNCION OBJETIVO =====================
function p = params_inventario()
persistent P
if isempty(P)
    P.D  = [12000 80000 150000 50000 100000];
    P.m  = [20 220 900 120 180];
    P.mu = [300 4000 7000 2000 3500];
    P.sg = [100 1200 2500 700 1100];
    P.Smax = [200 1200 2500 800 1500];
    P.Kmax = [500 2000 5000 1500 3000];
    P.I = 8000;
    P.Qmin = 10;          % valor positivo (el intervalo es (0, inf))
    P.Qmax = P.I/2;       % supuesto prudente en lugar de infinito
end
p = P;
end

function Z = fobj_inventario(Q, S, K, metodo)
% Ecuacion 8: Z = sum_i D_i/(Q_i m_i) * integral_{K_i+S_i}^{inf} (x-R_i) f(x) dx
% metodo 'analitico' (forma cerrada, rapida) o 'integral' (funcion integral de MATLAB)
if nargin < 4, metodo = 'analitico'; end
p = params_inventario();
R = K + S; Z = 0;
for i = 1:5
    if strcmp(metodo, 'integral')
        g = @(x) (x - R(i)).*exp(-(x - p.mu(i)).^2/(2*p.sg(i)^2))/(p.sg(i)*sqrt(2*pi));
        J = integral(g, R(i), Inf);
    else    % forma cerrada de E[(X-R)+] para X ~ N(mu, sigma)
        z = (R(i) - p.mu(i))/p.sg(i);
        J = p.sg(i)*(exp(-z^2/2)/sqrt(2*pi) - z*0.5*erfc(z/sqrt(2)));
    end
    Z = Z + p.D(i)/(Q(i)*p.m(i))*J;
end
end

function verificar()
Q = [2000 3000 1000 1500 1000]; S = [100 600 1000 400 700]; K = [400 1500 3000 1000 2000];
fprintf('Z forma cerrada = %.6f\n', fobj_inventario(Q, S, K, 'analitico'));
fprintf('Z con integral  = %.6f   (ambos deben coincidir, cerca de 1473.0238)\n', fobj_inventario(Q, S, K, 'integral'));
end

function q5 = q5_inventario(x)
% Ecuacion 9: sum(Q_i/2 + S_i) = I  ->  se despeja Q5 de las otras 14 variables
p = params_inventario();
q5 = 2*(p.I - sum(x(5:9)) - sum(x(1:4))/2);
end

function Z = fit_inventario(x)
% x = [Q1..Q4, S1..S5, K1..K5]
p = params_inventario();
q5 = q5_inventario(x);
if q5 < p.Qmin - 1e-9 || q5 > p.Qmax + 1e-9, Z = 1e6; return; end
Z = fobj_inventario([x(1:4) q5], x(5:9), x(10:14), 'analitico');
end

function [lb, ub] = limites_inv()
p = params_inventario();
lb = [repmat(p.Qmin,1,4), zeros(1,10)];
ub = [repmat(p.Qmax,1,4), p.Smax, p.Kmax];
end

%% ===================== POBLACION, CRUCE Y MUTACION =====================
function P = pob_inventario(n)
p = params_inventario(); [lb, ub] = limites_inv();
P = zeros(n, 14); k = 0;
while k < n
    x = lb + rand(1,14).*(ub - lb);
    q5 = q5_inventario(x);
    if q5 >= p.Qmin && q5 <= p.Qmax, k = k + 1; P(k,:) = x; end
end
end

function H = cruce_aritmetico(padres, nh)
% Cruce aritmetico: un unico lambda por hijo. Como la Ecuacion 9 es lineal,
% el hijo hereda la factibilidad de los padres.
np = size(padres,1); H = zeros(nh, size(padres,2));
for k = 1:nh
    ij = randperm(np, 2); lam = rand;
    H(k,:) = lam*padres(ij(1),:) + (1 - lam)*padres(ij(2),:);
end
end

function x = mut_inventario(x, frac)
% Mutacion gaussiana sobre una variable; sigma_m = frac * rango; Q5 se recalcula
[lb, ub] = limites_inv(); p = params_inventario();
j = randi(14); y = x;
y(j) = min(max(x(j) + frac*(ub(j) - lb(j))*randn, lb(j)), ub(j));
q5 = q5_inventario(y);
if q5 >= p.Qmin && q5 <= p.Qmax, x = y; end   % si Q5 queda inviable se descarta
end

function prob = problema_inventario(frac)
prob.fobj = @(x) fit_inventario(x);
prob.pob = @(n) pob_inventario(n);
prob.cruce = @(P, nh) cruce_aritmetico(P, nh);
prob.mut = @(x) mut_inventario(x, frac);
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

function s = fmt_ind(x)
x15 = [x(1:4) q5_inventario(x) x(5:14)];
s = ['Q=[' sprintf('%.0f ', x15(1:5)) '] S=[' sprintf('%.0f ', x15(6:10)) '] K=[' sprintf('%.0f ', x15(11:15)) ']'];
end

function imprimir_tabla(titulo, etiq, Rs)
fprintf('\n%s\n', titulo);
fprintf('| Configuracion | Tiempo (s) | Mejor Z | Media | Desv. est. | Mejor individuo |\n');
fprintf('|---|---|---|---|---|---|\n');
for k = 1:numel(Rs)
    r = Rs(k);
    fprintf('| %s | %.3f | %.4f | %.4f | %.4f | %s |\n', etiq{k}, r.t, r.mejor, r.media, r.std, fmt_ind(r.ind));
end
end