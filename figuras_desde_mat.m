% figuras_desde_mat.m  (en la misma carpeta que res_viajero.mat, res_ackley.mat y res_inventario.mat)
V = load('res_viajero.mat'); A = load('res_ackley.mat'); I = load('res_inventario.mat');
D = {V, A, I};
nom     = {'Viajero (min)','Ackley (d=3)','Inventarios (Z)'};
probs   = {'viajero','ackley','inventario'};
metodos = {'elitismo','torneo','ruleta'};
sel     = {'Elitismo','Torneo','Ruleta'};

% ---------- seleccion: media +- desviacion ----------
figure('Color','w','Position',[100 100 1000 340]);
for p = 1:3
    m = zeros(1,3); s = m;
    for k = 1:3
        m(k) = D{p}.S1.(metodos{k}).media;
        s(k) = D{p}.S1.(metodos{k}).std;
    end
    subplot(1,3,p); bar(m); hold on; errorbar(1:3,m,s,'k.','LineWidth',1.2);
    set(gca,'XTick',1:3,'XTickLabel',sel); grid on;
    ylabel('Mejor fitness (media \pm desv.)'); title(nom{p});
end
exportgraphics(gcf,'fig_seleccion.png','Resolution',300);

% ---------- convergencia promedio ----------
figure('Color','w','Position',[100 100 1000 340]);
for p = 1:3
    subplot(1,3,p); hold on;
    for k = 1:3
        y = D{p}.S1.(metodos{k}).curva;
        if p == 2, y = y - (3-exp(1)) + 1e-12; end
        plot(0:numel(y)-1, y, 'LineWidth',1.4);
    end
    if p == 2, set(gca,'YScale','log'); ylabel('f - (3-e)'); else, ylabel('Mejor fitness'); end
    grid on; xlabel('Iteración'); title(nom{p});
    if p == 1, legend(sel,'Location','northeast'); end
end
exportgraphics(gcf,'fig_convergencia.png','Resolution',300);

% ---------- iteraciones (calidad y tiempo) ----------
figure('Color','w','Position',[100 100 1000 340]);
for p = 1:3
    R = D{p}.S2;
    subplot(1,3,p);
    yyaxis left;  errorbar(D{p}.cfg.iterLista,[R.media],[R.std],'-o','LineWidth',1.2); ylabel('Mejor fitness');
    yyaxis right; plot(D{p}.cfg.iterLista,[R.t],'-s','LineWidth',1.2); ylabel('Tiempo (s)');
    grid on; xlabel('Iteraciones'); title(nom{p});
end
exportgraphics(gcf,'fig_iteraciones.png','Resolution',300);

% ---------- probabilidad de mutacion ----------
figure('Color','w','Position',[100 100 1000 340]);
for p = 1:3
    R = D{p}.S3.pm;
    subplot(1,3,p);
    errorbar(D{p}.cfg.PmLista,[R.media],[R.std],'-o','LineWidth',1.2); grid on;
    xlabel('P_m'); ylabel('Mejor fitness'); title(nom{p});
end
exportgraphics(gcf,'fig_pm.png','Resolution',300);

% ---------- magnitud de la mutacion ----------
figure('Color','w','Position',[100 100 700 340]);
R = A.S3.sig;
subplot(1,2,1); errorbar(A.cfg.sigLista,[R.media],[R.std],'-o','LineWidth',1.2);
set(gca,'XScale','log'); grid on; xlabel('\sigma_m'); ylabel('Mejor fitness'); title('Ackley (d=3)');
R = I.S3.sig;
subplot(1,2,2); errorbar(I.cfg.fracLista,[R.media],[R.std],'-o','LineWidth',1.2);
set(gca,'XScale','log'); grid on; xlabel('\sigma_m / rango'); ylabel('Mejor fitness'); title('Inventarios');
exportgraphics(gcf,'fig_sigma.png','Resolution',300);

% ---------- dimensiones ----------
figure('Color','w','Position',[100 100 700 340]);
yyaxis left;  semilogx(A.cfg.dimLista,[A.S4.rmse],'-o','LineWidth',1.2); ylabel('RMSE');
yyaxis right; semilogx(A.cfg.dimLista,[A.S4.t],'-s','LineWidth',1.2); ylabel('Tiempo (s)');
grid on; xlabel('Dimensiones d'); title('Ackley: efecto de la dimensión');
exportgraphics(gcf,'fig_dimensiones.png','Resolution',300);

% ---------- individuos completos de Ackley (para el anexo) ----------
for k = 1:numel(A.cfg.dimLista)
    fprintf('d = %d, f = %.4f\n', A.cfg.dimLista(k), A.S4(k).mejor);
    fprintf('%.4f ', A.S4(k).ind); fprintf('\n\n');
end