% figuras_extra.m  Diagrama de flujo (fig_flujo.png) y funcion de Ackley en 1 y 2 dimensiones (fig_ackley.png)

% ---------- Diagrama de flujo ----------
figure('Color','w','Position',[100 100 420 640]);
axes('Position',[0 0 1 1]); axis off; hold on;
txt = {'Generar población inicial','Evaluar fitness', ...
       'Selección (elitismo, torneo o ruleta)','Cruce', ...
       'Mutación con probabilidad P_m','Evaluar nuevos individuos', ...
       '¿Se alcanzó el número de iteraciones?','Devolver el mejor individuo'};
y = 0.93:-0.125:0.055;
for k = 1:8
    rectangle('Position',[0.15 y(k)-0.035 0.7 0.07],'Curvature',0.3,'FaceColor',[0.92 0.95 1]);
    text(0.5,y(k),txt{k},'HorizontalAlignment','center','FontSize',9);
    if k < 8
        quiver(0.5,y(k)-0.035,0,-(y(k)-y(k+1)-0.07),0,'k','MaxHeadSize',1.2);
    end
end
plot([0.85 0.95 0.95],[y(7) y(7) y(3)],'k-');
quiver(0.95,y(3),-0.10,0,0,'k','MaxHeadSize',1.5);
text(0.96,mean([y(7) y(3)]),'No','Rotation',90,'HorizontalAlignment','center');
text(0.52,y(7)-0.062,'Sí');
xlim([0 1]); ylim([0 1]);
exportgraphics(gcf,'fig_flujo.png','Resolution',300);

% ---------- Función de Ackley en 1 y 2 dimensiones ----------
a = 20; b = 0.2; c = 2*pi;
f1 = @(x) -a*exp(-b*abs(x)) - exp(cos(c*x)) + a + 3;
x = linspace(-10,10,2000);
figure('Color','w','Position',[100 100 900 380]);
subplot(1,2,1); plot(x,f1(x),'LineWidth',1.2); grid on;
xlabel('x'); ylabel('f(x)'); title('(a) Ackley, d = 1');
[X,Y] = meshgrid(linspace(-10,10,200));
Z = -a*exp(-b*sqrt((X.^2+Y.^2)/2)) - exp((cos(c*X)+cos(c*Y))/2) + a + 3;
subplot(1,2,2); surf(X,Y,Z,'EdgeColor','none'); colorbar;
xlabel('x_1'); ylabel('x_2'); zlabel('f(x)'); title('(b) Ackley, d = 2'); view(35,35);
exportgraphics(gcf,'fig_ackley.png','Resolution',300);
