# Algoritmos genéticos: viajero, Ackley e inventarios

Trabajo de Mateo Guerra Quintero, Bioingeniería, Unidad 3 (Optimización metaheurística), Universidad Católica de Oriente.

Implementación en MATLAB de un algoritmo genético aplicado a tres problemas de minimización: el problema del viajero (13 ubicaciones), la función de Ackley y un problema de inventarios con restricción de presupuesto.

## Archivos

- `ag_viajero.m`: problema del viajero (cruce OX, mutación por intercambio).
- `ag_ackley.m`: función de Ackley (cruce BLX-α, mutación gaussiana).
- `ag_inventario.m`: inventarios (cruce aritmético, mutación gaussiana).
- `figuras_desde_mat.m`: genera las figuras de resultados a partir de los `.mat`.
- `figuras_extra.m`: diagrama de flujo y función de Ackley en 1 y 2 dimensiones.
- `salida_*.txt`: tablas impresas por cada corrida real.

## Cómo ejecutarlo

Requiere MATLAB (`exportgraphics` necesita R2020a o posterior). Con los archivos en una misma carpeta:

    ag_viajero('verificar')      % comprueba la función objetivo
    ag_viajero                   % corrida real (10 corridas por escenario)
    ag_ackley
    ag_inventario
    figuras_desde_mat            % figuras de resultados
    figuras_extra                % diagrama de flujo y figura de Ackley

Cada script guarda sus resultados en `res_*.mat` e imprime las tablas en pantalla. Con `ag_...(2)` se hace una prueba rápida de 2 corridas.

## Configuración

Población de 50 individuos, 10 corridas por escenario y semilla `rng(2026)`.
