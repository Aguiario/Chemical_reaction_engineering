function Main_DME_DualType_EnhancedPlot()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ejemplo de código que combina:
    %  1) Búsqueda de L1, L2 bajo restricción L1 + L2 <= 18
    %  2) Resolución ODE de dos reactores secuenciales (dual) 
    %  3) Mejores gráficas (disposición y estética)
    %
    % NOTA: Ajustar parámetros de cinética, calores, coeficientes U, etc.,
    %       de acuerdo con el artículo y/o datos reales del proceso.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc; clear; close all;

    %% 1) Parámetros Iniciales
    F_CO = 100;       % Flujo inicial de CO (kmol/h)
    F_H2 = 200;       % Flujo inicial de H2 (kmol/h)
    T0   = 500;       % Temperatura inicial (K)
    P0   = 70;        % Presión inicial (bar)

    % Supuestos de cinética y transferencia (DEMOS):
    k1 = 0.1;         % Constante de reacción en el primer reactor
    k2 = 0.05;        % Constante de reacción en el segundo reactor
    U1 = 5;           % Coef. transferencia de calor en el primer reactor
    U2 = 3;           % Coef. transferencia de calor en el segundo reactor

    % Temperatura de refrigerante en cada reactor
    T_coolant1 = 450;  
    T_coolant2 = 470;

    % Rangos para L1 y L2 (en m)
    L1_vals = linspace(8, 12, 5); 
    L2_vals = linspace(5,  9,  5);

    %% 2) Variables para almacenar la mejor solución
    best_conversion = 0;
    best_L1 = 0;
    best_L2 = 0;

    % Contenedores para gráficas finales (guardar última corrida)
    z_total_saved    = [];
    F_total_saved    = [];
    T_total_saved    = [];
    P_total_saved    = [];
    X_CO_saved       = [];
    X_H2_saved       = [];
    F_DME_saved      = [];

    %% 3) Bucle: buscar la mejor combinación L1, L2
    for L1 = L1_vals
        for L2 = L2_vals

            % Checar restricción total
            if (L1 + L2) > 18
                continue;  % Descarta
            end

            % *********************
            % PRIMER REACTOR
            % *********************
            % Definimos ODEs p.ej: y = [FCO, FH2, T]
            % Ajustar la ecuación real:
            dydz1 = @(z,y) [
                -k1*y(1);                                    % Consumo de CO
                -2*k1*y(1);                                  % Consumo de H2 (ejemplo)
                0.3*y(1)*y(2) - U1*(y(3) - T_coolant1)        % Balance energía (ejemplo)
            ];

            % Cond. iniciales en z=0
            y0_1 = [F_CO; F_H2; T0];

            % Resolución ODE en [0, L1]
            [z1, y1] = ode45(dydz1, linspace(0, L1, 100), y0_1);

            % Salida del primer reactor
            F_CO_2in  = y1(end,1);  % CO que entra al 2do
            F_H2_2in  = y1(end,2);  % H2 que entra al 2do
            T1_out    = y1(end,3);

            % *********************
            % SEGUNDO REACTOR
            % *********************
            dydz2 = @(z, y) [
                -k2*y(1);                                    % Consumo de CO
                -1.5*k2*y(1);                                % Consumo de H2 (ejemplo)
                0.2*y(1)*y(2) - U2*(y(3) - T_coolant2)        % Balance energía (ejemplo)
            ];

            % Cond. iniciales en z=L1
            y0_2 = [F_CO_2in; F_H2_2in; T1_out];
            [z2, y2] = ode45(dydz2, linspace(L1, L1+L2, 100), y0_2);

            % *********************
            % Cálculo de conversiones
            % *********************
            F_CO_out = y2(:,1);
            F_H2_out = y2(:,2);

            X_CO = (F_CO - F_CO_out)./F_CO * 100;  % (%)
            X_H2 = (F_H2 - F_H2_out)./F_H2 * 100;  % (%)

            % Conversión final de CO
            X_CO_final = X_CO(end);

            % ¿Mejor?
            if X_CO_final > best_conversion
                best_conversion = X_CO_final;
                best_L1 = L1;
                best_L2 = L2;

                % Guardamos datos de la simulación para graficar
                z_total = [z1; z2];                 % Eje
                % Unimos y1, y2 en un solo contenedor
                y_total = [y1; y2];

                % Modelo sencillísimo de presión:
                % Por ejemplo: P(z) = P0 e^(-0.01*z).
                % Ajustado para no caer > 2 bar
                P_model = P0*exp(-0.01*z_total);
                P_limit = P0 - 2;      % no debe bajar de 68 bar
                P_model(P_model < P_limit) = P_limit;

                % Supongamos flujo de DME = (F_CO consumido)*0.95
                % (ficticio)
                F_DME = (F_CO - y2(:,1))*0.95; 

                z_total_saved = z_total;
                F_total_saved = y_total(:,1:2);  % CO,H2
                T_total_saved = y_total(:,3);
                P_total_saved = P_model;
                X_CO_saved    = X_CO;
                X_H2_saved    = X_H2;
                F_DME_saved   = F_DME;
            end
        end
    end

    %% 4) Impresión de la mejor configuración
    disp('=========================================================')
    disp('RESULTADOS:')
    disp([' - Mejor L1 = ', num2str(best_L1),' m'])
    disp([' - Mejor L2 = ', num2str(best_L2),' m'])
    disp([' - Max conversión de CO = ', num2str(best_conversion),' %'])
    disp([' - Longitud total = ', num2str(best_L1 + best_L2),' m'])
    disp('=========================================================')

    %% 5) Gráficas con mejor estética
    % Creamos una ventana de figura en pantalla completa:
    figure('Units','normalized','OuterPosition',[0 0 1 1]);

    % Eje Z final
    z_total = z_total_saved;
    F_CO_sol = F_total_saved(:,1);
    F_H2_sol = F_total_saved(:,2);
    T_sol    = T_total_saved;
    P_sol    = P_total_saved;
    F_DME_sol= F_DME_saved;
    X_CO_sol = X_CO_saved;
    X_H2_sol = X_H2_saved;

    % Subplot 1: Flujo de CO
    subplot(3,3,1)
    plot(z_total, F_CO_sol, 'r', 'LineWidth',1.5); grid on
    xlabel('Longitud (m)')
    ylabel('Flujo de CO (kmol/h)')
    title('Evolución de CO')
    legend('CO','Location','best')

    % Subplot 2: Flujo de H2
    subplot(3,3,2)
    plot(z_total, F_H2_sol, 'b', 'LineWidth',1.5); grid on
    xlabel('Longitud (m)')
    ylabel('Flujo de H2 (kmol/h)')
    title('Evolución de H2')
    legend('H2','Location','best')

    % Subplot 3: Temperatura
    subplot(3,3,3)
    plot(z_total, T_sol, 'g','LineWidth',1.5); grid on
    xlabel('Longitud (m)')
    ylabel('T (K)')
    title('Perfil de Temperatura')
    legend('Temperatura','Location','best')

    % Subplot 4: Presión
    subplot(3,3,4)
    plot(z_total, P_sol, 'm','LineWidth',1.5); grid on
    xlabel('Longitud (m)')
    ylabel('Presión (bar)')
    title('Perfil de Presión (modelo simple)')
    legend('Presión','Location','best')

    % Subplot 5: Flujo de DME
    subplot(3,3,5)
    plot(z_total, [zeros(length(z_total)-length(F_DME_sol),1); F_DME_sol], ...
         'c','LineWidth',1.5); 
    grid on
    xlabel('Longitud (m)')
    ylabel('Flujo DME (kmol/h)')
    title('Evolución de DME')
    legend('DME','Location','best')

    % Subplot 6: Conversión de CO
    subplot(3,3,6)
    % OJO: X_CO_sol solo está definido en la malla del 2do reactor
    % para simplificar la gráfica, repetimos 0 al inicio
    N1 = 100;  % del 1er reactor
    N2 = length(X_CO_sol);
    plot([z_total(1:N1); z_total(N1+1:end)], ...
         [zeros(N1,1); X_CO_sol], 'k','LineWidth',1.5);
    grid on
    xlabel('Longitud (m)')
    ylabel('X_{CO} (%)')
    title('Conversión de CO')
    legend('X_{CO}','Location','best')

    % Subplot 7: Conversión de H2
    subplot(3,3,7)
    plot([z_total(1:N1); z_total(N1+1:end)], ...
         [zeros(N1,1); X_H2_sol], 'Color',[0.5 0.2 0.9],'LineWidth',1.5);
    grid on
    xlabel('Longitud (m)')
    ylabel('X_{H2} (%)')
    title('Conversión de H2')
    legend('X_{H2}','Location','best')

    % Deja los últimos 2 subplots (3,8 y 3,9) vacíos o úsalos para más info
    subplot(3,3,8), axis off
    text(0,0.5,sprintf('Mejor L1=%.2f m\nMejor L2=%.2f m\nConv. CO=%.1f%%',...
         best_L1,best_L2,best_conversion),'FontSize',12);

    subplot(3,3,9), axis off
    %text(0,0.6,'Ejemplo de Modelo Simplificado','FontSize',14,'FontWeight','bold')
    %text(0,0.4,'Reemplazar con ecuaciones reales y cinética verificada','FontSize',10)
    
    F_DME = (F_CO - y2(:,1))*18;   % Ejemplo de “producción de DME”

    % ...
    % Subplot(5): Evolución de DME
    %subplot(7,1,5);
    %plot(z2, F_DME, '-c');

    % Supongamos que Mw_DME=46.07 g/mol
    Mw_DME = 46.07;        % g/mol
    final_F_DME_kmolh = F_DME(end);   % kmol/h (último valor)
    
    % Convertir a ton/día
    % 1) De kmol/h a kg/h => final_F_DME_kmolh * Mw_DME (g/mol)/1000
    % 2) Luego de kg/h a ton/día => multiplica por 24 y divide entre 1000
    FDME_ton_day = final_F_DME_kmolh * (Mw_DME / 1000) * 24;
    
    fprintf('Produccion final de DME: %.2f ton/dia\n', FDME_ton_day);

end
