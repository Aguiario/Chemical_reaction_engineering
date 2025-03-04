function cstr_control_sim()
    clear; clc; close all;

    %% Parámetros del reactor y cinética
    params.V       = 1.0;       % Volumen del reactor (m^3)
    params.k0      = 7;     % Factor pre-exponencial (1/min) (AUMENTADO)
    params.E       = 8000;      % Energía de activación (J/mol)
    params.R       = 8.314;     % Constante de gases (J/mol/K)
    params.T       = 350;       % Temperatura isoterma en el reactor (K)
    params.k_inh   = 0.5;       % Factor del inhibidor
    params.Cx_feed = 0.1;       % Concentración inhibidor (kmol/m^3)
    params.CA_feed = 1.0;       % Concentración de entrada de A (kmol/m^3)

    % Caudal de alimentación nominal
    F_nominal      = 2;       % (m^3/min) (AUMENTADO PARA MÁS FLUJO)

    % Parámetro de válvula
    params.Kv      = 2.0;       % Factor de la válvula (AUMENTADO)
    P2_nominal     = 10;        % Presión nominal (bar)

    %% Parámetros de control PID
    % Lazo interno (control de flujo F_A)
    Kc_F   = 3.0;      % Incrementamos la ganancia para mejorar respuesta
    tauI_F = 5.0;      % Ajustamos el tiempo integral
    tauD_F = 0.0;      

    % Lazo externo (control de concentración C_A)
    Kc_CA   = -5.0;    % Se ajusta para respuesta más rápida
    tauI_CA = 20.0;    % Se acorta el tiempo integral
    tauD_CA = 0.0;     

    % Inicialización de controladores
    ctrlF.integral = 0;
    ctrlF.pv_prev  = 0;
    ctrlCA.integral = 0;
    ctrlCA.pv_prev  = 0;

    %% Configuración de la simulación
    tSim = 50;         
    dt   = 0.1;        
    time  = (0:dt:tSim);

    % Estado inicial de la concentración de A en el reactor
    CA0     = 0.2;     % Se inicia con más reactivo
    xk      = CA0;    

    % Setpoints
    SP_CA = 0.5;       
    SP_F  = F_nominal; 

    % Perturbación en P2
    P2_profile = @(t) ( t<10 )*P2_nominal + ( t>=10 )*(P2_nominal + 2.0 );

    %% Simulación con integración ODE
    CA_history    = zeros(length(time),1);
    F_history     = zeros(length(time),1);
    SP_F_history  = zeros(length(time),1);
    P2_history    = zeros(length(time),1);
    
    CA_history(1)   = xk;      
    F_history(1)    = F_nominal;
    SP_F_history(1) = SP_F;
    P2_history(1)   = P2_profile(0);

    % Usamos ODE45 para mejorar la estabilidad
    options = odeset('RelTol',1e-6,'AbsTol',1e-8);
    
    for k = 1:length(time)-1
        tk   = time(k);
        tk1  = time(k+1);
        
        % Medición de C_A
        CA_meas = xk;
        
        % Lazo externo: Control de C_A -> Setpoint de F_A
        [SP_F_new, ctrlCA] = pid_controller(SP_CA, CA_meas, ctrlCA, Kc_CA, tauI_CA, tauD_CA, dt);
        SP_F_new = max(0.5, min(5.0, SP_F_new));  % Limitamos el rango
        
        % Suavizamos la referencia de F_A con filtro de primer orden
        alpha = 0.1;  
        SP_F = alpha * SP_F_new + (1 - alpha) * SP_F;

        % Lazo interno: Control de F_A
        F_meas = getFlow(SP_F, P2_profile(tk), params);
        [MV_F, ctrlF] = pid_controller(SP_F, F_meas, ctrlF, Kc_F, tauI_F, tauD_F, dt);
        MV_F = max(0.0, min(100.0, MV_F));  

        % Flujo real F_A
        F_actual = getFlow(MV_F, P2_profile(tk), params);
        
        % Integración del modelo usando ODE45
        [~, CA_next] = ode45(@(t,C) cstr_ode(t, C, F_actual, params), [tk tk1], xk, options);
        xk = CA_next(end);

        % Guardar datos
        CA_history(k+1)   = xk;
        F_history(k+1)    = F_actual;
        SP_F_history(k+1) = SP_F;
        P2_history(k+1)   = P2_profile(tk1);
    end
    
    %% Graficar resultados
    figure; 
    subplot(3,1,1)
       plot(time, P2_history, 'LineWidth',2); grid on
       xlabel('Tiempo (min)'); ylabel('P_2 (bar)');
       title('Perturbación en la presión de alimentación')
    subplot(3,1,2)
       plot(time, CA_history, 'b-','LineWidth',2); hold on;
       yline(SP_CA,'r--','SP C_A','LineWidth',1.5);
       grid on; xlabel('Tiempo (min)'); ylabel('C_A (kmol/m^3)');
       title('Concentración de A en el reactor')
    subplot(3,1,3)
       plot(time, F_history,'LineWidth',2); hold on;
       plot(time, SP_F_history,'r--','LineWidth',1.5);
       grid on; xlabel('Tiempo (min)'); ylabel('F_A (m^3/min)');
       legend('F_A real','SP F_A','Location','Best');
       title('Caudal de alimentación de A')
    
end

%% Subfunción: Dinámica del reactor
function dCA_dt = cstr_ode(t, CA, F, p)
    V     = p.V;
    k0    = p.k0;
    E     = p.E;
    R     = p.R;
    T     = p.T;
    k_inh = p.k_inh;
    Cx    = p.Cx_feed;
    CAin  = p.CA_feed;

    % Cinética de reacción
    rA = -k0 * exp(-E/(R*T)) * (CA / (1 + k_inh * Cx));
    
    % Balance de materia
    dCA_dt = ( F*(CAin - CA) + V*rA ) / V;
end

%% Subfunción: Relación entre válvula y flujo F_A
function F = getFlow(valveSignal, P2, p)
    Kv = p.Kv;
    frac = valveSignal / 100;
    F = Kv * frac * sqrt( max(P2,0) );
end

%% Subfunción: Control PID
function [MV, ctrl] = pid_controller(SP, PV, ctrl, Kc, tauI, tauD, dt)
    error = SP - PV;
    
    P = Kc * error;
    
    if tauI > 1e-6
        ctrl.integral = ctrl.integral + (Kc/tauI)*error*dt;
    end
    I = ctrl.integral;
    
    D = 0;
    if tauD > 1e-6
        dPV = (PV - ctrl.pv_prev)/dt;
        D = -Kc * tauD * dPV;
    end
    
    MV = P + I + D;
    ctrl.pv_prev = PV;
end
