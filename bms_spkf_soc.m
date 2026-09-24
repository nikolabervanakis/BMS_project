
% Estimacija stanja napunjenosti (SoC) Li-ion celije pomocu SPKF-a
% Referenca: Plett ECE5550 (Lekcija 6.6)


clear; clc; clf;

%% 1. Parametri baterije i simulacije
Cap_Ah   = 2.5;                     % nominalni kapacitet celije [Ah]
Q_total  = Cap_Ah * 3600;           % kapacitet u As
R0       = 0.035;                   % unutarnji omski otpor [Ohm]
dt       = 1;                       % period uzorkovanja [s]
maxIter  = 300;                     % trajanje simulacije (5 min)

% Aproksimacija nelinearne OCV-SoC krivulje
ocv_func = @(s) 3.0 + 0.9*s + 0.25*s.^2 + 0.05*log(s + 1e-4) - 0.03*log(1.001 - s);

%% 2. Postavke SPKF-a i tezinski faktori
Nx  = 1;                            % dimenzija stanja (SoC)
Nxa = 3;                            % augmentirano stanje [x, w, v]
Nz  = 1;                            % dimenzija mjerenja (napon)

h = sqrt(3);
Wmx(1) = (h*h - Nxa)/(h*h);
Wmx(2) = 1/(2*h*h);
Wcx = Wmx;
Wmxz = [Wmx(1), repmat(Wmx(2), [1, 2*Nxa])]';

%% 3. Sumovi i inicijalizacija filtara
SigmaW = 1e-6;                      % kovarijanca procesnog suma
SigmaV = 0.0004;                    % kovarijanca suma voltmetra (~20 mV)
SigmaX = 0.05;                      % pocetna nesigurnost procjene

xtrue = 0.85;                       % stvarni pocetni SoC (85%)
xhat  = 0.50;                       % pocetna procjena (50% - test konvergencije)

% Polja za spremanje podataka
xstore      = zeros(maxIter, 1);
xhatstore   = zeros(maxIter, 1);
SigmaXstore = zeros(maxIter, 1);
Istore      = zeros(maxIter, 1);
Vstore      = zeros(maxIter, 1);

%% 4. Petlja estimacije
for k = 1:maxIter
    % Profil struje opterecenja (gradska voznja, semafor, gas, regeneracija)
    if k < 50
        Ik = 2.0 + 1.0*sin(k/5);
    elseif k < 100
        Ik = 0.0;
    elseif k < 200
        Ik = 5.0 + 2.0*randn;
    elseif k < 250
        Ik = -2.5;                  % punjenje
    else
        Ik = 1.5;
    end
    
    % SPKF Step 1a: State estimate time update
    xhata = [xhat; 0; 0];
    Pxa = blkdiag(SigmaX, SigmaW, SigmaV);
    sPxa = chol(Pxa, 'lower');
    
    % Generiranje 7 sigma tocaka
    X = xhata(:, ones([1, 2*Nxa+1])) + h*[zeros([Nxa, 1]), sPxa, -sPxa];
    
    % Propagacija kroz fiziku (Coulomb counting)
    Xx = X(1,:) - (dt / Q_total)*Ik + X(2,:);
    xhat = Xx * Wmxz;
    
    % SPKF Step 1b: Covariance of prediction 
    Xs  = (Xx(:, 2:end) - xhat*ones([1, 2*Nxa])) * sqrt(Wcx(2));
    Xs1 = Xx(:, 1) - xhat;
    SigmaX = Xs*Xs' + Wcx(1)*(Xs1*Xs1');
    
    % Stvarna celija (simulacija)
    w = sqrt(SigmaW) * randn;
    v = sqrt(SigmaV) * randn;
    xtrue = xtrue - (dt / Q_total)*Ik + w;
    ztrue = ocv_func(xtrue) - R0*Ik + v;
    
    % SPKF Step 1c: Output estimate 
    Z = ocv_func(Xx) - R0*Ik + X(3,:);
    zhat = Z * Wmxz;
    
    % SPKF Step 2a: Estimator gain matrix 
    Zs  = (Z(:, 2:end) - zhat*ones([1, 2*Nxa])) * sqrt(Wcx(2));
    Zs1 = Z(:, 1) - zhat;
    SigmaXZ = Xs*Zs' + Wcx(1)*(Xs1*Zs1');
    SigmaZ  = Zs*Zs' + Wcx(1)*(Zs1*Zs1');
    Lx = SigmaXZ / SigmaZ;
    
    % SPKF Step 2b: Measurement state update
    xhat = xhat + Lx*(ztrue - zhat);
    xhat = max(0, min(1, xhat));     % zasicenje na fizikalne granice [0, 1]
    
    % SPKF Step 2c: Measurement covariance update 
    SigmaX = SigmaX - Lx*SigmaZ*Lx';
    
    % Spremanje rezultata
    xstore(k)      = xtrue;
    xhatstore(k)   = xhat;
    SigmaXstore(k) = SigmaX;
    Istore(k)      = Ik;
    Vstore(k)      = ztrue;
end

%% 5. Proracun pogreske estimacije
eval_idx = 20:maxIter;              % evaluacija nakon prijelazne pojave
rmse = sqrt(mean((xstore(eval_idx) - xhatstore(eval_idx)).^2)) * 100;
max_err = max(abs(xstore(eval_idx) - xhatstore(eval_idx))) * 100;

fprintf('RMSE u stacionarnom stanju: %.3f %%\n', rmse);
fprintf('Maksimalna greska:          %.3f %%\n', max_err);

%% 6. Crtanje grafova
figure(1); clf;

subplot(3, 1, 1);
plot(1:maxIter, Istore, 'r', 'LineWidth', 1.2); grid on;
ylabel('Struja [A]');
title('Profil struje opterecenja');

subplot(3, 1, 2);
plot(1:maxIter, Vstore, 'm', 'LineWidth', 1); grid on;
ylabel('Napon [V]');
title('Izmjereni napon na klemama (sa sumom senzora)');

subplot(3, 1, 3);
plot(1:maxIter, xstore*100, 'k-', 'LineWidth', 2); hold on;
plot(1:maxIter, xhatstore*100, 'b--', 'LineWidth', 1.8);
plot(1:maxIter, (xhatstore + 3*sqrt(SigmaXstore))*100, 'g:', 'LineWidth', 1.2);
plot(1:maxIter, (xhatstore - 3*sqrt(SigmaXstore))*100, 'g:', 'LineWidth', 1.2);
grid on;
xlabel('Vrijeme [s]'); ylabel('SoC [%]');
legend('Stvarni SoC', 'SPKF procjena', '\pm 3\sigma granice', 'Location', 'best');
title('Rezultat estimacije stanja napunjenosti');