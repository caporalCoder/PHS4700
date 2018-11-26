%%% Devoir 4

%% Function principal  du devoir 
function [tps, fTrain, Itrain] = Devoir4(vtrainkmh, fAvion)
    global c_son
    global v_avion
    global r_avion
    global i_son
    global r_train
    
    %%% Donnees sur l'avion
    v_avion = transpose(300 * [cos(pi/18) 0 sin(pi/18)]); %km/h
    r_avion = transpose([0 0 0]);
    i_son = 160; %dB

    %%% Donnees sur le train
    r_train = transpose([10 10 0]); %km

    %%% Donnees sur l'air et le son
    phi = 10; %centrigrade 
    taux_humidite = 0.7; 

    c_son = (331.3 + 0.606 * phi); % m/s
    
    tps = 0;
    currentIntensity = i_son;
    
    fTrain = [];
    Itrain = [];
    
    % Iterate while the intensity is greater than or egal to 20 dB
    while currentIntensity >= 20
        
        current_r_avion = r_avion + v_avion * tps;
        current_r_train = r_train + vtrainkmh * tps;
        % [delta_t1, pos2_after]
        [delta_t1, pos2_after] = computeDeltaTPosition(current_r_avion, current_r_train, vtrainkmh);
        %[delta_t2, pos1_after] = computeDeltaTPosition(current_r_train, current_r_avion, v_avion);
        
        u = pos2_after - current_r_avion;
        u = u / norm(u);
        v1 = (c_son - dot(vtrainkmh, u));
        currentFrequency = v1 / norm(v1) * fAvion;
        
        d = norm(u) - 100; %each 100m
        
        currentIntensity = 160 - 20 * log10(d /10) - A(v1);
        
        fTrain = [fTrain, currentFrequency];
        Itrain = [Itrain, currentIntensity];
        % Each time, add delta_t to the current time
        tps = tps + 1;
    end
    
end

%% Coeficient de viscosite
function coef_viscosite = A(v)
    coef_viscosite = 0.8 + 0.0041 * v; % dB/km
end

%% Calcul de la frequence et de l'intensite.
function [delta_t, pos] = computeDeltaTPosition(current_r_1, current_r_2, v_2)
    global c_son
    
    position_diff = current_r_2 - current_r_1;
    
    a1 = c_son^2 - norm(v_2)^2;
    b_t1 = dot(position_diff, v_2) ;
    c_t1 = norm(position_diff);
    
    delta_t = (b_t1 + sqrt(b_t1^2 + a1 * c_t1))/ a1;
    
    pos = current_r_2 + v_2 * delta_t;
    
end

