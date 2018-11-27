%%% Devoir 4

%% Function principal  du devoir 
function [tps fTrain Itrain] = Devoir4(vtrainkmh, fAvion)
    global c_son
    global v_avion
    global r_avion
    global i_son
    global r_train
    
    %fTrain = [];
   % Itrain = [];
    %%% Donnees sur l'avion
    v_avion = transpose(300 * [cos(pi/18) 0 sin(pi/18)]) / 3.6; %m/s
    r_avion = transpose([0 0 0]);
    i_son = 160; %dB

    %%% Donnees sur le train
    r_train = transpose([10 10 0]) * 1000; %m
    v_train= vtrainkmh / 3.6;
    %%% Donnees sur l'air et le son
    phi = 10; %centrigrade 
    taux_humidite = 0.7; 

    c_son = (331.3 + 0.606 * phi); % m/s
    
    t = 0;
    currentIntensity = i_son;
    
    % Calcul de l'intensite au au temps t= 0
    d = (r_train + vtrainkmh * t) - (r_avion + v_avion * t);
    u_initial = d / norm(d);
    
    currentIntensity = i_son - norm(d) / 1000 * A(fAvion);% - 20 * log10(norm(d) /10);
    delta_t = 1;
    t = t + delta_t;
    % Iterate while the intensity is greater than or egal to 20 dB
    while currentIntensity >= 20
        
        d = (r_train + vtrainkmh * t) - (r_avion + v_avion * t);
        u = d / norm(d);
        
        currentIntensity = i_son - norm(d) / 1000 * A(fAvion); % - 20 * log10(norm(d) /10);
        %currentIntensity
        if currentIntensity >= 20
            Itrain(t) = currentIntensity;
            currentFrequency = fAvion * (c_son - dot(v_train, u))/(c_son - dot(v_avion, u));
            fTrain(t) = currentFrequency;
        end
        % Each time, add delta_t to the current time
        t = t + delta_t;
    end
    %disp("end");
    q_initial = [transpose(r_train); transpose(r_avion)];
    d = pdist(q_initial, 'euclidean');
    tps = d / (c_son - dot(v_train, u_initial));

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

