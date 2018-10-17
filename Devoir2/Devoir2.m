
function [But tf rf vf ] = Devoir2(ri,vi,wi)
    % Cette fonction permet de determiner si le but est compter,
    %       , le temps de la fin de la simulation, la position du 
    %       ballon a la fin de la simulation et de la vitesse
    %       angulaire.
    % ri : est la position du centre de masse du ballon
    % vi : est la vitesse du centre de masse du ballon
    % wi : est la vitesse angualire du centre de masse du ballon
    global terrain_min_x
    global terrain_max_x
    global terrain_min_y
    global terrain_max_y
    global terrain_sol
    global but_min_y
    global but_max_y
    global but_hauteur   
    global rayon_ballon

    global masse_ballon
    global gravitation
    global aire_ballon
    global Ro_air
    global Mu
    % Definition des bornes du terrain
    terrain_min_x = 0;
    terrain_max_x = 120;
    terrain_min_y= 0;
    terrain_max_y = 90;
    terrain_sol = 0;

    % Definition de la region du but
    but_min_y= 41.34;
    but_max_y = 48.36;
    but_hauteur = 2.44;

    gravitation = 9.8; %m/s^2

    % Definition des informations du ballon
    rayon_ballon = 0.110;
    masse_ballon = 0.45;

    % Definition des informations actuelles du ballon
    aire_ballon = pi * rayon_ballon^2;

    %Propriete de l'environnement
    Ro_air = 1.2754; %kg/m^3
    Mu = 1.8 * 10^(-5); %kg(mxs)
    % Condition initiale de la simulation 
    % q0 = [ri_x ri_y ri_z, vi_x vi_y vi_z, w_x, w_y, w_z]
    q0 = [ri, vi, wi];
    
    % Intervale de temps de la simulation
    delta_t = 0.01;  %% a ajuster
    
    % Nombre d'iteration de la simulation
    nb_iteration = 10; %% a ajuster
    
    % Temps de demarrage de la simulation
    t0 = 0;
    
    precision_minimale = 0.001;
    q_tf_1 = SEDRK4t0(q0,t0,delta_t);
    n = 1;
    % Realisation de la simulation 
    while abs(norm(q_tf_1(1:3))-norm(q0(1:3))) > precision_minimale  %|| n =  2 : nb_iteration + 1
        n = n +  1;
        delta_t = delta_t/2;
        q0 = q_tf_1;
        q_tf_1 = SEDRK4t0(q0,t0,delta_t);%,g)
        t0 = t0 + delta_t;
        But = FinSimulation(q_tf_1(1:3));
        tf = t0;
        rf = q0(1:3);
        vf = q0(4:6);
        if(n > 10)            
            break;
        end
    end
    
    
    
end


function goal= FinSimulation(positionBallon)
    % Cette fontion permet de terminer le status d'une simulation
    % positionBallon : Represente la position du centre de masse
    %                  du ballon.
    % Cas o� le but est compte
    global terrain_min_x
    global terrain_max_x
    global terrain_min_y
    global terrain_max_y
    global terrain_sol
    global but_min_y
    global but_max_y
    global but_hauteur 
    global rayon_ballon
    
    % Cas o� le ballon entre dans le but
    if (positionBallon(1) <= terrain_min_x || positionBallon(1) >= terrain_max_x) && ...
           (positionBallon(2) >= but_min_y && positionBallon(2) <= but_max_y)   
        goal = 1;
        return;      
    end
    
    % Cas o� le ballon touche le sol
    if positionBallon(3) < (terrain_sol + rayon_ballon)
        goal = 0;
        return;
    end
    
    % Cas o� le ballon touche un des montants
    if (positionBallon(1) <= (terrain_min_x + rayon_ballon) || positionBallon(1) >= (terrain_max_x - rayon_ballon)) && ...
            (positionBallon(2) == but_min_y || positionBallon(2) == but_max_y ...
             || positionBallon(3) == but_hauteur)    
        goal = -1;
        return;
    end
    
    % Cas o� le ballon sort du terain
    if positionBallon(1) < terrain_min_x || positionBallon(1) > terrain_max_x || ...
            positionBallon(2) < terrain_min_y || positionBallon(2) > terrain_max_y
        goal = -2;
        return;
    end
    return;    
end


function qs=SEDRK4t0(q0,t0,DeltaT)%,g)
    % Solution equations differentielles par methode de RK4
    % Equation a resoudre : dq/dt=g(q,t)
    % avec
    % qs : solution [q(to+DeltaT)]
    % q0 : conditions initiales [q(t0)]
    % DeltaT : intervalle de temps
    % g : membre de droite de ED.
    % C?est un m-file de matlab
    % qui retourne la valeur de g au temps choisi
    k1=g(q0,t0);
    k2=g(q0+(k1*DeltaT/2),t0+DeltaT/2);
    k3=g(q0+k2*DeltaT/2,t0+DeltaT/2);
    k4=g(q0+k3*DeltaT,t0+DeltaT);
    qs=q0+DeltaT*(k1+2*k2+2*k3+k4)/6;
end

function res=g(q0, ~)%t0
    global masse_ballon;
    res = [transpose(Forces(q0)/masse_ballon), transpose(q0(4:6)), transpose([0 0 0])];
end

function res=Forces(q0)
    global rayon_ballon

    global masse_ballon
    global gravitation
    global aire_ballon
    global Ro_air
    global Mu
    norm_v = norm(q0(4:6)); 
    norm_w = norm(q0(7:9));
    
    %Poids
    
    Fg = [0, 0, -masse_ballon * gravitation];

    %Force de frottement visqueux
    Re = Ro_air * norm_v * aire_ballon / Mu;
    if Re < 100000
        C_vis = 0.235 * norm_v;
    elseif Re < 135000 || Re > 100000
        C_vis = 0.235 * norm_v - 0.125 * norm_v *((Re - 100000)/35000);
    else
        C_vis = 0.110 * norm_v;        
    end
    F_vis = -aire_ballon * Ro_air * C_vis * q0(4:6);

    %Effet Magnus
    C_m = 0.1925 * ((norm_w * rayon_ballon)/(norm_v))^(1/4);
    cross_w_v = cross( q0(7:9) , q0(4:6) );
    F_m = Ro_air * C_m * aire_ballon  * (norm_v^2) * (cross_w_v/norm(cross_w_v));

    res = Fg + F_vis + F_m;
end

