
% Definition des bornes du terrain
terrain_min_x = 0;
terrain_max_x = 120;
terrain_min_y= 0;
terrain_max_y = 31.32;
terrain_sol = 0;

% Definition de la region du but
but_min_y= 41.34;
but_max_y = 48.36;
but_hauteur = 2.44;

gravitation = 9.8 %m/s^2

% Definition des informations du ballon
rayon_ballon = 0.11;
masse_ballon = 0.45;

% Definition des informations actuelles du ballon
aire_ballon = pi * rayon_ballon^2

%Propriete de l'environnement
Ro_air = 1.2754 %kg/m^3
Mu = 1.8 * 10^(−5) %kg(mxs)
function [But tf rf vf ] = Devoir2(ri,vi,wi)
    % Cette fonction permet de determiner si le but est compter,
    %       , le temps de la fin de la simulation, la position du 
    %       ballon a la fin de la simulation et de la vitesse
    %       angulaire.
    % ri : est la position du centre de masse du ballon
    % vi : est la vitesse du centre de masse du ballon
    % wi : est la vitesse angualire du centre de masse du ballon
    
    % Condition initiale de la simulation 
    % q0 = [ri_x ri_y ri_z, vi_x vi_y vi_z, w_x, w_y, w_z]
    q0 = [ri, vi, wi];
    
    % Intervale de temps de la simulation
    delta_t = 0.1;  %% a ajuster
    
    % Nombre d'iteration de la simulation
    nb_iterator = 10; %% a ajuster
    
    % Temps de demarrage de la simulation
    t0 = 0;
    
    
    % Realisation de la simulation 
    for n =  2 : nb_iteration + 1
        q0 = SEDRK4t0(q0,t0,delta_t)%,g)
        t0 = t0 + delta_t;
    end
    
    
    
end


function fin = FinSimulation(positionBallon)
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
    if (positionBallon(0) <= terrain_min_x || positionBallon(0) >= terrain_max_x) && ...
           (positionBallon(1) >= but_min_y && positionBallon(1) <= but_max_y)       
        fin = true;
        return;      
    end
    
    % Cas o� le ballon touche le sol
    if positionBallon(2) <= (terrain_sol + rayon_ballon)
        fin = true;
        return;
    end
    
    % Cas o� le ballon touche un des montants
    if (positionBallon(0) <= (terrain_min_x + rayon_ballon) || positionBallon(0) >= (terrain_max_x - rayon_ballon)) && ...
            (positionBallon(1) == but_min_y || positionBallon(1) == but_max_y ...
             || positionBallon(2) == but_hauteur)    
        fin = true;
        return;
    end
    
    % Cas o� le ballon sort du terain
    if positionBallon(0) < terrain_min_x || positionBallon(0) > terrain_max_x || ...
            positionBallon(1) < terrain_min_y || positionBallon(1) > terrain_max_y
        fin = true;
        return;
    end
    
    fin = false;
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
    k2=g(q0+k1*DeltaT/2,t0+DeltaT/2);
    k3=g(q0+k2*DeltaT/2,t0+DeltaT/2);
    k4=g(q0+k3*DeltaT,t0+DeltaT);
    qs=q0+DeltaT*(k1+2*k2+2*k3+k4)/6;
end

function res=g(q0, t0)
    global masse_ballon;
    res = [Forces(q0)/masse_ballon, q0(4:6)]
end

function res=Forces(q0)
    norm_v = norm(q0(4:6)) 
    norm_w = norm(q0(7:9))
    global masse_ballon
    %Poids
    global gravitation
    Fg = [0, 0, -masse_ballon * gravitation]

    %Force de frottement visqueux
    global aire_ballon
    global Ro_air
    global Mu
    C_vis = 0
    Re = Ro_air * norm_v * aire_ballon / Mu
    if Re < 100000
        C_vis = 0.235 * norm_v
    elseif Re < 135000 || Re > 100000
        C_vis = 0.235 * norm_v - 0.125 * norm_v *((Re - 100000)/35000)
    else
        C_vis = 0.110 * norm_v        
    end
    F_vis = -aire_ballon * Ro_air * C_vis * norm_v * q0(4:6)

    %Effet Magnus
    C_m = 0.1925 * ((norm_w * rayon_ballon)/(norm_v))^(1/4)
    cross_w_v = cross( q0(7:9) , q0(4:6) )
    F_m = Ro_air * C_m * aire_ballon  * (norm_v^2) * (cross_w_v/norm(cross_w_v))

    res = Fg + F_vis + F_m

end