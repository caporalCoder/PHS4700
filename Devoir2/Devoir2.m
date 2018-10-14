
% Definition des bornes du terrain
terrain_min_x = 0;
terrain_max_x = 120;
terrain_min_y= 0;
terrain_max_y = 31.32;
terrain_sol = 0;

% Definition de la region du but
but_min_y= 12;
but_max_y = 19.32;
but_hauteur = 2.44;

% Definition des informations du ballon
rayon_ballon = 0.11;


% Definition des informations actuelles du ballon


function [But tf rf vf ] = Devoir2(ri,vi,wi)
    % Cette fonction permet de determiner si le but est compter,
    %       , le temps de la fin de la simulation, la position du 
    %       ballon a la fin de la simulation et de la vitesse
    %       angulaire.
    % ri : est la position du centre de masse du ballon
    % vi : est la vitesse du centre de masse du ballon
    % wi : est la vitesse angualire du centre de masse du ballon



end


function fin = FinSimulation(positionBallon)
    % Cette fontion permet de terminer le status d'une simulation
    % positionBallon : Represente la position du centre de masse
    %                  du ballon.
    % Cas où le but est compte
    
    global terrain_min_x
    global terrain_max_x
    global terrain_min_y
    global terrain_max_y
    global terrain_sol
    global but_min_y
    global but_max_y
    global but_hauteur   
    global rayon_ballon
    
    % Cas où le ballon entre dans le but
    if (positionBallon(0) == terrain_min_x || positionBallon(0) == terrain_max_x) && ...
           (positionBallon(1) > but_min_y && positionBallon(1) < but_max_y)       
        fin = true;
        return;      
    end
    
    % Cas où le ballon touche le sol
    if positionBallon(2) == (terrain_sol + rayon_ballon)
        fin = true;
        return;
    end
    
    % Cas où le ballon touche un des montants
    if (positionBallon(0) == (terrain_min_x + rayon_ballon) || positionBallon(0) == (terrain_max_x - rayon_ballon)) && ...
            (positionBallon(1) == but_min_y || positionBallon(1) == but_max_y ...
             || positionBallon(2) == but_hauteur)    
        fin = true;
        return;
    end
    
    % Cas où le ballon sort du terain
    if positionBallon(0) < terrain_min_x || positionBallon(0) > terrain_max_x || ...
            positionBallon(1) < terrain_min_y || positionBallon(1) > terrain_max_y
        fin = true;
        return;
    end
    
    fin = false;
    return;    
end


function qs = SEDRK4t0(q0,t0,DeltaT,g)
    % Solution equations differentielles par methode de RK4
    % Equation a resoudre : dq/dt=g(q,t)
    % avec
    % qs : solution [q(to+DeltaT)]
    % q0 : conditions initiales [q(t0)]
    % DeltaT : intervalle de temps
    % g : membre de droite de ED.
    % C?est un m-file de matlab
    % qui retourne la valeur de g au temps choisi
    k1=feval(g,q0,t0);
    k2=feval(g,q0+k1*DeltaT/2,t0+DeltaT/2);
    k3=feval(g,q0+k2*DeltaT/2,t0+DeltaT/2);
    k4=feval(g,q0+k3*DeltaT,t0+DeltaT);
    qs=q0+DeltaT*(k1+2*k2+2*k3+k4)/6;
end
