
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
    delta_t = 0.1;  %% a ajuster
    
    % Nombre d'iteration de la simulation
    %nb_iteration = 1; %% a ajuster
    
    % Temps de demarrage de la simulation
    t0 = 0.0;
    
    precision_minimale = [transpose([0.001, 0.001, 0.001]) transpose(Inf(1,3)) transpose(Inf(1,3))];
    m=1;
    % Solution avec m=1
    qs1=SEDRK4t0(q0,t0,delta_t);
    But = FinSimulation(qs1(1:3));
    t2=t0;
    while But < -2
        qs1=SEDRK4t0(qs1,t2,delta_t);
        t2=t2+delta_t;
        But = FinSimulation(qs1(1:3));
    end;
    [conv Err]=ErrSol(qs1,q0,precision_minimale);
    qs2=qs1;
    % Iteration avec m>1
    while not(conv)
        delta_t=delta_t/2;
        m=m+1;
        t2=t0;        
        qs2=q0;
        But = -3;
        trajectory = [q0(1:3)];
        while But < -2
            qs2=SEDRK4t0(qs2,t2,delta_t);
            t2=t2+delta_t;
            But = FinSimulation(qs2(1:3));
            trajectory=[trajectory; qs2(1:3)];
        end;
        [conv Err]=ErrSol(qs2,qs1,precision_minimale);
        qs1=qs2;
        if m>10
            break;
        end;
    end;
    qs=qs2+Err/15;
    tf = t2;
    rf = qs(1:3);
    vf = qs(4:6);
    [x,y,z] = sphere();

    s = surf(x*rayon_ballon+ rf(1), y*rayon_ballon + rf(2), z*rayon_ballon + rf(3));
    s.EdgeColor = 'w';
    s.FaceColor = 'w';
    hold on;
    a1 = plot3(trajectory(:, 1),trajectory(:, 2),trajectory(: ,3), '-.r*');
    DessinerTerrain();
    axis('equal');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    legend(a1, 'Trajectoire Ballon');
    
    
    
end
function DessinerTerrain()
    global terrain_min_x
    global terrain_max_x
    %global terrain_min_y
    global terrain_max_y
    global terrain_sol
    global but_min_y
    global but_max_y
    global but_hauteur  
    patch(terrain_max_x*[0, 1, 1, 0], terrain_max_y*[0, 0, 1, 1], terrain_sol*[1, 1, 1, 1], 'green');
    patch(terrain_min_x* [1, 1, 1, 1], (but_max_y - but_min_y)*[0, 1, 1, 0] + but_min_y, but_hauteur*[0, 0, 1, 1], 'white');
    patch(terrain_max_x* [1, 1, 1, 1], (but_max_y - but_min_y)*[0, 1, 1, 0] + but_min_y, but_hauteur*[0, 0, 1, 1], 'white');
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
    goal = -3;
     % Cas o� le ballon touche un des montants
    if (intersectLineSphere([[0 but_min_y 0], [0 but_min_y but_hauteur]], positionBallon, rayon_ballon) || ...
        intersectLineSphere([[0 but_max_y 0], [0 but_max_y but_hauteur]], positionBallon, rayon_ballon) || ...
        intersectLineSphere([[0 but_min_y but_hauteur] [0 but_max_y but_hauteur]], positionBallon, rayon_ballon))
        goal = -1;
        return;
    end
     % Cas o� le ballon touche le sol
    if positionBallon(3) <= (terrain_sol + rayon_ballon)
        goal = 0;
        return;
    end
    % Cas o� le ballon entre dans le but
    if (positionBallon(1) <= terrain_min_x || positionBallon(1) >= terrain_max_x) && ...
           (positionBallon(2) >= but_min_y && positionBallon(2) <= but_max_y && positionBallon(3) <= but_hauteur)   
        goal = 1;
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
    res = [transpose(q0(4:6)), transpose(Forces(q0)/masse_ballon), transpose([0 0 0])];
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
    Re = Ro_air * norm_v * rayon_ballon / Mu;
    if Re < 100000
        C_vis = 0.235 * norm_v;
    elseif Re < 135000 && Re > 100000
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

function [conv Err]=ErrSol(qs1,qs0,epsilon)
    % Verification si solution convergee
    %   conv      : variable logique pour convergence
    %               Err<epsilon pour chaque elements
    %   Err       : Difference entre qs1 et qs0 
    %   qs1       : nouvelle solution
    %   qs0       : ancienne solution
    %   epsilon   : pr?cision pour chaque variable
    Err = qs1-qs0;
    conv = all(abs(Err) < epsilon);
    conv = conv(1);
end

function res=intersectLineSphere(line, sphereCoordinates, rayon)
    vecLine = line(4:6) - line(1:3);
    unitVecLine = vecLine / norm(vecLine);
    P1C = sphereCoordinates - line(1:3);           % Line from one point to center
    v = norm(cross(unitVecLine, P1C));
    res = (v <= rayon);
end





