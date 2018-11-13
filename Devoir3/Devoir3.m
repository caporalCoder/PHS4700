%Devoir3


function [Coup tf vbaf vbof wbof rbaf rbof ]=Devoir3(vbal,wboi,tl)
    

    global mBoite
    global hBoite
    global RayonBoite
    global AxeCylindre
    global vInitialeBoite
    global rBoite
    global RayonBalle
    global mBalle
    global rBalle
    global k
    global aBalle
    global aBoite
    global coefficientRestitution
    global gravitation

    gravitation = 9.8;%m/s^2

    mBoite = 0.075; % kg
    hBoite = 0.15; %m
    RayonBoite = hBoite/sqrt(6); %m
    AxeCylindre = [0 0 1];
    vInitialeBoite = [0 0 0];
    rBoite = [3 0 10]; %m
    RayonBalle = 0.0335;
    mBalle = 0.058; %kg
    rBalle = [0 0 2]; %m
    k = 0.1; %kg/((m^2)s)
    aBalle = pi * (RayonBalle^2);
    aBoite = (RayonBoite^2) + (hBoite^2);
    coefficientRestitution = 0.5;
    precision_minimale = [[0.001 0.001 0.001] [Inf(1,3)] [Inf(1,3)];
    % pos_x, pos_y, pos_z v_x, v_y, v_z w_x w_y w_z tl
    q0Balle=[vbal(1) vbal(2) vbal(3) rBalle(1) rBalle(2) rBalle(3) 0 0 0 1 0 0 0 1 0 0 0 1];

    % Initialisation de l'etat initial de la boite.
    q0Boite=[vInitialeBoite(1) vInitialeBoite(2) vInitialeBoite(3) rBoite(1) rBoite(2) rBoite(3) wboi(1) wboi(2) wboi(3) 1 0 0 0 1 0 0 0 1];
    delta_t = 0.01;
    m=1;
    % Solution avec m=1
    t0 = 0;
    qsBoite=SEDRK4t0(q0Boite,t0,delta_t, 0);
    qsBalle=SEDRK4t0(q0Balle,t0,delta_t, 1);
    [Coup, normale] = FinSimulation(qsBalle(4:6)-qsBoite(4:6));
    t2=t0;
    while Coup < 0
       if t2 >= tl 
           qsBalle=SEDRK4t0(qsBalle,t2,delta_t, 1);
       end
       qsBoite=SEDRK4t0(qsBoite,t2,delta_t, 0);
       t2=t2+delta_t;
       [Coup, normale] = FinSimulation( qsBalle(4:6)-qsBoite(4:6));
    end
    [conv Err]=ErrSol(qsBoite,q0Boite,precision_minimale);
    [conv2 Err2]=ErrSol(qsBalle,q0Balle,precision_minimale);

    conv = and(conv, conv2);

    qs2Boite=qsBoite;
    qs2Balle=qsBalle;
    % Iteration avec m>1
    while not(conv)
        delta_t=delta_t/2;
        m=m+1;
        t2=t0;        
        qs2Boite=q0Boite;
        qs2Balle=q0Balle;
        Coup = -1;
        trajectoryBoite = [q0Boite(4:6)];
        trajectoryBalle = [q0Balle(4:6)];

        while Coup < 0
            qs2Boite=SEDRK4t0(qs2Boite,t2,delta_t, 0);
            if t2 >= tl 
                qs2Balle=SEDRK4t0(qs2Balle,t2,delta_t, 1);
            end
            t2=t2+delta_t;
            [Coup, normale] = FinSimulation(qs2Balle(4:6)-qs2Boite(4:6));
            trajectoryBoite = [qs2Boite(4:6)];
            trajectoryBalle = [qs2Balle(4:6)];
        end

        
        [conv Err]=ErrSol(qs2Boite,qsBoite,precision_minimale);
        [conv2 Err2]=ErrSol(qs2Balle,qsBalle,precision_minimale);    
        conv = and(conv, conv2);
        
        qsBoite=qs2Boite;
        qsBalle=qs2Balle;
        if m>10
            break;
        end
    end
    qsBoite=qs2Boite+Err/15;
    qsBalle=qs2Balle+Err2/15;
end

function F = ForceFortementVisqueuse(A, v)
    global k
    F = -k * A * v;
end

function qs=SEDRK4t0(q0,t0,DeltaT, Boite_ou_Balle) %0 pour boite et 1 pour balle
    % Solution equations differentielles par methode de RK4
    % Equation a resoudre : dq/dt=g(q,t)
    % avec
    % qs : solution [q(to+DeltaT)]
    % q0 : conditions initiales [q(t0)]
    % DeltaT : intervalle de temps
    % g : membre de droite de ED.
    % C?est un m-file de matlab
    % qui retourne la valeur de g au temps choisi
    k1=g(q0,t0, Boite_ou_Balle);
    k2=g(q0+(k1*DeltaT/2),t0+DeltaT/2,Boite_ou_Balle);
    k3=g(q0+k2*DeltaT/2,t0+DeltaT/2, Boite_ou_Balle);
    k4=g(q0+k3*DeltaT,t0+DeltaT, Boite_ou_Balle);
    qs=q0+DeltaT*(k1+2*k2+2*k3+k4)/6;
end

function res=g(q0, t0, Boite_ou_Balle)
    global mBalle
    global mBoite
    % Rxx = rotx(t0*q0(7));
    % Ryy = roty(t0*q0(8));
    % Rzz = rotz(t0*q0(9));
    % MatRot = Rzz * Ryy * Rxx;
    if Boite_ou_Balle == 0
        acceleration= ForcesBoite(q0)/mBoite;
    else %Boite_ou_Balle == 1
        acceleration= ForcesBalle(q0)/mBalle;
    end
    res = [acceleration q0(1:3), [0,0,0], [q0(8)*q0(16)-q0(9)*q0(13),q0(8)*q0(17)-q0(9)*q0(14),q0(8)*q0(18)-q0(9)*q0(15)], [q0(9)*q0(10)-q0(7)*q0(16),q0(9)*q0(11)- q0(7)*q0(17),q0(9)*q0(12)- q0(7)*q0(18)], [q0(7)*q0(13)-q0(8)*q0(10),q0(7)*q0(14)- q0(8)*q0(11),q0(7)*q0(15)- q0(8)*q0(12)] ];
end

function [conv, Err]=ErrSol(qs1,qs0,epsilon)
    % Verification si solution convergee
    %   conv      : variable logique pour convergence
    %               Err<epsilon pour chaque elements
    %   Err       : Difference entre qs1 et qs0 
    %   qs1       : nouvelle solution
    %   qs0       : ancienne solution
    %   epsilon   : pr?cision pour chaque variable
    Err = qs1-qs0;
    conv = all(abs(Err) < epsilon);
end

function [vBoiteF, vBalleF] = vitesseApresCollision(n, vBoite, vBalle)
    global coefficientRestitution
    global mBalle
    global mBoite
    %n dirige vers l'interieur de la balle
    vitesseRelativeInitiale = dot(n, vBoite - vBalle);
    j = -((1+coefficientRestitution)*vitesseRelativeInitiale)/((1/mBalle)+(1/mBoite));
    vBalleF = vBalle + (j* n/mBalle);
    vBoiteF = vBoite - (j* n/mBoite);

    
end

function res=ForcesBoite(q0)
    global mBoite
    global gravitation
    global aBoite % Aire Boite
        
    Fg = [0  0 -mBoite*gravitation];
    F_vis = ForceFortementVisqueuse(aBoite, q0(1:3));

    res = Fg + F_vis ;
end

function res=ForcesBalle(q0)
    global mBalle
    global gravitation
    global aBalle % Aire Balle
    
    Fg = [0, 0, -mBalle * gravitation];
    F_vis = ForceFortementVisqueuse(aBalle, q0(1:3));

    res = Fg + F_vis ;
end

function [Coup, normale] = FinSimulation(diffCentreMasse)
    global RayonBalle
	global RayonBoite
    global hBoite
    %Surface
	condition1 = (-hBoite/2 < diffCentreMasse(3)) && (diffCentreMasse(3) < hBoite/2) && (sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2) <= (RayonBalle + RayonBoite));
	% Face Inferieur
	condition2 = (sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2) < RayonBoite) && ((hBoite/2 - RayonBalle) <= diffCentreMasse(3)) && ((-hBoite/2 + RayonBalle) >= diffCentreMasse(3));
	% Face Superieur
	condition3 = (sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2) < RayonBoite) && (diffCentreMasse(3) <= (hBoite/2 + RayonBalle)) && (diffCentreMasse(3) >= (hBoite/2 - RayonBalle));
	% Arrete Inferieur
	condition4 = ((hBoite/2 - RayonBalle) <= diffCentreMasse(3) && diffCentreMasse(3) <= -hBoite/2) && RayonBoite <= sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2) && sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2) <= (RayonBoite + sqrt(RayonBalle^2 - (abs(diffCentreMasse(3)) - hBoite/2)^2));
	% Arrete Superieur
	condition5 = (hBoite/2 <= diffCentreMasse(3) && diffCentreMasse(3) <= (hBoite/2 + RayonBalle)) && (RayonBalle <= sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2)) && sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2) <= (RayonBoite + sqrt(RayonBalle^2 - (abs(diffCentreMasse(3)) - hBoite/2)^2));

	% sol
	condition6 = (diffCentreMasse <= RayonBalle);	
	Coup = -1;
	normale = [0 0 0];

	if (condition6)
		Coup = 0;
		normale = [0 0 0];
	elseif (condition1)
		%
		Coup = 1;
		k = RayonBoite / (sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2));
		normale = [k*diffCentreMasse(1) k*diffCentreMasse(2) diffCentreMasse(3)];
	elseif (condition2)
		%
		Coup = 1;
		k = RayonBoite / (sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2));
		normale = [diffCentreMasse(1) diffCentreMasse(2) -hBoite/2];
	elseif (condition3)
		%
		Coup = 1;
		k = RayonBoite / (sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2));
		normale = [diffCentreMasse(1) diffCentreMasse(2) hBoite/2];
	elseif (condition4)
		%
		Coup = 1;
		k = RayonBoite / (sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2));
		normale = [k * diffCentreMasse(1) k * diffCentreMasse -hBoite/2];
	elseif (condition5)
		%
		Coup = 1;
		k = RayonBoite / (sqrt(diffCentreMasse(1)^2 + diffCentreMasse(2)^2));
		normale = [k * diffCentreMasse(1) k * diffCentreMasse hBoite/2];
	end

	normale = transpose(normale);

end