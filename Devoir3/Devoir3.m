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
global posBalleDepart
global k
global aBalle
global aBoite
global coefficientRestitution


mBoite = 0.075; % kg
hBoite = 0.15; %m
RayonBoite = hBoite/sqrt(6); %m
AxeCylindre = [0 0 1];
vInitialeBoite = [0 0 0];
rBoite = [3 0 10]; %m
%w = constante;

RayonBalle = 0.0335;

mBalle = 0.058; %kg
posBalleDepart = [0 0 2]; %m

k = 0.1; %kg/((m^2)s)
aBalle = pi * (RayonBalle^2);
aBoite = (RayonBoite^2) + (hBoite^2);

coefficientRestitution = 0.5;

max_error = [transpose([0.001 0.001 0.001]) transpose(Inf(1,3)) transpose(Inf(1,3))];
% pos_x, pos_y, pos_z v_x, v_y, v_z w_x w_y w_z tl
q0Balle=[posBalleDepart(1) posBalleDepart(2) posBalleDepart(3) vbal(1) vbal(2) vbal(3) 0 0 0];

% Initialisation de l'etat initial de la boite.
q0Boite=[rBoite(1) rBoite(2) rBoite(3) vInitialeBoite(1) vInitialeBoite(2) vInitialeBoite(3) wboi(1) wboi(2) wboi(3)];

% Deplacer la boite jusqu'au moment de lancement de la balle
t = 0;
dT = 0.001; % A ajuster
while (t < tl)
    next_t = t + dT;
    q0Boite = SEDR4t0E(qBoite, t, next_t, max_error);
    %
    t = next_t;
end

%Solution
DeltaT = 0.001;
m=1;
 % Solution avec m=1
 qs1=SEDRK4t0(q0,t0,delta_t);
  
 But = FinSimulation(qs1(1:3));
 t2=t0;
 while But < -2
     qs1=SEDRK4t0(qs1,t2,delta_t);
     t2=t2+delta_t;
     But = FinSimulation(qs1(1:3));
 end
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
     end
     
     [conv Err]=ErrSol(qs2,qs1,precision_minimale);
     qs1=qs2;
     if m>10
         break;
     end
 end
 qs=qs2+Err/15;

end


function MI = MomentInertieSphere(m, r)
    Ic = 2 * m * r^2 / 3;
    MI = [Ic 0 0; 0 Ic 0; 0 0 Ic];
end

function MI = MomentInertieCylindre(m, r, l)
    Ixy = m * r^2 / 2 + m * l^2 / 12;
    Iz = m * r^2
    MI = [Ixy 0 0; 0 Ixy 0; 0 0 Iz];
end

function F = ForceFortementVisqueuse(k, A, v)
    F = -k * A * v;
end

function normal = VecteurNormale(p1, p2)
    normale = cross([(p2-p1),0],[0,0,1]);
    normale = normale(1:2);
end

function [qs, Err] = CalculTrajectoire(q0, t0, tf, epsilon)

end

function qs=SEDRK4t0(q0,t0,Deltat,g)
    % Solution equations differentielles par methode de RK4
    % Equation a resoudre : dq/dt=g(q,t)
    % avec
    %   qs        : solution [q(to+Deltat)]
    %   q0        : conditions initiales [q(t0)]
    %   Deltat    : intervalle de temps
    %   g         : membre de droite de ED. 
    %               C'est un m-file de matlab
    %               qui retourne la valeur de g au temps choisi
    k1=feval(g,q0,t0);
    k2=feval(g,q0+k1*Deltat/2,t0+Deltat/2);
    k3=feval(g,q0+k2*Deltat/2,t0+Deltat/2);
    k4=feval(g,q0+k3*Deltat,t0+Deltat);
    qs=q0+Deltat*(k1+2*k2+2*k3+k4)/6;
end

function res=g(q0, ~)%t0
    global masse_ballon;
    res = [transpose(q0(4:6)), transpose(Forces(q0)/masse_ballon), transpose([0 0 0])];
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
