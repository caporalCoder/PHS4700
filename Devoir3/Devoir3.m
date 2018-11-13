%Devoir3


function [Coup tf vbaf vbof wbof rbaf rbof ]=Devoir3(vbal,wboi,tl)
    

global mBoite
global hBoite

global RayonBoite

global AxeCylindre

global vInitialeBoite

global rCM

global RayonBalle

global mBalle

global posBalleDepart

global k
global aBalle
global aBoite

global coefficientRestitution


mBoite = 0.075; % kg
hBoite = 0.15; %m

RayonBoite = hBoite/sqrt(6);

AxeCylindre = [0 0 1];

vInitialeBoite = 0;

rCM = [3 0 10]; %m
%w = constante;

RayonBalle = 0.0335;

mBalle = 0.058; %kg

posBalleDepart = [0 0 2]; %m

k = 0.1; %kg/((m^2)s)
aBalle = pi * (RayonBalle^2);
aBoite = (RayonBoite^2) + (hBoite^2);

coefficientRestitution = 0.5;

end
