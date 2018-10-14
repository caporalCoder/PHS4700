%
% Definir donnees pour les cinq tests
%
rbi([1:3],1)=[0.2;89.8;0.11];
vbi([1:3],1)=[5.3;-21;16.5];
wbi([1:3],1)=[0;0;6.3]; 

rbi([1:3],2)=[0.2;89.8;0.11];
vbi([1:3],2)=[5.3;-21;16.5];
wbi([1:3],2)=[0;-5;-6.3];

rbi([1:3],3)=[0.2;89.8;0.11];
vbi([1:3],3)=[5.3;-21;16.5];
wbi([1:3],3)=[0;0;-6.5];

rbi([1:3],4)=[0.2;89.8;0.11];
vbi([1:3],4)=[5.3;-21;18.3];
wbi([1:3],4)=[0;0;-6.3];

rbi([1:3],5)=[120-0.2;89.8;0.11];
vbi([1:3],5)=[-5.3;-21;16.5];
wbi([1:3],5)=[0;-3;6.3]; 
%
% Rouler les cinq tests
%
for test=1:5
  fprintf('\nSimulation %3d\n',test);
  fprintf('rbi (m)  = (%8.4f,%8.4f,%8.4f)  \n',rbi(1,test),rbi(2,test),rbi(3,test));
  fprintf('vbi (m/s)= (%8.4f,%8.4f,%8.4f)  \n',vbi(1,test),vbi(2,test),vbi(3,test));
  fprintf('wbi (r/s)= (%8.4f,%8.4f,%8.4f)  \n',wbi(1,test),wbi(2,test),wbi(3,test));
  [But tf rbf vbf]=Devoir2(rbi(:,test),vbi(:,test),wbi(:,test));
  if But == 1
    fprintf('But ! \n');
  elseif But == 0
    fprintf('Le ballon touche le sol en premier\n');
  elseif But == -1
    fprintf('Le ballon touche les montants du but en premier\n');
  elseif But == -2
    fprintf('Le ballon sort du terrain en premier\n');
  end
  fprintf('Vitesse finale de la balle     (%8.4f,%8.4f,%8.4f)  m/s\n',vbf(1),vbf(2),vbf(3));
  fprintf('Position finale de la balle    (%8.4f,%8.4f,%8.4f)  m \n',rbf(1),rbf(2),rbf(3));
  fprintf('\n\n');
end
