%Jay Maini 101037537
%MNPA
set(0, 'DefaultFigureWindowStyle','docked')
clear all
close all
global G C F

%Circuit parameters
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
Ro = 1000;

L1 = 0.2;
a = 100;
C1 = 0.25;
G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
Go = 1/Ro;

% Note: There is an error within the G matrix causing incorrect results
% that I was not able to diagnose. See graphs.

G = [G1 -1 -G1 0 0 0 0;
    1 0 0 0 0 0 0;
    -G1 0 G1+G2 -1 0 0 0;
    0 0 1 -G3 0 0 0;
    0 0 0 0 G4 -1 -G4;
    0 0 0 0 -G4 0 G4-Go;
    0 0 0 a 0 -G4-Go 0];

C = [C1 0 -C1 0 0 0 0;
    0 0 0 0 0 0 0;
    -C1 0 C1 0 0 0 0;
    0 0 0 -L1 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];

F = [0;1;0;0;0;0;0];

% [V1 Iin V2 I3 V4 Icc Vo]

%DC Case

s = 0; %no capacitive/inductive effects
Vin = linspace(-10, 10, 100);
V3 = zeros(size(Vin));
Vo = zeros(size(Vin));
for i=1:length(Vin)
   V = (G+s.*C)\(F.*Vin(i)); 
   V3(i) = V(4)*G3;
   Vo(i) = V(7);
end


figure(1)
subplot(2,2,1)
plot(Vin,Vo)
hold on
plot(Vin,V3)
title('DC Case: Vin vs Vo & V3')
xlabel('Voltage')
ylabel('Gain (Vo/Vi)')
legend('Vo','V3')

f = linspace(0.01, 100, 10000);
w = 2*pi*f;
Vo = zeros(size(f));
for i = 1:length(f)
    s = j*w(i);
    V = (G+s*C)\F;  % Solving
    Vo(i) = V(7);
end
subplot(2,2,2)
plot(w,abs(Vo));
title('AC Sweep');
xlabel('\omega (rad/s)');
ylabel('Voltage');

subplot(2,2,3)
%normal distribution of Capacitance
C1_dist = C1 + 0.05 * randn(1000,1);
hist(C1_dist)

subplot(2,2,4)
Vo = zeros(size(C1_dist));
%Solved distribution of gain 
for i = 1:length(C1_dist)
    C1 = C1_dist(i);
    %recalculate C matrix
    C = [C1 0 -C1 0 0 0 0;
    0 0 0 0 0 0 0;
    -C1 0 C1 0 0 0 0;
    0 0 0 -L1 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];
    
    V = (G+s.*C)\F;
    Vo(i) = abs(V(7));
end
hist(Vo)


