clc, clear all, close all;

% Pendulo con observador (posibilidad de usar sin observador). Bien implementado con dos observadores que nos permite mejorar las osculaciones. La unica diferencia con la v2 es que tiene los dos observadores y tiene una mejor dinamica, aunque son muy parecidos

m=.1;
Fricc=0.1; 
l=1.6;
g=9.8;
M=1.5;

Atc=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(m+M)/(l*M) 0];
Btc=[0; 1/M; 0; 1/(l*M)];
Ctc=[1 0 0 0; 0 0 1 0]; 
Dtc=[0; 0];


Ts=1e-2; %Tiempo de muestreo

%Sistema para masa 0.1----------------------------------------------------
%Calculo de sistema discreto y controlador LQR
sys1=ss(Atc,Btc,Ctc,Dtc);
sys1_d=c2d(sys1,Ts,'zoh');

A1=sys1_d.a;
B1=sys1_d.b;
C1=sys1_d.c; 

%Amplio el sistema
Aamp1=[A1,zeros(4,1);-C1(1,:)*A1, eye(1)];
Bamp1=[B1;-C1(1,:)*B1];

Q1=1*diag([1 1000 1000 1 .01]);    R1=1000;
% Q1=1*diag([.1 1000 1000 .1 .0005]);    R1=100;
Q1=1*diag([1 1000 1000 1 .0001]);    R1=1;

Kamp1=dlqr(Aamp1,Bamp1,Q1,R1);
PolosLazoCerrado1= eig(Aamp1-Bamp1*Kamp1) %polos de lazo cerrado

K1=Kamp1(1:4);
KI1=-Kamp1(5);

%%Sistema para masa 0.1*10-------------------------------------------------
%Calculo de sistema discreto y controlador LQR

%Segundo controlador para cuando cambie la masa:
m=.1*10;

A2_tc=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(m+M)/(l*M) 0];
B2_tc=[0; 1/M; 0; 1/(l*M)];
C2_tc=[1 0 0 0; 0 0 1 0]; 
D2_tc=[0; 0];

sys2=ss(A2_tc,B2_tc,C2_tc,D2_tc);
sys2_d=c2d(sys2,Ts,'zoh');

A2=sys2_d.a;
B2=sys2_d.b;
C2=sys2_d.c;

Aamp2=[A2,zeros(4,1);-C2(1,:)*A2, eye(1)];
Bamp2=[B2;-C2(1,:)*B2];

%Diseño con LQR del segundo controlador
Q2=1*diag([1 1000 1000 .1 .0001]);    R2=1000;%poquitas oscilaciones y buen tiempo de llegada

Kamp2=dlqr(Aamp2,Bamp2,Q2,R2);
PolosLazoCerrado2=eig(Aamp2-Bamp2*Kamp2) 

K2=Kamp2(1:4);
KI2=-Kamp2(5);

% %Observador------------------------------------------------------------
% Co=B1';
% Ao=A1';
% Bo=C1';
% 
% Qo=1e2*diag([.1 100 10000 .1]);    Ro=.1;
% Ko=dlqr(Ao,Bo,Qo,Ro);
% eig(Ao-Bo*Ko)
% %--------------------------------------------------------------------
%Observador 1------------------------------------------------------------
Co1=B1';
Ao1=A1';
Bo1=C1';

Qo1=1e2*diag([10 1 10000 1]);    Ro1=1;
Qo1=1e2*diag([10 1 10000 1]);    Ro1=0.01;

Ko1=dlqr(Ao1,Bo1,Qo1,Ro1);
polosObs1=eig(Ao1-Bo1*Ko1) 
%Observador 2------------------------------------------------------------
Co2=B2';
Ao2=A2';
Bo2=C2';

Qo2=1e2*diag([1 10 10000/2 .01]);    Ro2=10000;
% Qo2=1e2*diag([1 10 100 .01]);    Ro2=10;

Ko2=dlqr(Ao2,Bo2,Qo2,Ro2);
polosObs2=eig(Ao2-Bo2*Ko2)
%--------------------------------------------------------------------

%Simulacion
tsim=140;
h=1e-4;%paso

T_switch=tsim/2;
pasos=round(tsim/h);
t=0:h:(tsim);
Kmax=tsim/Ts;
m=.1;
ref=5*square(2*pi*t/(2*T_switch))+5;

masa=zeros(1,length(t));
for i=1:length(t)
    if i<=round(T_switch/h)
        masa(i)=0.1;
    elseif i>round(T_switch/h)
        masa(i)=1;
    end
end

%condiciones iniciales
delta(1)=0;        %x1
delta_p(1)=0;      %x2
theta(1)=pi;       %x3
theta_p(1)=0;      %x4
psi(1)=0;

x=zeros(5,pasos);
x(1,1)=delta(1);
x(2,1)=delta_p(1);
x(3,1)=theta(1);
x(4,1)=theta_p(1);
x(5,1)=psi(1);

xobs(1,1)=delta(1);
xobs(2,1)=delta_p(1);
xobs(3,1)=theta(1)*0;
xobs(4,1)=theta_p(1);

x_ts=x((1:4),1);
v_ts=x(5,1);
ua(1)=0;
z=1;
Xop=[0; 0; pi; 0];
phi_pp=0;

for i=1:1:Kmax
    x_ts=x((1:4),z);
    v_ts=v_ts+ref(z)-C1(1,:)*x_ts;
    ys=C1*(x(1:4,z));
    if masa(z)<0.5
        K=K1;
        Ki=KI1;
    else
        K=K2;
        Ki=KI2;
    end
%     u=-K(1:4)*(x_ts(1:4)-Xop)+Ki*v_ts; %Sin observador
    u=-K(1:4)*(xobs(1:4)-Xop)+Ki*v_ts; %Con observador
    uu=u;
    
    Alin=0.1;
    if abs(u)<Alin
        u=0;
    else
        u=sign(u)*(abs(u)-Alin);
    end
    
    for j=1:round(Ts/h)
        ua(z)=uu;%accion de control sobre la muestra
        %Evolucion del sistema en un Ts
        estados=x((1:4),z);
        delta_pp=(u-Fricc*estados(2)-masa(z)*l*phi_pp*cos(estados(3))+masa(z)*l*sin(estados(3))*estados(4)^2)/(M+masa(z));
        phi_pp=(g*sin(estados(3))-delta_pp*cos(estados(3)))/l;
        x1_p=estados(2);
        x2_p=delta_pp;
        x3_p=estados(4);
        x4_p=phi_pp;
        x_p_actual=[x1_p; x2_p; x3_p; x4_p];
        x((1:4),z+1)=x((1:4),z)+h*x_p_actual;
        z=z+1;
    end
    yobs=C1*(xobs);
    e=ys-yobs;
    
    
    if masa(z)<0.5
        xobs=A1*(xobs-Xop)+B1*u+Ko1'*e+Xop;
    else
        xobs=A2*(xobs-Xop)+B2*u+Ko2'*e+Xop;
    end
    
end

figure(1)
subplot(3,2,1)
hold on;
grid on;
plot(t(1:length(x(1,:))),x(1,:));
plot(t(1:length(ref)),ref,'k--');
xlim([0 tsim]);
title('Distancia');
xlabel('Tiempo');
legend({'Salida','Referencia'});

subplot(3,2,2)
hold on;
grid on;
plot(t(1:length(x(2,:))),x(2,:));
xlim([0 tsim]);
title('Velocidad');
xlabel('Tiempo');
legend({'Salida'});

subplot(3,2,3)
hold on;
grid on;
plot(t(1:length(x(3,:))),x(3,:)*(180/pi));
xlim([0 tsim]);
title('Angulo');
xlabel('Tiempo');

subplot(3,2,4)
hold on;
grid on;
plot(t(1:length(x(4,:))),x(4,:));
xlim([0 tsim]);
title('Velocidad angular');
xlabel('Tiempo');

subplot(3,1,3)
hold on;
grid on;
plot(t(1:length(ua)),ua); 
zonammuerta=ones(2,length(ua));
zonammuerta(1,:)=zonammuerta(1,:)*Alin;
zonammuerta(2,:)=zonammuerta(2,:)*-Alin;
plot(t(1:length(ua)),zonammuerta,'k');
xlim([0 tsim]);
title('Accion de control');
xlabel('Tiempo');
legend({'Accion de control','Zona muerta'});

figure (2)
plot(x(3,:)*(180/pi),x(4,:)*(180/pi));
title('Plano de fase');
grid on
xlabel('Angulo');
ylabel('Velocidad angular');