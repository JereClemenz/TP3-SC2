clc, clear all, close all;

%Motor en tiempo discreto, mas observador e integrador

%Curvas compartidas------------------------------------------------------
ruta_archivo = 'D:\Facultad Jeremías\2023 primer semestre\Sistemas de control 2\TPS\TP3 pucheta\Motor\Curvas_Medidas_Motor_2023.xls';
hoja_trabajo = 'Hoja1';
[num, txt, raw] = xlsread(ruta_archivo, hoja_trabajo);
tiempo = num(:, 1); % Primera columna de datos numéricos
angulo = num(:, 2); % Segunda columna de datos numéricos
velocidadAngular = num(:, 3);
CorrienteDeArm = num(:, 4);
TensionAplicada = num(:, 5);
Torque = num(:, 6);
%-----------------------------------------------------------------------

%Parametros de motor
Laa=0.56e-3;
J=0.0019;
Ra=1.35;
Bm=0.000792; %Pagina 357 Benjamin Kuo
Ki=.1;
Km=.1;

Ts=1e-4; %Tiempo de muestreo

A=[-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B=[1/Laa; 0; 0];
C=[0 0 1; 0 1 0];
D=[0];

sys1=ss(A,B,C,D);
sys_d=c2d(sys1,Ts,'zoh');

A=sys_d.a;
B=sys_d.b;
C=sys_d.c;


%Agrego un integrador
Aamp=[A,zeros(3,1);-C(1,:)*A, 1];
Bamp=[B;-C(1,:)*B];
Camp=[C(1,:) 0];

%Verifico controlabilidad
rank([Bamp Aamp*Bamp Aamp^2*Bamp Aamp^3*Bamp Aamp^4*Bamp]);%--> 4

%Controlador LQR
Q=1*diag([1 0.001 1 10]);    R=10e2;
Q=1*diag([1 1 1 0.01]);    R=0.5;
Kamp=dlqr(Aamp,Bamp,Q,R);
K=Kamp(1:3);
Ki=-Kamp(4);

%Observador----------------------------------------------------------------

Co=B';
Ao=A';
Bo=C';

Qo=1*diag([1 1 1]);    Ro=1e4;
Ko=dlqr(Ao,Bo,Qo,Ro);

%Simulación del control:
Tsim=10;
T_switch=5;
h=1e-5;
Kmax=Tsim/Ts;
pasos=round(Tsim/h);
t=0:h:(Tsim);

%Entradas del sistema
flag=1;
contador=0;
ref=zeros(1,round(Tsim/h+1));
for i=1:(Tsim/h+1)
    if flag==1
        ref(1,i)=pi/2;
        contador=contador+1;
        if contador==round(5/h)
            flag=0;
            contador=0;
        end
    else
        ref(1,i)=-pi/2;
        contador=contador+1;
        if contador==round(5/h)
            flag=1;
            contador=0;
        end
    end
    
end
% figure(2)
% plot(t,ref);
% title('Referencia');
flag=1;
contador=0;
tLin=zeros(1,round(Tsim/h+1));
for i=1:(Tsim/h+1)
    if flag==1
        %     tLin(1,i)=0.103;
        tLin(1,i)=0;
        contador=contador+1;
        if contador==round(5/h)
            flag=0;
            contador=0;
        end
    else
%         tLin(1,i)=0.103;
        tLin(1,i)=1.75;
%         tLin(1,i)=1.15;
        contador=contador+1;
        if contador==round(5/h)
            flag=1;
            contador=0;
        end
    end
    
end
% figure(3)
% plot(t,tLin);
% title('Torque de entrada');

x(1,1)=0;
x(2,1)=0;
x(3,1)=0;
x(4,1)=0;

xobs(1,1)=0;
xobs(2,1)=0;
xobs(3,1)=0;

x_ts=x((1:3),1);
v_ts=x(4,1);
z=1;

for i=1:1:Kmax
    x_k=x_ts;
    v_k=v_ts;
%     u=-K(1:3)*x_k(1:3)+Ki*v_k; %Sin observador
    u=-K(1:3)*xobs(1:3)+Ki*v_k; %Con observador
    uant=u;
    
    %Alinealidad
    Alin=0;
    if abs(u)<Alin
        u=0;
    else
        u=sign(u)*(abs(u)-Alin);
    end
    
    ys=C*x(1:3,z); 
    for j=1:1:Ts/h 
        ua(z)=uant;
        x1p=-Ra*x(1,z)/Laa-Km*x(2,z)/Laa+u/Laa;
        x2p=Ki*x(1,z)/J-Bm*x(2,z)/J-tLin(z)/J;
        x3p=x(2,z);
        x_p=[x1p; x2p; x3p];
        
        x((1:3),z+1)=x((1:3),z)+h*x_p;
        z=z+1;
    end
    yhat=C*xobs;
    e=ys-yhat;
    xobs=A*xobs+B*u+Ko'*e;
    v_ts=v_ts+ref(z)-C(1,:)*x_ts;
    x_ts=x((1:3),z);
end


figure(1)
subplot(2, 2, 1);
hold on;
grid on;
plot(t(1:length(x(1,:))),x(1,:));
plot(tiempo,CorrienteDeArm,'g');
hold off
xlim([0 Tsim]);
% legend({'Corriente','Corriente no optimizada'});
legend({'Corriente'});
title('Corriente');
xlabel('Tiempo');

subplot(2, 2, 2);
hold on;
grid on;
plot(t(1:length(x(2,:))),x(2,:));
plot(tiempo,velocidadAngular,'g');
hold off
xlim([0 Tsim]);
% legend({'Velocidad angular','Velocidad angular no optmizada'});
legend({'Velocidad angular'});
title('Velocidad angular');
xlabel('Tiempo');
ylabel('rad/s');

subplot(2, 2, 3);
hold on;
grid on;
plot(t(1:length(x(3,:))),x(3,:));
% plot(tiempo,angulo,'g');
plot(t(1:length(ref)),ref,'k--');
hold off
xlim([0 Tsim]);
legend({'Salida','Posicion no optimizada','Referencia'});
legend({'Salida','Referencia'});
title('Posicion angular');
xlabel('Tiempo');
ylabel('Ángulo');


subplot(2, 2, 4);
zona_muerta=ones(2,length(ua));
zona_muerta(1,:)=zona_muerta(1,:)*Alin;
zona_muerta(2,:)=zona_muerta(2,:)*-Alin;
hold on;
plot(t(1:length(ua)),ua);
plot(tiempo,TensionAplicada,'g');
plot(t(1:length(ua)),zona_muerta,'k--');
hold off
xlim([0 Tsim]);
grid on;
title('Acción de control');
xlabel('Tiempo');
ylabel('Voltaje');
% legend({'Acción de control','Accion de control no optimizada','Zona muerta'});
legend({'Acción de control','Zona muerta'});

figure (2)
hold on;
grid on;
plot(x(3,:),x(2,:));
hold off
xlabel('Pos. angular');
ylabel('Vel. angular');
title('Plano de fase');