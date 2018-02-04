clc;
clear;

T = 300;
m0=9.11*10^-31; %in kg
W=1*10^-7;
L=2*10^-7;
mn=0.26*m0;
kB=1.38*10^-23;
tmn=0.2*10^-12;
vth=sqrt(2*kB*T/mn);%thermal velocity
N=5; %number of particles
deltaT=5*10^-15;
tTotal=1000*deltaT;
mfp=tmn*vth; %mean free path

%Initial velocity 
vx=zeros(1,N);
vy=zeros(1,N);

%Initial position of each particles
xPos=zeros(1,N);
yPos=zeros(1,N);

%Initial delta position
deltaPx=zeros(1,N);
deltaPy=zeros(1,N);

%Initial angle 
theta=zeros(1,N);

for i = 1:N
    x=rand*L;
    y=rand*W;
    xPos(i)=xPos(i)+x;
    yPos(i)=yPos(i)+y;
    theta(i) = theta(i) + rand*2*pi;
    vx(i)=vth*cos(theta(i));
    vy(i)=vth*sin(theta(i));
    deltaPx(i)=deltaPx(i)+vx(i)*deltaT;
    deltaPy(i)=deltaPy(i)+vy(i)*deltaT;
end


for t = 0 : deltaT : tTotal    
     
    for i=1:N
        if yPos(i)+deltaPy(i)>W||yPos(i)+deltaPy(i)<0
            theta(i)=2*pi - theta(i);
            vy(i)=vth*sin(theta(i));
            deltaPy(i)=vy(i)*deltaT;
        end
    end
       
    xPos=xPos+deltaPx;
    %Periodic boundary condition in x direction
    Ix=xPos>L;
    xPos(Ix)=xPos(Ix)-L;
    Ix=xPos<0;
    xPos(Ix)=xPos(Ix)+L;
    
    
    yPos=yPos+deltaPy;
    %Boundary condition in y direction
    Iy=yPos>W;
    yPos(Iy)=yPos(Iy)-2*(yPos(Iy)-W);
    Iy=yPos<0;
    yPos(Iy)=-yPos(Iy);
    
    figure(1);
    plot(xPos,yPos,'.');
    hold on;
    xlim([0 L]);
    ylim([0 W]);
    
    KEsum=0;
    for i = 1:N
        KEsum = KEsum + (1/2)*mn*vth^2;
    end
    KEavg = KEsum /N; 
    T=KEavg/kB;
    
    figure(2);
    xlabel('Time(s)');
    ylabel('Temperature(K)');
    
    plot(t,T,'.r');
    xlim([0 tTotal]);
    hold on;
    pause(0.1)
    
end

hold off;

