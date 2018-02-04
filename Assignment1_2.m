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
N=10; %number of particles
deltaT=5*10^-15;
tTotal=1000*deltaT;

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
    while x>0.4*L && x<0.6*L &&(y<0.4*W||y>0.6*W)
        x=rand*L;
        y=rand*W;
    end
       
    xPos(i)=xPos(i)+x;
    yPos(i)=yPos(i)+y;
    theta(i) = theta(i) + rand*2*pi;
    sigma = sqrt(kB*T/mn);
    vx(i) = vth/sqrt(2)*randn;
    vy(i) =vth/sqrt(2)*randn;
    v(i) = sqrt(vx(i).*vx(i) + vy(i).*vy(i));
 
    deltaPx(i)=deltaPx(i)+vx(i)*deltaT;
    deltaPy(i)=deltaPy(i)+vy(i)*deltaT;    
end

%Part2 question a:histogram of initial velocity for each particle
figure(1); 
xlabel('v(m/s)');
vAvg = mean(v);
hist(v);
xp=xPos;
yp=yPos;
vyp=vy;
vxp=vx;
xPosp=xPos;
yPosp=yPos;
dt=zeros(1,N);
sumP=zeros(1,N);
sumt=zeros(1,N);
sumNP=zeros(1,N);
sumNt=zeros(1,N);

for t = 0 : deltaT : tTotal  

    dt=dt+deltaT;
    for i=1:N
        P = 1-exp(-deltaT/tmn);
        if P > rand()
            vx(i) = vth/sqrt(2)*randn;
            vy(i) =vth/sqrt(2)*randn;
            deltaPx(i)=deltaPx(i)+vx(i)*deltaT;
            deltaPy(i)=deltaPy(i)+vy(i)*deltaT;
        end
    end
    
    for i=1:N
        if yPos(i)+deltaPy(i)>W||yPos(i)+deltaPy(i)<0
            vy(i)=-vy(i);
            deltaPy(i)=vy(i)*deltaT;
        end
        if xPos(i)>0.4*L && xPos(i)<0.6*L
            if yPos(i)+deltaPy(i)<0.4*W||yPos(i)+deltaPy(i)>0.6*W
                vy(i)=-vy(i);
                deltaPy(i)=vy(i)*deltaT;
            end
        end
        if yPos(i)<0.4*W||yPos(i)>0.6*W
            if xPos(i)+deltaPx(i)>0.4*L && xPos(i)+deltaPx(i)<0.6*L
                vx(i)=-vx(i);
                deltaPx(i)=vx(i)*deltaT;
            end
        end
        
    end
       
    xPos=xPos+deltaPx;
    yPos=yPos+deltaPy;
    
    %sum of the free path for each particle
    for i=1:N
        if vx(i)~=vxp(i) || vy(i)~=vyp(i)
            FP(i)=sqrt((xPos(i)-xPosp(i))^2 + (yPos(i)-yPosp(i))^2);
            tc(i)=dt(i);
            sumP(i)=sumP(i)+FP(i);
            sumNP(i)=sumNP(i)+1;
            sumt(i)=sumt(i)+tc(i);
            sumNt(i)=sumNt(i)+1;
            dt(i)=0;
            xPosp(i)=xPos(i);
            yPosp(i)=yPos(i);
        end
    end
    
    %Periodic boundary condition in x direction
    Ix=xPos>L;
    xPos(Ix)=xPos(Ix)-L;
    Ix=xPos<0;
    xPos(Ix)=xPos(Ix)+L;
    
%     Iy=(yPos>0.6*W | yPos<0.4*W);
%     Ix=(xPos>0.4*L & xPos<0.6*L);
%     Ixy=Iy & Ix;
%     
%     Ixx = xPos(Ixy)>xp(Ixy);
%     xPos(Ixx)=xPos(Ixx)-2*(xPos(Ixx)-0.4*L)
%     Ixx = xPos(Ixy)<xp(Ixy);        
%     xPos(Ixx)=xPos(Ixx)+2*(0.6*L-xPos(Ixx));
%     Iyy = yPos(Ixy)>yp(Ixy);
%     yPos(Iyy)=yPos(Iyy)-2*(yPos(Iyy)-0.6*W)
%     Iyy = yPos(Ixy)<yp(Ixy);        
%     yPos(Iyy)=yPos(Iyy)+2*(0.4*W-yPos(Iyy));
%     
%    
%     
%     %Boundary condition in y direction
%     Iy=yPos>W;
%     yPos(Iy)=yPos(Iy)-2*(yPos(Iy)-W);
%     Iy=yPos<0;
%     yPos(Iy)=-yPos(Iy);
    
    %Part3 question a:2-D plot of particle trajectories
    figure(2);
    plot(xPos,yPos,'.');
    hold on;
    xlim([0 L]);
    ylim([0 W]);
    
    %block definition
    line([0.8*10^-7 0.8*10^-7],[W 0.6*W]);
    line([1.2*10^-7 1.2*10^-7],[W 0.6*W]);
    line([0.8*10^-7 0.8*10^-7],[0 0.4*W]);
    line([1.2*10^-7 1.2*10^-7],[0 0.4*W]);
    line([0.8*10^-7 1.2*10^-7],[0.6*W 0.6*W]);
    line([0.8*10^-7 1.2*10^-7],[0.4*W 0.4*W]);
         
    KEsum=0;
    for i = 1:N
        v_Squared = vx(i)^2+vy(i)^2;
        KEsum = KEsum + (1/2)*mn*v_Squared;
    end
    KEavg = KEsum /N; 
    T=KEavg/kB;
    
    %Part 2, question c: Temperature plot
    figure(3);
    xlabel('Time(s)');
    ylabel('Temperature(K)');
    plot(t,T,'.r');
    xlim([0 tTotal]);
    hold on;
    pause(0.1)
    
    xp=xPos;
    yp=yPos;
    
end

%mean free path for each particle
mfp=zeros(1,N);
mtc=zeros(1,N);
for i=1:N
    mfp(i)=sumP(i)/sumNP(i);
    mtc(i)=sumt(i)/sumNt(i);
end

%Part 2, question d: MFP and tmn
%mean free path for the system
MeanFreePath = mean(mfp);
MeanCollisionTime = mean(mtc);
fprintf('Mean free path is %0.15f m\n',MeanFreePath);
fprintf('Mean collision time is %0.20f s\n',MeanCollisionTime);

%Part3 question c: electron density map 
P=zeros(200,100);
xPos = xPos.*10^9;
yPos = yPos.*10^9;
for i=1:N
    for j = 1:200
        for k = 1:100
            if xPos(i) > j && xPos(i)< (j+1) && yPos(i)>k && yPos(i)<(k+1)
                P(j,k) = P(j,k)+1;
            end
        end
    end
end
figure(4); 
surf(P);

%Part3 question d: temperature map 
Temp=zeros(20,10);
xPos = xPos.*0.1;
yPos = yPos.*0.1;
for i=1:N
    for j = 1:20
        for k = 1:10
            if xPos(i) > j && xPos(i)< (j+1) && yPos(i)>k && yPos(i)<(k+1)
                v_Squared = vx(i)^2+vy(i)^2;
                T = (1/2)*mn*v_Squared/kB;
                Temp(j,k) = Temp(j,k)+T; 
            end
        end
    end
    KE=0;
end
figure(5); 
surf(Temp);
           
hold off;

