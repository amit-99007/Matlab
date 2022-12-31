clear all;
close all;
clc
for run = 1:100                       %number of run
LB = [-10 -10];                           %parameter lower bound
UB = [10 10];                             %parameter upper bound
n= 100;                                  %number of particle
Kmax = 50;                              %maximum iterations
r=rand(n,2);                            %random number
pos0=[];

pos0(:,1)= LB(1,1) + (UB(1,1)-LB(1,1))*r(:,1);
pos0(:,2)= LB(1,2) + (UB(1,2)-LB(1,2))*r(:,2);
x=pos0;                                 %initial position                    
v = zeros(n,2);                         %initial velocity
%k=1;

for k=1:Kmax
for i=1:n   
    fv(i,1)= ((x(i,1)^2+x(i,2)^2)/50-cos(x(i,1)/(i)^(1/2))*cos(x(i,2)/(i)^(1/2))+1 );
    if k==1
        pbestval(i,1) = fv(i,1);
        pbest(i,:) = x(i,:);               %finding pbest
    elseif fv(i,1)< pbestval(i,1)          %checking condition if function value is less than pbest value
            pbestval(i,1)=fv(i,1);         %assign new pbest value as function value 
            pbest(i,:) = x(i,:);           %and correspondingly pbest will be at that position 
    end   
end

if k==1
    [gbestval(k,1),index] = min(fv);      %minimum f out and position 
    gbest(k,:) = pbest(index,:);          %finding gbest
elseif min(fv)< gbestval(k-1,1)           %if previous gbest value is greater than min of fv 
    [gbestval(k,1),index] = min(fv);      %updating gbestval to min of fv    
    gbest(k,:) = pbest(index,:);          %updating gbest to pbest
elseif min(fv) >= gbestval(k-1,1)         %if previous gbest value is less than min of fv
    gbest(k,:)=gbest(k-1,:);              %update gbest as previous one
    gbestval(k,1)=gbestval(k-1,1);        %update gbest value as previous gbest value
end
    
w=0.9 - (k/Kmax).*(0.9-0.4);              %calculation of w
c1=2.5 + (k/Kmax)*(0.5-2.5);              %calculation of c1
c2=0.5 + (k/Kmax)*(2.5-0.5);              %calculation of c2

%updating velocity and position
for i = 1:n
        v(i,1) = w*v(i,1) + c1*rand(1)*(pbest(i,1)-x(i,1)) + c2*rand(1)*(gbest(k,1)- x(i,1)); %updating velocity 1
        v(i,2) = w*v(i,1) + c1*rand(1)*(pbest(i,2)-x(i,2)) + c2*rand(1)*(gbest(k,2)- x(i,2)); %updating velocity 2
        x(i,:) = v(i,:) + x(i,:);         %updating positions
end
end
end
plot(1:Kmax,gbestval,'k*');                %plotting gbest value v/s k

