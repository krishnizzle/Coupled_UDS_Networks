function output = coupled_EIrate(varargin)
%This will make a coupled network. To get different amounts of SPA, change
%Jee2. To get different amounts of SPI, change JEtoE. 
p = inputParser; 
addParamValue(p,'tau_e',10)    %excitatory neuron membrane relaxation time (ms)
addParamValue(p,'tau_i',5)     %inhibitory neuron membrane relaxation time (ms)
addParamValue(p,'tau_a1',300)   %adaptation variable relaxation time (ms)
addParamValue(p,'tau_a2',300)   %adaptation variable relaxation time (ms)

addParamValue(p,'g_e',1*6)       %gain for firing rate of excitatory neurons (mV)
addParamValue(p,'g_i',5*6)       %gain for inhibitory (mV)

addParamValue(p,'g_e2',1*6)       %gain for firing rate of excitatory neurons (mV)
addParamValue(p,'g_i2',5*6)       %gain for inhibitory (mV)

%within network coupling
addParamValue(p,'Jee',5/6)       %synapse strength from excitatory to excitatory 
addParamValue(p,'Jii',0.5/6)     %synapse from inhibitory to inhibitory 
addParamValue(p,'Jie',10/6)      %synapse from excitatory to inhibitory
addParamValue(p,'Jei',1/6)       %synapse from inhibitory to excitatory
addParamValue(p,'Jee2',5/6)       %synapse strength from excitatory to excitatory 
addParamValue(p,'Jii2',0.5/6)     %synapse from inhibitory to inhibitory 
addParamValue(p,'Jie2',10/6)      %synapse from excitatory to inhibitory
addParamValue(p,'Jei2',1/6)       %synapse from inhibitory to excitatory

%from one network to another 
addParamValue(p,'JEtoE',0.7/6)        %synapse strength from E cells of Network 1 to Network 2

addParamValue(p,'beta1',1.1)      %adaptation strength
addParamValue(p,'beta2',1.1)
addParamValue(p,'theta_e1',6/(6*15))
addParamValue(p,'theta_i1',30/(6*15))
addParamValue(p,'theta_e2',6/(6*15))   %threshold for excitatory neurons
addParamValue(p,'theta_i2',30/(6*15))    %threshold for inhibitory neuorns 
addParamValue(p,'sig1',4.5)       %standard deviation of noise centered at zero
addParamValue(p,'sig2',3) 

addParamValue(p,'dt',0.2)       %time step for newton's method in ms
addParamValue(p,'dtk',0.8)      %time step for Runge Kutta in ms

addParamValue(p,'totaltime',100)     %total time in seconds 

parse(p,varargin{:});
p=p.Results;
dt = p.dt; dtk = p.dtk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Newton's Method 
%create time vector of 10 seconds
t = (1:dt:p.totaltime*1000);  %in ms
time = (1:10:p.totaltime*1000);
%initialize the rates and adaptation variables 
r_e1 = zeros(size(t)); r_e2 = r_e1;
r_i1 = zeros(size(t)); r_i2 = r_i1;
a_1 = zeros(size(t)); a_2 = a_1+2;

re1 = zeros(size(time)); re2 = re1; ri2= re1; ri1 = re1; a1 = re1; a2 = re2;
j=1;
%loop over time steps 
for i = 1:length(t)-1
    if i == round(length(t)/2)
        disp('HALFWAY')
    end
   %This is simple Newton's method 
   r_e1(i+1) = r_e1(i) + dt*dre1(p.g_e,p.theta_e1,p.tau_e,r_e1(i),p.Jee,p.Jei,r_i1(i),a_1(i),p.sig1);
   r_i1(i+1) = r_i1(i) + dt*dri(p.g_i,p.theta_i1,p.tau_i,r_i1(i),p.Jie,p.Jii,r_e1(i),p.sig1);
    a_1(i+1)  = a_1(i)   + dt*da(p.tau_a1,a_1(i),p.beta1,r_e1(i)); 
    
    r_e2(i+1) = r_e2(i) + dt*dre2(p.g_e2,p.theta_e2,p.tau_e,r_e2(i),p.Jee2,p.Jei2,r_i2(i),a_2(i),p.sig2,p.JEtoE,r_e1(i));
    r_i2(i+1) = r_i2(i) + dt*dri(p.g_i2,p.theta_i2,p.tau_i,r_i2(i),p.Jie2,p.Jii2,r_e2(i),p.sig2);
    a_2(i+1)  = a_2(i)   + dt*da(p.tau_a2,a_2(i),p.beta2,r_e2(i));
    
    if mod(i-1,10./dt)==0
        re1(j) = r_e1(i+1);
        re2(j) = r_e2(i+1);
        ri1(j) = r_i1(i+1);
        ri2(j) = r_i2(i+1);
        a1(j)  = a_1(i+1);
        a2(j)  = a_2 (i+1);
        j=j+1;
    end
end

output.r_e1 = re1;
output.r_e2 = re2;
output.r_i1 = ri1;
output.r_i2 = ri2;
output.a1 = a1;
output.a2 = a2;
output.time = time/1e3;

%

figure('units','normalized','outerposition',[0 0 1 1]);
s1 = subplot(2,1,1);
plot(output.time,output.r_i1,'b');
hold on;
plot(output.time,output.r_e1,'r');
hold on;
yyaxis right;
plot(output.time,output.a1,'k');
title({'Network 1:';strcat( 'EtoE:',num2str(p.JEtoE));strcat('beta2:',num2str(p.beta2));strcat('thetaE:',num2str(p.theta_e1))});

s2 = subplot(2,1,2);

plot(output.time,output.r_i2,'b');
hold on;
plot(output.time,output.r_e2,'r');
hold on;
yyaxis right;
plot(output.time,output.a2,'k');
title('Network 2');

linkaxes([s1, s2]);
%}
end


function dr_e1 = dre1(ge,theta_e,tau_e,re,jee,jei,ri,a,sig)

v = jee*re - jei*ri - a + normrnd(0,sig);
phi = ge*heaviside(v-theta_e)*(v-theta_e);
dr_e1 = 1/tau_e*(-re + phi);

end

function dr_e2 = dre2(ge,theta_e,tau_e,re2,jee,jei,ri,a,sig,JEE,re1)

v = jee*re2 - jei*ri - a + normrnd(0,sig) + JEE*re1;
phi = ge*heaviside(v-theta_e)*(v-theta_e);
dr_e2 = 1/tau_e*(-re2 + phi);

end

function dr_i = dri(gi,theta_i,tau_i,ri,jie,jii,re,sig)

v = jie*re - jii*ri + normrnd(0,sig);
phi = gi*heaviside(v-theta_i)*(v-theta_i);
dr_i = 1/tau_i*(-ri + phi);

end

function d_a = da(tau_a,a,beta,re)

d_a = 1/tau_a * (-a + beta*re);

end