%Written: March 2020
%Author: Amanda N. Laubmeier

%At time of use, this code was compatible with MATLAB R2014b, 
% version 8.4.0.150421, on an SMP Linux operating system.

%Results from the optimization (fmincon/multistart) may be impacted by 
%MATLAB version and OS. This code is made available without warranty.

%Queries about the code can be directed to anjlaubmeier (via gmail),
% but no guarantee of assistance is offered.

%This code can compute an optimal initial abundances of a theoretical
%   community of beetles and spiders with a common (aphid) prey, informed by 
%   empirical observations. The optimization problem is repeated under 
%   different assumptions for habitat use (diversity or overlap in predator 
%   foraging areas) and levels of intraguild interference or 
%   intraspecific competition. The performance of optimal communities is 
%   also evaluated under different temperatures.

%REFERENCES:
%MODEL, OPTIMIZATION CRITERIA: Laubmeier et. al., 2020. Ecosphere.
%ADDITIONAL DATA: Curtsdotter et. al., 2019. Journal of Animal Ecology.

%Must be able to read in data from doi.org/10.5061/dryad.41b5b06:
%   'BM FoodWeb Fulltime ___.csv' - species body masses
%   'Density numDate ___ Fulltime ___.csv' - predator abundances
%   'Density numDate Aphid Shifted Truncated ___.csv' - aphid abundances
%   'FoodWeb.csv' - binary interaction matrix
%   'Templog numDate comma ___.csv' - recorded field temperatures

function [opt1,opt2,opt3,div,tmp1,tmp3]=OP1results
close all, clear all,

%baseline parameter values
a0=.17; phi=1; h0=.12; x00=4*exp(19.75); 
v0=1; c0=0; b0=.4; 

Tspan=15:1:45;  %consider temperatures from 15 to 45 C
tBd=[0,30]; %times to consider: until day 30

opts = optimoptions('fmincon','Algorithm','interior-point');
%optimization options

name='results-'; %saving filename

%% Habitat Use Results

 RP=0; RB=0; PB=0; BB=0; PP=0; %assume complete overlap - no difference in foraging area
 [I,TV,W,y0,Ropt,rs,~]=OP1_initializingFun(RP,RB,PB,BB,PP); %initialize foodweb

[opt1,~]=OP1_resultsGen(Tspan,a0,phi,v0,h0,c0,b0,Ropt,x00,...
    rs,tBd,I,TV,W,y0,opts); %compute optimal predator community

RP=.4; RB=.4; PB=.8; BB=0; PP=0; %assume complete overlap within predator type
[I,TV,W,y0,Ropt,rs,~]=OP1_initializingFun(RP,RB,PB,BB,PP); %initialize foodweb
[opt2,~]=OP1_resultsGen(Tspan,a0,phi,v0,h0,c0,b0,Ropt,x00,...
    rs,tBd,I,TV,W,y0,opts); %compute optimal predator community

RP=.4; RB=.4; PB=.8; BB=.3; PP=.3; %assume every predator has its own foraging area
[I,TV,W,y0,Ropt,rs,~]=OP1_initializingFun(RP,RB,PB,BB,PP); %initialize foodweb
[opt3,~]=OP1_resultsGen(Tspan,a0,phi,v0,h0,c0,b0,Ropt,x00,...
    rs,tBd,I,TV,W,y0,opts); %compute optimal predator community

save([name,'habitat'],'opt1','opt2','opt3') 
%save in case you want to run other things without this


%% Interference Results

%this part is unnecessary if you haven't commented out above
RP=.4; RB=.4; PB=.8; BB=.3; PP=.3; %use low overlap / high foraging diversity case
[I,TV,W,y0,Ropt,rs,~]=OP1_initializingFun(RP,RB,PB,BB,PP); %initialize foodweb

cspan=0:.1:1; %c0 to consider
bspan=0:1:10; %b0 to consider - chosen to scale FR with c0
bs={'0','1','2','3','4','5','6','7','8','9','10'}; %readable names for c0
cs={'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'};% and b0

div=zeros(length(bspan),length(cspan)); %initialize for results

for j=1:length(bspan) %loop over b0 (interference)
    for k=1:length(cspan) %loop over c0 (competition)
        [relresies,~]=OP1_resultsGen(Tspan,a0,phi,v0,h0,cspan(k),bspan(j),Ropt,x00,...
            rs,tBd,I,TV,W,y0,opts); %compute optimal predator community
        save([name,'-b0-',bs{j},'-c0-',cs{k}],'relresies') %save in case you want raw results
          div(j,k)=OP1_diversityCounter(relresies); %quantify diversity
    end
end
save([name,'interference'],'bspan','bs','cspan','cs','div')
%save in case you want to run other things without this


%% Temperature Results
%THIS PORTION NEEDS opt1,opt3 from above - if running separately, make
%sure to load in those predator communities

day=3; %lenth of temperature chunks to save (3 days)
tf=30; %length of temperature trials to consider (30 days)
trial=5000;  %number of trials to repeat
tBd=[0,29]; %bounds to use for the averages (30 days starting at 0)

[lib,list]=OP1_tempgen(day); %reshape the temperature data
[c,e]=hist(list,14:46); %bin the temperatures
pb=[c(1)+c(2), c(3:end-2), c(end-1)+c(end)]; %store <15 and >45 together
pab=pb/sum(pb); %compute relative frequency of temperatures

% UNIFORM HABITAT - EVERYONE SHARES
RP=0; RB=0; PB=0; BB=0; PP=0;
[I,TV,W,~,Ropt,rs,~]=OP1_initializingFun(RP,RB,PB,BB,PP);

pKeep=[1,2,3,4,6,7]+1; %indices for preds in opt1
Ntest1=OP1_genPred(pKeep,opt1,pab); %opt1 abundances to test
cases=length(Ntest1(1,:)); %how many initial abundances I have (12)

J=zeros(trial,cases); %initialize
for r=1:trial %for all my replicate trials
    Tf=OP1_monthtemp(lib,day,tf); %build a temperature from my library
    for p=1:cases %then for all the predator abundances to test
        [Jt,~,~]=OP1_costTvar(I,TV,W,a0,phi,v0,h0,c0,b0,...
            Ropt,x00,rs,Ntest1(:,p),tBd,Tf); %compute the cost
        J(r,p)=Jt; %and store it
    end
end
tmp1=J;

%MIXED HABITAT
RP=.4; RB=.4; PB=.8; BB=.3; PP=.3;
[I,TV,W,~,Ropt,rs,~]=OP1_initializingFun(RP,RB,PB,BB,PP);

pKeep=[1,2,3,4,6,7,8]+1; %indices for preds in opt3
Ntest3=OP1_genPred(pKeep,opt3,pab); %opt3 abundances to test
cases=length(Ntest3(1,:)); %how many initial abundances I have (12)

J=zeros(trial,cases); %initialize
for r=1:trial %for all my replicate trials
    Tf=OP1_monthtemp(lib,day,tf); %build a temperature from my library
    for p=1:cases %then for all the predator abundances to test
        [Jt,~,~]=OP1_costTvar(I,TV,W,a0,phi,v0,h0,c0,b0,...
            Ropt,x00,rs,Ntest3(:,p),tBd,Tf); %compute the cost
        J(r,p)=Jt; %and store it
    end
end
tmp3=J;

save([name,'temperature'],'tmp1','tmp3','Ntest1','Ntest3')
%save in case you want to run other things without this


end

%% FUNCTIONS RELATED TO COST/OPTIMIZATION
function [relresies,Js]=OP1_resultsGen(Tspan,...
    a0,phi,v0,h0,c0,b0,Ropt,x00,rs,tBd,I,TV,W,y0,opts)
%function takes in temperatures (Tspan) and necessary inputs for
%   OP1_minField
%function outputs the optimal predator communities at those temperatures
%   (relresies) and associated costs (Js)
relresies=zeros(length(y0(2:end)),length(Tspan)); %initialize
Js=zeros(length(Tspan),1); %initialize
for k=1:length(Tspan) %loop over temperatures
    T=Tspan(k)+273.25; %convert to Kelvin
    [yOpt,J]=OP1_minField(a0,phi,v0,h0,c0,b0,Ropt,x00,T,...
        rs,tBd,I,TV,W,y0,opts); %find the optimal community at current temp
    relresies(:,k)=yOpt; %save results
    Js(k)=J; %save results
end
end

function [yOpt,J]=OP1_minField(a0,phi,v0,h0,c0,b0,Ropt,x00,T,...
    r,tBd,I,V,W,N0,opts)
%function takes in necessary inputs for OP1_costfun and optimization
%   options (opts) at a fixed temperature
%function outputs the solution of an optimal predator community (yOpt) and
%   its cost (J)
x0=x00*exp(-.69/(8.917e-5*T)); %compute metabolic constant at current temp
th0=N0(2:end); %true predator community
Acon=ones(1,length(th0)); %for Ax=b, use A "I sum up over all predators"
Bcon=Acon*th0; %for Ax=b, use b "the abundance in true community"
thLB=zeros(size(th0)); %lower bound is 'none of a predator'
thUB=ones(size(th0))*Bcon; %upper bound is 'total size of true community'
minProb=createOptimProblem('fmincon','objective', ...
    @(th) OP1_costFun(I,V,W,a0,phi,v0,h0,c0,b0,Ropt,x0,r,[N0(1);th],tBd),...
    'x0',th0,'lb',thLB, 'ub',thUB,'Aeq',Acon,'beq',Bcon,'options',opts);
%create an optimization problem with constraints lisited above
ms=MultiStart;
%create a multistart problem
[yOpt,J]=run(ms,minProb,100); %solve the optimization 100 times
yOpt=round(100*(yOpt./sum(yOpt))); %report 'percent of optimal predator community'

%option: plot the solutions you get
% [A, H, R,X]=OP1_structure(I,W);
% [a,ah,ba,x]=OP1_parameterization(a0,phi,v0,h0,b0,Ropt,x0,A,H,R,V,X);
% [t,ys]=ode45(@OP1_DE,tBd(1):1:tBd(2),[N0(1); yOpt],[],a,ah,ba,r,x);
% figure, hold on,
% plot(t,ys(:,1),'r--'), xlabel('time'), ylabel('aphid pop'),
% figure, hold on,
% for k=1:9
%    subplot(3,3,k), hold on,
%    plot(t,ys(:,k+1),'k-'), xlabel('time'), ylabel('pred pop'),
%    title(spec{k+1}),
% end
end

function [J,t,ys]=OP1_costFun(I,V,W,a0,phi,v0,h0,c0,b0,Ropt,x0,r,N0,tBd)
%function takes input from "OP1_initializingFun" (IM,TV,W,r), parameter
%   values (a0,phi,v0,h0,c0,b0,Ropt,x0), initial abundance (N0 - search over
%   this), and the  times to consider (tBd)
%function outputs the cost (J), and solutions (ys, at times t)
[A, H, R,C,X]=OP1_structure(I,W); %generate important matrices
[a,ah,c,ba,x]=OP1_parameterization(a0,phi,v0,h0,c0,b0,Ropt,x0,A,H,C,R,V,X); 
%parameterize model
ba=ba-diag(diag(ba)); %since we have c, remove interference penalty with self
[t,ys]=ode45(@OP1_DE,tBd(1):1:tBd(2),N0',[],a,ah,c,ba,r,x); %solve the ODE
J=sum(ys(:,1))/length(t); %compute average daily prey abundance
end

function [J,t,ys]=OP1_costTvar(I,V,W,a0,phi,v0,h0,c0,b0,...
    Ropt,x0,r,N0,tBd,Tf)
%function takes input from "OP1_initializingFun" (IM,TV,W,r), parameter
%   values (a0,phi,v0,h0,c0,b0,Ropt,x0), initial abundance (N0 - search over
%   this), and the  times to consider (tBd)
%function outputs the cost (J), and solutions (ys, at times t) USING TIME
%VARYIN TEMPERATURE
[A, H, R,C,X]=OP1_structure(I,W); %generate important matrices
[a,ah,c,ba,x]=OP1_parameterization(a0,phi,v0,h0,c0,b0,...
    Ropt,x0,A,H,C,R,V,X); %parameterize model
ba=ba-diag(diag(ba)); %since we have c, remove interference penalty with self
[t,ys]=ode45(@OP1_DETvar,tBd,N0',[],a,ah,c,ba,r,x,Tf);%solve the ODE
J=sum(ys(:,1))/length(t);%compute average daily prey abundance
end

%% FUNCTIONS THAT RUN THE MODEL
function dy=OP1_DE(~,y,a,ah,c,ba,r,x)
%function takes input current population (y) and parameters from
%   "OP1_parameterization"(a,ah,c,ba,r,x)
%function outputs step for ODE solver
F=ones(size(y))+(ah'+ba)*y+c'.*y; %compute functional response
effy=y./F; %scale predators by response
dy=y.*(r-a*effy-x'); %compute ODE term

% FINDING PARAMETER RANGES - ran this to see what c0, b0 are comparable
% Cs=c'.*y; Cs=Cs(2:end); %compute c0 contribution to predator FRs
% Bs=ba*y; Bs=Bs(2:end); %combute b0 contribution to predator FRs
% comps=[Cs,Bs]; %store in a matrix
% mean(comps), %report average over all predators
%Results: check contribution to FR close to equiv for baseline preds
%mean for c0=1, b0=10: 4.2772e+00   4.4233e+00
%for c0=.5 b0=5 : 2.1386e+00   2.2117e+00
%for c0=0.001; b0=.01 : 4.2772e-03   4.4233e-03
%range will be c0=[0,1] and equivalently b0=[0,10]
end
function dy=OP1_DETvar(t,y,a,ah,c,ba,r,x,Tf)
%function takes input current population (y) and parameters from
%   "OP1_parameterization"(a,ah,c,ba,r,x) AS WELL AS CURRENT TIME (t)
%function outputs step for ODE solver DEPENDENT ON VARYING TEMPERATURE
ind1=round(t*24*60/15)+1; %find the nearest 'day'
t1=(ind1-1)*15/24/60; %and the index before that
T=Tf(ind1)+(Tf(ind1+1)-Tf(ind1))*(t-t1)+273.25; 
%interpolate the current temperature
x=x*exp(-.69/(8.917e-5*T)); %compute metabolic death at this temperature
F=ones(size(y))+(ah'+ba)*y+c'.*y; %compute functional response
effy=y./F; %scale predators by response
dy=y.*(r-a*effy-x'); %compute ODE
end

function [a,ah,c,ba,x]=OP1_parameterization(a0,phi,v0,h0,c0,b0,Ropt,x0,...
    A,H,C,R,V,X)
%function takes in parameters (a0, phi, v0, h0, c0, b0, Ropt, x0) from the
%   ATN model along with matrices computed in "OP1_structure" (A,H,C,R,V,X)
%function outputs parameterized terms for attacks (aij), handling+attacks
%   (hij), competition (c), interference+attacks (ba), death (x)
Rmat=repmat(Ropt,length(Ropt),1); %stretch out Ropt values
Rm=R./Rmat; %compute Wj/Wi/Ropt
v=ones(size(V))-v0.*V; %compute overlap matrix
rat=(Rm.*exp(1-Rm)).^phi; %compute success curve in aij
a=a0.*A.*rat.*v; %compute aij
ah=h0.*a0.*H.*rat; %compute aij*hij
ba=b0*a; %compute b0*aij
c=c0*C; %compute cj
x=x0*X; x(1)=0; %compute death rate xj, but ignore for prey
end

function [A, H, R, C, X]=OP1_structure(I,W)
%function takes in interaction matrix (I) and masses (W)
%function outputs the base of ATN model (without parameters): for aij 
%   (A - outside parens, R - inside parens), hij (H), cj (C), xj (X)
W2=W.^(1/2); 
W4=W.^(1/4); %powers appearing in model
A=(W4'*W4).*I; %encounter term for aij. I makes this zero if no consumption occurs
H=repmat(W2',1,length(I)).*I; %for handling time
R=((W.^(-1))'*W).*I; %success term for aij
X=W.^(3/4); %metabolic death term
C=W.^(1/2); %competition term
end

%% BOOK KEEPING FUNCTIONS

function Ntest=OP1_genPred(pKeep,opt,pab)
%function takes in predator indices to use (pKeep), an optimal predator
%community (opt), and frequency of observed temperatures (pab)
%function outputs a group of predator communities to test, using different
%   weightings / selections of predators
[~,~,~,y0,~,~,~]=OP1_initializingFun(1,1,1,1,1); %find true abundances
BMA=sum(y0(2:end)); %total predator abundance

dt=BMA/sum(y0(pKeep)); %average (preserving true abundance) for preds in opt1
N0t=zeros(size(y0)); %initialize
N0t(1)=y0(1); %store prey abundance
N0t(pKeep)=dt*y0(pKeep)'; %store 'equal abundance' using preds from opt1

abun=opt*pab'/sum(opt*pab'); %scale opt1 by relative frequency of temps
N0w=[y0(1);abun.*sum(y0(2:end))]; %store temperature-scaled predator distribution

N0P=zeros(length(y0),9); %initialize empty vectors
for p=1:9 %for all predators
    N0P(1,p)=y0(1); %store the prey
    N0P(1+p,p)=BMA; %and put all available abundance into this predator
end
Ntest=[y0,N0t,N0w,N0P]; %this is the vector of initial abundances to test

end

function Div=OP1_diversityCounter(resies)
%function takes in predator community over a range of temperatures (resies)
%function outputs the average number of species (Div)
binres=resies~=0; %where is the result nonzero 
%(includes all temperature trials)
diversity=sum(binres); %sum up over all nonzero results
Div=sum(diversity)./length(diversity); 
%divide your sum over the number of trials
end

%% FUNCTIONS THAT HANDLE EXTERNAL DATA
function Tf=OP1_monthtemp(lib,day,tf)
%function takes in a library of temperature chunks (lib) of length (day)
%   and a final time to consider (tf)
%function outputs a sample stream of temperatures (Tf) which uses chunks
%   from lib to populate tf days of data
tall=length(lib(:,1)); %how much data I have
long=length(lib(1,:));
ext=tf/day; %how many chunks I need to fill
ind=randi([1,long],ext,1); %random indices to use
Tf=zeros(tall*ext,1); %initialize
count=1;
for j=1:(ext) %until I have all chunks
    Tf(count:count+tall-1)=lib(:,ind(j)); 
    %populate my temperature with the random chunks
    count=count+tall; %move forward
end
end

function [Chunky,Smooth]=OP1_tempgen(day)
%function takes as input the desired length of chunks (day)
%function outputs the temperature from complete fields cut into chunks of
%   that many days (Chunky) and also stored in a single vector (Smooth)
fields={'JC','JO','KC','KO','MC','MO','OC','OO','SC','SO'}; %names of fields
Chunky=[]; Smooth=[]; %initialize
for j=[1:8,10] %loop over fields
    %skip 9, because logger in SC goes missing partway through season
    [chunky,smooth]=OP1_choptemp(fields{j}, day); %cut up the field's temperatures
    Chunky=[Chunky,chunky]; %store with prior fields
    Smooth=[Smooth;smooth]; %store with prior fields
end
end

function [chunky,smooth]=OP1_choptemp(field, day)
%function takes in current field and number of days
%function outputs the temperatures (smooth) and recording times (d) in that
%   field and a matrix (chunk) with temperatures cut into chunks of length 'day'
A=load(['Templog numDate comma ',field,'.csv']); %load raw data
smooth=A(:,2); %temperatures
l=length(smooth); %length of full temps
dind=96; %how many indices in a day
ch=day*dind; %how long a 'chunk' of temperatures will be
chunky=[]; count=1; %initialize
while count+dind+ch-1 < l %while you can still take new chunks
    chunk=smooth(count:count+ch-1); %cut out a chunk of temperatures
    chunky=[chunky,chunk]; %store the temperature chunks
    count=count+dind; %move forward
end
end

function [IM,TV,W,y0,Ropt,rs,spec]=OP1_initializingFun(RP,RB,PB,BB,PP)
%function takes in overlap parameters for prey-spiders (RP), prey-beetles
%   (RB), beetles-spiders (PB), beetles-beetles (BB), and spiders-spiders (PP)
%function outputs quantities needed for ATN model solutionos: interaction
%   matrix (IM), overlap matrix of nu values (TV), vector of masses (W),
%   initial abundancecs (y0), optimal predator-prey body mass ratios (Ropt),
%   prey growth rates (rs), and names for all species (spec).

fields={'JC','JO','KC','KO','MC','MO','OC','OO','SC','SO'}; %field refs
specs={'Aphid','Bembidion','Harpalus','Poecilus','Pterostichus','Other_Carabid',...
    'Linyphiidae','Lycosidae','Tetragnathidae','Other_Spider'}; %pred refs
rr=.3402; %prelim estimate - aphid growth rate from Wooton 2020 experiments

Wtot=zeros(1,10); ytot=zeros(10,1); %initialize vectors

for k=1:length(fields) %loop over fields
    Wfull=load(['BM FoodWeb Fulltime ',fields{k},'.csv']); %load the masses
    W=[Wfull(1),Wfull(6:14)']; %save masses for prey, beetles, and spiders
    y0=zeros(length(W),1); %initialize vector
    for j=2:length(specs) %loop over predators
        dats=load(['Density numDate ',specs{j},' Fulltime ',fields{k},'.csv']); %load abundances
        pop=dats(:,2); %discard times
        averagePop=mean(pop); %compute average abundance
        y0(j)=averagePop; %save this predator's average abundance
    end
    aphs=load(['Density numDate Aphid Shifted Truncated ',fields{k},'.csv']); %load abundances
    y0(1)=aphs(1,2); %save initial abundance
    Wtot=Wtot+W; %add up the masses we've found
    ytot=ytot+y0; %add up the initial abundances we've found
end
Wtot=Wtot./length(fields); %compute average mass over all fields

[WB,Bind]=sort(-1*Wtot(2:6)); %arrange beetles by size (easier to view)
[WS,Sind]=sort(-1*Wtot(7:end)); %arrange spiders by size (easier to view)
W=[Wtot(1),-1*WB,-1*WS]; %re-pack all species

ytot=ytot./length(fields); %compute average init. abund. over all fields
yB=ytot(2:6); yS=ytot(7:end); %separate beetles/spiders
y0=[ytot(1);yB(Bind);yS(Sind)]; %pack using ordering from above (by mass)

specs={'Aphid','Bembidion','Harpalus','Poecilus','Pterostichus','Other Carabid',...
    'Linyphiidae','Lycosidae','Tetragnathidae','Other Spider'}; %readable species names

spec={specs{Bind(1)+1},specs{Bind(2)+1},specs{Bind(3)+1},specs{Bind(4)+1},specs{Bind(5)+1},...
    specs{Sind(1)+6},specs{Sind(2)+6},specs{Sind(3)+6},specs{Sind(4)+6}}; 
%pack using ordering from above (by mass)

TVB=BB*ones(5,5)-BB*eye(5,5); %overlap matrix for beetles
TVP=PP*ones(4,4)-PP*eye(4,4); %overlap matrix for spiders
TVBP=PB*ones(5,4); %overlap matrix for beetles-spiders
TVr=[RB,RB,RB,RB,RB,RP,RP,RP,RP]; %overlap matrix for prey
TV=[0,TVr;TVr(1:5)',TVB,TVBP;TVr(6:end)',TVBP',TVP]; %pack all overlap

IMfull=load('FoodWeb.csv'); %load interaction matrix
IMhalf=[IMfull(1,:);IMfull(6:14,:)]; %save the part for prey, beetles, spiders
IM=[IMhalf(:,1),IMhalf(:,6:14)];
IM(9,2)=1; %real foodweb introduced inconsistencies
%the only non-consumption is Bembidion won't eat tetragnathidae, 
%allow this link to make comparison possible between groups
Ropt=W./W(1); %set optimal predator-prey body mass ratio to "common prey is optimal"
rs=zeros(length(W),1); rs(1)=rr; %pack only prey growth rate
end

