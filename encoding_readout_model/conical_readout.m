function [choice,m1,m2,region] = conical_readout(r1,r2,alpha,p_conf)

if nargin<4
    p_conf = 'ideal';
end

% alpha is expressed in radians and it is the angle indicating the aperture
% of the cone

% beta makes sense between -pi/4 and pi/4;
beta=.5*(alpha-pi/2);

xp = r1-.5;
yp = r2-.5;

m1 = -1*(cos(beta)/sin(beta));
m2 = -1*(sin(beta)/cos(beta));

choice=nan*ones(size(r1));

if m1<0
    above_sb1 = yp-m1*xp>0;
    above_sb2 = yp-m2*xp>0;
elseif m1>0
    above_sb1 = yp-m1*xp<0;
    above_sb2 = yp-m2*xp>0;
end


% choice(above_sb1 & above_sb2) = 1;
% choice(~above_sb1 & ~above_sb2) = 0;

% define region which trials belong to
% 0: outside the cone
% 1: region where choice=0
% 2: region where choice=1
region = zeros(size(r1));
region(~above_sb1 & ~above_sb2)=1;
region(above_sb1 & above_sb2)=2;

% behavior on the boundaries
region(abs(yp-m1*xp)<10e-15)=0;
region(abs(yp-m2*xp)<10e-15)=0;

% decoded stimulus
sp = zeros(size(r1));
sp(xp+yp>0)=1;

% for the nice figure: .6 .4 - .9 .1

switch p_conf
    case 'ideal'
        choice(region==0 & sp==0)=randsample([0,1],nnz(region==0 & sp==0),true,[.5 .5]);
        choice(region==0 & sp==1)=randsample([0,1],nnz(region==0 & sp==1),true,[.5 .5]);
        choice(region==1)=randsample([0,1],nnz(region==1),true,[1 0]);
        choice(region==2)=randsample([0,1],nnz(region==2),true,[0 1]);
    case 'strong'
        choice(region==0 & sp==0)=randsample([0,1],nnz(region==0 & sp==0),true,[.6 .4]);
        choice(region==0 & sp==1)=randsample([0,1],nnz(region==0 & sp==1),true,[.4 .6]);
        choice(region==1)=randsample([0,1],nnz(region==1),true,[.9 .1]);
        choice(region==2)=randsample([0,1],nnz(region==2),true,[.1 .9]);
    case 'nonsense'
        choice(region==0 & sp==0)=randsample([0,1],nnz(region==0 & sp==0),true,[.4 .6]);
        choice(region==0 & sp==1)=randsample([0,1],nnz(region==0 & sp==1),true,[.6 .4]);
        choice(region==1)=randsample([0,1],nnz(region==1),true,[.9 .1]);
        choice(region==2)=randsample([0,1],nnz(region==2),true,[.1 .9]);
end








