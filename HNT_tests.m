%% HNT tests

function [v_max,khat]=VbarStat(yt)
% H&H 2012 CUSUM statistics to detect change point in panel data.
% yt - data of interest
% v_max - test statistic
% khat - detected change point
[T,N]=size(yt);
for i=1:N
    dmy(:,i)=yt(:,i)-mean(yt(:,i));
end
StT=sum(dmy,1);
vbar=nan(T,1);

% renormalization
[zt,~,~]=CP_recenter(yt);
zt=zt-repmat(mean(zt,1),T,1);
z_mat=(zt.^2)-repmat(mean(zt.^2,1),T,1);
bw = floor(log(T));
V=mycovnw(z_mat,bw,1,'hacc_f');
w_bar=sum(sum(V));


for i=1:T
    Sti=sum(dmy(1:i,:),1);
    Zti=(Sti-(i/T)*StT);
    vbar(i,1)=sum((Zti.^2)-(i*(T-i)/(T^2)))/sqrt(w_bar);
end

	[v_max,khat]=max(abs(vbar));

end

%% critical values derived using the bootstrap method (1)
function [cv_bboot]=gfmbootstrap_cv(data,rep)
% Bootstrap the empirical critical values from the data
% data - data of interest
% rep - number of replications
[T,N]=size(data);% T>=500, e.g. block_len=50

[data,~,~]=CP_recenter(data);
%% block bootstrap (moving block-non-overlapping) critical values
Vnbb_con=zeros(rep,1);

parfor i=1:rep
    Bsample=gfm_boots(data);
    [vnt_max,~]=VbarStat(Bsample);
    Vnbb_con(i)=vnt_max;
end

cv_bboot=zeros(3,1);
cv_bboot(1,1)=quantile(Vnbb_con,0.90);
cv_bboot(2,1)=quantile(Vnbb_con,0.95);
cv_bboot(3,1)=quantile(Vnbb_con,0.99);

end




function [boot_data]=gfm_boots(data)
% this function boostrap a panel dataset by using the generalized dynamic factor model-parametric boostrap method (Cho, 2016)
% the number of common factors are estimated by Bai & Ng (2002)
[T,N]=size(data);
%sd=std(data);% SD of each column
%ndata=data./repmat(sd,100,1);

ndata=nan(T,N);
bw = floor(log(T));
V=mycovnw(data,bw,1,'hacc_b');
lr=diag(V);
for i=1:N
ndata(:,i)=data(:,i)/sqrt(lr(i));
end

DEMEAN=1;%=2 standardized data
rmax=10;
BN_method=1; % range=1:7

%% step 1-3, get common factors
 
[rhat, ~, ~]=nbplog(ndata,rmax,BN_method,DEMEAN);%nbpiid(x,rmax,BN_method,DEMEAN);

if rhat>0
[eigv,~] = eig((ndata'*ndata)/T);
lam=eigv(:,1:rhat)*sqrt(N);
[coeff,~]=pca(data);
cf=data*coeff(:,1:rhat)/sqrt(N);
%reminder=ndata-(lam*cf')';

%% step 4.1, draw independent common factors
ind=[];
m=floor(log(T))-1;
while length(ind)<T
    I = randsample(T-m,1);%-m
    I = I:I+m-1;
    ind=[ind, I];
end
ind=ind(1:T);
bf=cf(ind,:);
bc=(lam*bf')';

%% step 4.2, boostrap reminder sequence using Time-frequency Toggle-boostrap method
bi=mvnrnd(zeros(1,N),eye(N),100+T);
bi=bi(101:end,:);
%% combine components and get boostrapped data

boot_data=bc+bi;

for i=1:N
boot_data(:,i)=boot_data(:,i)*sqrt(lr(i));
end

else

bi=mvnrnd(zeros(1,N),eye(N),100+T);
bi=bi(101:end,:);

boot_data=bi;
for i=1:N
boot_data(:,i)=boot_data(:,i)*sqrt(lr(i));
end

end

end

%% critical values derived using the bootstrap method (2)

function [cv_bboot]=newbootstrap_cv(data,rep)
% Bootstrap the empirical critical values from the data
% data - data of interest
% rep - number of replications
[T,N]=size(data);% T>=500, e.g. block_len=50

[data,~,~]=CP_recenter(data);
%% block bootstrap (moving block-non-overlapping) critical values
Vnbb_con=zeros(rep,1);

parfor i=1:rep
    Bsample=new_boots(data);
    [vnt_max,~]=vnt_stat(Bsample);
    Vnbb_con(i)=vnt_max;
end

cv_bboot=zeros(3,1);
cv_bboot(1,1)=quantile(Vnbb_con,0.90);
cv_bboot(2,1)=quantile(Vnbb_con,0.95);
cv_bboot(3,1)=quantile(Vnbb_con,0.99);

end





function [boot_data]=new_boots(data)
% this function boostrap a panel dataset with estimated sigma_bar and c_star (our paper)
% the number of common factors are estimated by Bai & Ng (2002)
[T,N]=size(data);
data=demean(data);

%de-volatize
bw = floor(log(T));
%estimate the number of common factors

DEMEAN=1;%=2 standardized data
rmax=10;
BN_method=1; % range=1:7
[rhat, ~, ~]=nbplog(data,rmax,BN_method,DEMEAN);%nbpiid(x,rmax,BN_method,DEMEAN);


%%find sigma_bar and c_star
% if there is a common factor(s)
if rhat>0

[eigv,~] = eig((data'*data)/T);
lam=eigv(:,1:rhat)*sqrt(N);

[coeff,score]=pca(data);
cf_hat=data*coeff(:,1:rhat)/sqrt(N); 
reminder=data-(lam*cf_hat')';


sig4=0;
for i=1:N
    %sig4=sig4+mean(reminder(:,i).^2)^2;
    sig4=sig4+(mycovnw(reminder(:,i),bw,1,'hacc_b')^2);
end
sigma_bar=(sig4/N)^(1/2);


% construct bootstrap data
Vcf_hat=mycovnw(cf_hat,bw,1,'hacc_b');
cf_sim=mvnrnd(zeros(1,rhat),(Vcf_hat)^(1/2),100+T);
cf_sim=cf_sim(101:end,:);

edata=normrnd(0,sqrt(sigma_bar),[(T+100),N]);
edata=edata(101:end,:);

lam_sim=nan(1,rhat);
for i=1:rhat
c_star=sum(abs(lam(:,i)).^2)/sqrt(N);
lam_sim(i)=sqrt(c_star*sqrt(N)/N);
end

cfdata=0;
for i=1:rhat
    cfdata=cfdata+repmat(cf_sim(:,i),1,N)*lam_sim(i);
end
boot_data=cfdata+edata;

% if there is no common factor, i.e. c_star=0
else

sig4=0;
for i=1:N
    %sig4=sig4+mean(data(:,i)^2)^2;
    sig4=sig4+(mycovnw(data(:,i),bw,1,'hacc_b')^2);
end
sigma_bar=(sig4/N)^(1/2);

boot_data=normrnd(0,sqrt(sigma_bar),[(T+100),N]);
boot_data=boot_data(101:end,:);

end

end



%% utility

function [r_yt,dmy,kv]=CP_recenter(yt)
% recentered r_yt
[T,N]=size(yt);

%% find change point khat
%dmy=yt-repmat(mean(yt,1),T,1);
dmy=1;

r_yt=nan(T,N);
kv=nan(N,1);
for i=1:N
    [~, ~, khat] = CUSUMtest(yt(:,i));
    khat=round(khat);
    r_yt(1:khat,i)=yt(1:khat,i)-repmat(mean(yt(1:khat,i),1),khat,1);
    r_yt((khat+1):end,i)=yt((khat+1):end,i)-repmat(mean(yt((khat+1):end,i),1),T-khat,1);
    kv(i)=khat;
end

end




function [Stat, stat, ind] = CUSUMtest(series)
% one function performs non-parametric CUSUM test to detect change point.
n=length(series);
Tn=zeros(n,1);
for i=1:n
    tn=(sum(series(1:i,1))-(i/n)*sum(series))/sqrt(n);
    Tn(i,1)=tn;
end

stat=zeros(n,1);
t=1:n;
t=t./n;
w=sqrt(t.*(1-t));
w=w';


Var=mycovnw(series,floor(log(n)),1,'hacc_b');
for i=1:n
    stat(i,1)=1/sqrt(Var)*abs(Tn(i,1))./w(i,1);
end

[Stat,ind]=max(stat);

if ind<=0.15*n
   ind=0.15*n;
end

if ind>=0.85*n
   ind=0.85*n;
end    

end




function x=demean(X);

[T,N]=size(X);
m=mean(X);
x=X-repmat(m,T,1);

end


function pos=minindc(x);
ncols=size(x,2);
nrows=size(x,1);
pos=zeros(ncols,1);
seq=seqa(1,1,nrows);
for i=1:ncols;
dum=min(x(:,i));
dum1= seq .* ( (x(:,i)-dum) ==0);
pos(i)=sum(dum1);
end;





function V=mycovnw(data,nlag,demean,type)
T=size(data,1);
if nargin==1
    nlag=min(floor(1.2*T^(1/3)),T);
    demean=true;
elseif nargin==2
    demean=true;    
end    
if isempty(nlag)
    nlag=min(floor(1.2*T^(1/3)),T);
end
if isempty(demean)
    demean=true;
end
if ~ismember(demean,[0 1]) 
    error('DEMEAN must be either logical true or false.')
end
%if floor(nlag)~=nlag || nlag<0 
%    error('NLAG must be a non-negative integer.')
%end
nlag=floor(nlag);
if ndims(data)>2
    error('DATA must be a T by K matrix of data.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if demean
    data=data-repmat(mean(data),T,1);
end
 
% NW weights
switch lower(type) 
    case{'hacc_b'}
   a=[1/(nlag):1/(nlag):1];
   a=[0,a];
   w=1-a;
% Truncated weights
    case{'hacc_t'} 
w=ones(1,nlag+1);
% Parzen weights
    case{'hacc_p'} 
   a=[1/(nlag):1/(nlag):1];
   a=[0,a];
   w=ones(1,nlag+1);
   for i=1:nlag+1
       if a(1,i)<=0.5
          w(1,i)=1-6*a(1,i)^2+6*abs(a(1,i))^3;
       else
          w(1,i+1)=2*(1-abs(a(1,i)))^3;
       end
   end
    case{'hacc_th'} 
   a=[1/(nlag):1/(nlag):1];
   a=[0,a];
   w=ones(1,nlag+1);
   for i=1:nlag+1
       w(1,i)=(1+cos(pi*a(1,i)))/2;
   end  
% Quadratic Spectral weights
    case{'hacc_qs'}
   a=[1/(nlag):1/(nlag):1];
   a=[0,a];
   w=ones(1,nlag+1);
   for i=1:nlag
       w(1,i+1)=(25/(12*pi^2*a(1,i+1)^2))*(sin(6*pi*a(1,i+1)/5)/(6*pi*a(1,i+1)/5)-cos(6*pi*a(:,i+1)/5));
   end
% Flat Top weights
    case{'hacc_f'}
    a=[1/(nlag):1/(nlag):1];
    a=[0,a];
    w=ones(1,nlag);
    for i=1:nlag
        if a(1,i)<=0.05
        w(1,i)=1;
        else
        w(1,i)=exp(-1*exp(-1/(a(1,i)-0.05)^2)/(a(1,i)-1)^2);
        end
    end
   w=[w,0];
end

% Start the covariance
V=data'*data/T;
for i=1:nlag
    Gammai=(data((i+1):T,:)'*data(1:T-i,:))/(T-i);
    GplusGprime=Gammai+Gammai';
    V=V+w(i+1)*GplusGprime;
end

end





function [ic1, chat,Fhat]=NBPLOG(x,kmax,jj,DEMEAN);
T=size(x,1);
N=size(x,2);
NT=N*T;
NT1=N+T;
CT=zeros(1,kmax);
ii=1:1:kmax;
if jj ==1;CT(1,:)=log(NT/NT1)*ii*NT1/NT;end;
if jj==2; CT(1,:)=(NT1/NT)*log(min([N;T]))*ii;end;
GCT=min([N;T]);
if jj==3; CT(1,:)=ii*log(GCT)/GCT; end;
if jj==4; CT(1,:)=2*ii/T; end;
if jj==5; CT(1,:)=log(T)*ii/T;end;
if jj==6; CT(1,:)=2*ii*NT1/NT; end;
if jj==7; CT(1,:)=log(NT)*ii*NT1/NT;end;

 if DEMEAN ==2;
 X=standard(x);
 end;

if DEMEAN ==1;
 X=demean(x);
 end;
if DEMEAN==0;
  X=x;;
  end;

IC1=zeros(size(CT,1),kmax+1);
Sigma=zeros(1,kmax+1);
XX=X*X';
[Fhat0,eigval,Fhat1]=svd(XX');
for i=kmax:-1:1;
Fhat=Fhat0(:,1:i);
lambda=Fhat'*X;
chat=Fhat*lambda;
ehat=X-chat;
Sigma(i)=mean(sum(ehat.*ehat/T));
IC1(:,i)=log(Sigma(i))+CT(:,i);
end;
Sigma(kmax+1)=mean(sum(X.*X/T));
IC1(:,kmax+1)=log(Sigma(kmax+1));
ic1=minindc(IC1')';
ic1=ic1 .*(ic1 <= kmax);
Fhat=[];
Fhat=Fhat0(:,1:ic1);
lambda=Fhat'*X;
chat=Fhat*lambda;
end




function x=standard(y);
T=size(y,1);
N=size(y,2);
my=repmat(mean(y),T,1);
sy=repmat(std(y),T,1);
x=(y-my)./sy;

end