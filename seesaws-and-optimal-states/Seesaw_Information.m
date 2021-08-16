clear all
clc


%% Setup problem

nX=5; % number of preparation inputs
nY=4; % number of measurement inputs
nO=2; % number of measurement outcomes

q=ones(1,nX)/nX; % prior distribution


coef=[1 1 1 1; 1 1 1 -1; 1 1 -1 0; 1 -1 0 0;-1 0 0 0];


D=8; % dimension of the message 

%% See-saw
    
for y=1:nY
        B(:,:,y,1)=RandomDensityMatrix(D);
        B(:,:,y,2)=RandomDensityMatrix(D);
end


for K=1:20

cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best
variable rho(D,D,nX) hermitian semidefinite
variable sig(D,D) hermitian semidefinite


obj=0;
for x=1:nX
    trace(rho(:,:,x))==1;
    q(x)*rho(:,:,x) <= sig;
    for y=1:nY
        for b=1:nO
        obj=obj+coef(x,y)*(-1)^(b-1)*trace(rho(:,:,x)*B(:,:,y,b));
        end
    end
end

subject to 
trace(sig)<=2/5;

maximise(real(obj))

cvx_end
cvx_optval




%%%%%%%%%%%%%%%%%%%%%%%%%



cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best
variable B(D,D,nY,nO) hermitian semidefinite

obj=0;
for x=1:nX
    for y=1:nY
        for b=1:nO
        obj=obj+coef(x,y)*(-1)^(b-1)*trace(rho(:,:,x)*B(:,:,y,b));
        end
    end
end


for y=1:nY
    s=0;
    for b=1:nO
        s=s+B(:,:,y,b);
    end
    s==eye(D);
end


maximise(real(obj))

cvx_end
cvx_optval
end

