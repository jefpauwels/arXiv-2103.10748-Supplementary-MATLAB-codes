clear
clear cvx
clc 

%% Basic parameters

prec = 10^(-7); % precision stop condition
format long

d=4; %dimension Bob's system
N=5; %number of states
m=4; %number of measurements

MaxIter = 200; % maximal number of seesaw iterations
Nloops = 10; % number of times to execute seesaw

Ns = 4; %alphabet size

%% Witness

% Gallego

for r = 1:N
	for o= 1:N-1
		if r <= N - o 
			c(r,o) =1;
		elseif o == N -r +1
		c(r,o) = -1;
		end
	end
end

% RAC

% c = [ 1 1 ; 1 -1; -1 1; -1 -1];

%% Seesaw 

optvals = zeros(Nloops,1);
for loop = 1:Nloops

% tech parameters
dat= [ 0 ];
k=0; % initialise true

i=1; %counter for data list



for y = 1:m
    for s = 1:Ns
        U = RandomUnitary(d);
        E(:,:,y,s) = 2*U*diag([1 1 0 0])*U'- eye(d); %adapt to alphabet size
    end
end


while k == 0
    i = i +1;

cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best
variable rho(d,d,N,Ns) hermitian semidefinite

obj = 0;
for x = 1:N 
    for y = 1:m
        for s = 1:Ns
            obj=obj+c(x,y)*(trace(rho(:,:,x,s)*E(:,:,y,s)));
        end
    end
end


maximise(obj)

subject to


for x= 1:N
    for xx = 1:N
        if x>xx
    if Ns==2
        trace(rho(:,:,x,1)) + trace(rho(:,:,x,2)) == 1
        rho(:,:,x,1) + rho(:,:,x,2)  == rho(:,:,xx,1) + rho(:,:,xx,2) 
    elseif Ns==3
        trace(rho(:,:,x,1)) + trace(rho(:,:,x,2)) +  trace(rho(:,:,x,3)) ==1
        rho(:,:,x,1) + rho(:,:,x,2) +rho(:,:,x,3)  == rho(:,:,xx,1) + rho(:,:,xx,2) +rho(:,:,xx,3)
    elseif Ns==4
    trace(rho(:,:,x,1)) + trace(rho(:,:,x,2)) +  trace(rho(:,:,x,3))+ trace(rho(:,:,x,4))==1
    rho(:,:,x,1) + rho(:,:,x,2) +rho(:,:,x,3) +rho(:,:,x,4)  == rho(:,:,xx,1) + rho(:,:,xx,2) +rho(:,:,xx,3) +rho(:,:,xx,4)
    end
        end
    end
end


cvx_end

    
cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best
variable E(d,d,m,Ns) hermitian


obj = 0;
for x = 1:N 
    for y = 1:m
        for s = 1:Ns
            obj=obj+c(x,y)*(trace(rho(:,:,x,s)*E(:,:,y,s)));
        end
    end
end

maximise(obj)
subject to

for y = 1:m
    for s = 1:Ns
        E(:,:,y,s) >= -eye(d)
       E(:,:,y,s) <= eye(d)
    end
end


cvx_end


dat=[dat; cvx_optval]


if dat(i) - dat(i-1) < prec
    k =1; optvals(loop)= dat(i);
elseif i > MaxIter
    k =1;
elseif  isnan(dat(i)) 
    k=1; 
else
  k= 0;
end

end

end


