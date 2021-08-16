clearvars
clear cvx
clc 

%% Basic parameters 

prec = 10^(-7); %precision stop condition
format long

d=9; %dimension Bob's system
N=9;
m=2;
MaxIter = 120;
Nloops = 10;

Ns = 3;

%% trit RAC witness

t = 0;
for x1=1:3
    for x2=1:3
        t = t + 1;
        counter(x1,x2) = t;
    end
end

for x1 = 1:3
    for x2 = 1:3
        for y = 1:2
            for b = 1:3
                if y==1
                    if x1 == b
                        c(counter(x1,x2),y,b) = 1;
                    end
                elseif y == 2
                    if x2 == b
                        c(counter(x1,x2),y,b) = 1;
                    end
                end
            end
        end
    end
end
                  
c=c/18;

%% seesaw loop

optvals = zeros(Nloops,1);
for loop = 1:Nloops

% tech parameters
dat= [ 0 ];
k=0; % initialise true

i=1; %counter for data list



for y = 1:m
    for s = 1:Ns
        U = RandomUnitary(d);
        M(:,:,y,s,1) = U*diag([1 0 0 0 0 0 0 0 0])*U'; %adapt to alphabet size
        U = RandomUnitary(d);
        M(:,:,y,s,2) = U*diag([1 0 0 0 0 0 0 0 0])*U';
        M(:,:,y,s,3) = eye(d) - M(:,:,y,s,2) - M(:,:,y,s,1);
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
        for b = 1:3
            for s = 1:Ns
                obj=obj+c(x,y,b)*(trace(rho(:,:,x,s)*M(:,:,y,s,b)));
            end
        end
    end
end


maximise(real(obj))

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
variable M(d,d,m,Ns,3) hermitian semidefinite



obj = 0;
for x = 1:N 
    for y = 1:m
        for b = 1:3
            for s = 1:Ns
                obj=obj+c(x,y,b)*(trace(rho(:,:,x,s)*M(:,:,y,s,b)));
            end
        end
    end
end

maximise(real(obj))
subject to

for y = 1:m
    for s = 1:Ns
      M(:,:,y,s,1) + M(:,:,y,s,2) + M(:,:,y,s,3) == eye(d)
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


