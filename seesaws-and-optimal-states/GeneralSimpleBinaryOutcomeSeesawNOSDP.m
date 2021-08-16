%% Simple seesaw for scenarios with binary outcome measurements (no SDPs)
% Inputs:
% - c = correlation witness
% - d = dimension of the quantum message
% - r = a flag, 0 corresponding to complex strategies, 1 to real strategies
% - NoLoops = number of times to run the seesaw algorithm
% - prec = stopping precision of the algorithm
function [witnessval rho M] = GeneralSimpleBinaryOutcomeSeesawNOSDP(c,d,r,NoLoops,prec)

N = size(c,1);
m= size(c,2);

% measurements ranks

for y = 1:m
    rankpos(y) = 0;
    for x = 1:N
        if c(x,y) > 0
            rankpos(y)= rankpos(y) +1;
        end
    end
    dg{y} = zeros(1,d);
        for l = 1:min(rankpos(y),d)
            dg{y}(l) = 1;
    end
end
    
% seesaw loop 
loopobj = zeros(NoLoops,1);

for loop = 1: NoLoops
k=0;
dat = [ 0 ];
i = 1;

%initialize
 for y = 1:m   
    U = RandomUnitary(d,r);  %make real if needed (add flag 1)
    M{y,1} = U*diag(dg{y})*U';
    M{y,2} = eye(d) - M{y,1};
    E{y} = 2*M{y,1} - eye(d);
    E{y} = (E{y}+E{y}')/2;
 end

 for x = 1:N
    W{x} = zeros(d,d);
    for y = 1:m
        W{x} = W{x} + c(x,y)*E{y};
    end
end



while k == 0
    i = i +1;
    

for x = 1:N
    [V{x} L{x}] = eig(W{x});
    psi{x} = V{x}(:,d);
    rho{x} = psi{x}*psi{x}';
    rho{x} = (rho{x}+rho{x}')/2;
end

for y = 1:m
    T{y} = zeros(d,d);
    for x = 1:N
        T{y} = T{y} + c(x,y)*rho{x};
    end
    [VM{y} LM{y}] = eig(T{y});
end


for y=1:m
    M{y,1} = zeros(d,d);
    for l = 1:d
        if LM{y}(l,l) > 0
          M{y,1} = M{y,1} + VM{y}(:,l)*VM{y}(:,l)';
        end
    end
    E{y} = 2*M{y,1} - eye(d);
    E{y} = (E{y}+E{y}')/2;
end

for x = 1:N
    W{x} = zeros(d,d);
    for y = 1:m
        W{x} = W{x} + c(x,y)*E{y};
    end
end

S = 0;
for x = 1:N
    S = S + real(trace(W{x}*rho{x}));
end



dat=[dat; S];

if dat(i) - dat(i-1) < prec
    k =1; loopobj(loop) = S; 
else
  k= 0;

end


end


end
witnessval = max(loopobj);
end