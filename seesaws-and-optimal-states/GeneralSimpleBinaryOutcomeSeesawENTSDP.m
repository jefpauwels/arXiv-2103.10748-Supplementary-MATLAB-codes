%% Simple seesaw for scenario with binary outcome measurements and shared entanglement 
% Inputs:
% - c = correlation witness
% - dA = dimension of the quantum message
% - dB = dimension of the shared entanglement
% - r = a flag, 0 corresponding to complex strategies, 1 to real strategies
% - NoLoops = number of times to run the seesaw algorithm
% - prec = stopping precision of the algorithm
% - ranks should be the best guess for the ranks of the measurement
% - varargin is an optional parameter for qubit-qubit entanglement of a
% particular entangled state (theta). If you don't want to use 

function [witnessval rho M] = GeneralSimpleBinaryOutcomeSeesawENTSDP(c,dA,dB,r,NoLoops,prec,ranks,varargin)  

if size(varargin) > 0
    theta = varargin{1};
    chi = [cos(theta)^2 0 0 cos(theta)*sin(theta);0 0 0 0; 0 0 0 0; cos(theta)*sin(theta) 0 0 sin(theta)^2];
else
    fm = 0;
end


N = size(c,1);
m= size(c,2);
d=dA*dB;


% for y = 1:m
%     rankpos(y) = 0;
%     for x = 1:N
%         if c(x,y) > 0
%             rankpos(y)= rankpos(y) +1;
%         end
%     end
%     dg{y} = zeros(1,d);
%         for l = 1:min(rankpos(y),d)
%             dg{y}(l) = 1;
%     end
% end
%     

% for y = 1:m
%     NoPos(y) = 0;
%         for x = 1:N
%             if c(x,y) > 0
%                 NoPos(y) = NoPos(y) + 1;
%             end
%         end
%     dg{y} = zeros(1,d);
%         for l = 1:max(1,round(NoPos(y)/N)*d)
%             dg{y}(l) = 1;
%         end
% end

%custom rank

    
for y = 1:m
    dg{y} = zeros(1,d);
        for l = 1:ranks(y)
            dg{y}(l) = 1;
        end
end


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
    

cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best
variable rho(d,d,N) hermitian semidefinite

S = 0;
for x = 1:N
    S = S + real(trace(W{x}*rho(:,:,x)));
end

maximise(S)

for x = 1:N
    trace(rho(:,:,x)) == 1
    if fm == 0
        for xx = 1:N
            if xx > x
                 PartialTrace(rho(:,:,x),1,[dA,dB]) ==  PartialTrace(rho(:,:,xx),1,[dA,dB])
            end
        end
    else
        PartialTrace(rho(:,:,x),1,[dA,dB]) == PartialTrace(chi,1,[dA,dB]);
    end
end


cvx_end

for y = 1:m
    T{y} = zeros(d,d);
    for x = 1:N
        T{y} = T{y} + c(x,y)*rho(:,:,x);
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

S = 0;
for x = 1:N
    S = S + real(trace(W{x}*rho(:,:,x)));
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