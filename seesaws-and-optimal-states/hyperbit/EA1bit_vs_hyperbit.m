%{
Proof that EA 1-bit models outperform hyperbit models
 Supplement to Correlations in entanglement-assisted prepare-and-measure scenarios
 by A. Tavakoli. J. Pauwels, E. Woodhead, S. Pironio
 2020
%}


%% Optimal QM strategy for the witness I5 defined in eq.(8) of Phys.Rev.Lett.105:230501,2010
%{
We load the optimal EA 1-bit strategy that we have found, defined by
    psi = shared entangled state in R^4 \otimes R^4
    A{x} = Alice's binary observable for input x (x=1,..,5)
    B{y,a} = Bob's binary observable for input y (y=1,..,4) and communicated bit a (a=1,..2) from Alice
%}
load qmstrategy

% We define the coefficients of I5
I5 = [1 1 1 1; 1 1 1 -1; 1 1 -1 0; 1 -1 0 0; -1 0 0 0]; 

% We compute the correlations E{x,y} according to the QM strategy
Id = eye(4); % the idendity matrix
E=zeros(5,4);
for x=1:5
    for y=1:4
        E(x,y) = psi'*kron((Id+A{x})/2,B{y,1})*psi + psi'*kron((Id-A{x})/2,B{y,2})*psi;
    end
end

% We compute the value of I_5
I5qm = sum(sum(I5.*E)),
% This yields I5qm = 9.0343

%% Maximal value of I5 for hyperbit models
% requires CVX

I5hyperbit = 0;
% ij = 0 --> Bob discards the hyperbit and returns b=1 for input y=j
% ij = 1 --> Bob measures the hyperbit for input y=j
for i1=0:1
    for i2=0:1
        for i3=0:1
            for i4=1:1
                I5h=[i1 i2 i3 i4; i1 i2 i3 -i4; i1 i2 -i3 0; i1 -i2 0 0; -i1 0 0 0];
                Obj=[zeros(5,5) I5h; I5h' zeros(4,4)];
                cvx_solver mosek
                cvx_begin quiet sdp 
                    variable G(9,9) symmetric semidefinite
                    maximize 3*(1-i1)+2*(1-i2)+1*(1-i3)+0.5*trace(Obj*G)
                    diag(G)<=[1 1 1 1 1 1 1 1 1]';
                cvx_end
                if cvx_optval>I5hyperbit I5hyperbit = cvx_optval; end
            end
        end
    end
end
% we find
I5hyperbit,
% I5hyperbit = 9<9.0343


%% Attempt to reconstruct a hyperbit model from the QM strategy following Phys. Rev. A 85, 022331 (2012)

% We build the vectors X and Y as in eq. (A8) and randomize Alice's message
% as in the beginning of Appendix A
for x=1:5
    X{x} = 1/sqrt(2)*[kron(A{x},Id)*psi; -kron(A{x},Id)*psi]; 
end
for y=1:4
    for s=1:2
        Y{y,s} = 1/sqrt(2)*[kron(Id,B{y,s})*psi; kron(Id,B{y,-s+3})*psi];
    end
end
I = 1/sqrt(2)*[psi; psi]; % unit vector

% Let us do some consistency checks:
% 1) We verify that the norm of the vectors is <=1:
norms=[];
for x=1:5
    norms=[norms norm(X{x})];
end
for y=1:4
    for s=1:2
        norms=[norms norm(Y{y,s})];
    end
end
norms=[norms norm(I)],

% 2) We verify that Alice's vectors X are orthogonal to the idendity vector I
orth=[];
for x=1:5
    orth=[orth dot(X{x},I)]; 
end
orth,

% 3) We verify that this vector model recovers the same correlations E{x,y} as the qm model
Ev=zeros(5,4);
for x=1:5
    for y=1:4
        Ev(x,y) = dot((I+X{x})/2,Y{y,1})+dot((I-X{x})/2,Y{y,2});
    end
end
E-Ev,

% The correlations <B(x,y,a)> given in the second line of (A9) are
EB=zeros(5,4,2);
for x=1:5
    for y=1:4
        for a=1:2
            EB(x,y,a) = dot(I+(-2*a+3)*X{x},Y{y,a});
        end
    end
end

% As in the last line of (A9), let us write Y as c I + d Yperp, where the
% Yperp are normalized
for y = 1:4
    for a = 1:2
        c{y,a} = dot(I,Y{y,a});
        Yperp{y,a} = Y{y,a}-c{y,a}*I;
        Yperp{y,a} = Yperp{y,a}/norm(Yperp{y,a});
        d{y,a} =  dot(Y{y,a},Yperp{y,a});
    end
end

% we can verify the last line of (A9)
EB2=zeros(5,4,2);
for x=1:5
    for y=1:4
        for a=1:2
            EB2(x,y,a) = c{y,a}+d{y,a}*dot((-2*a+3)*X{x},Yperp{y,a});
        end
    end
end
EB2-EB,

% However, we have that 1-d{y,a}/(1-|c{y,a}|) is negative for most y,a
pr=[];
for y=1:4
    for a=1:2
        pr=[pr 1-d{y,a}/(1-abs(c{y,a}))];
    end
end
pr,

