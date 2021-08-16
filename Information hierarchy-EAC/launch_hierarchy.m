clear all
clc

% SDP hierarchy for informationally restricted correlations. Launched from
% this file exclusively. Input the scenario, the witness and the monomial
% list (relaxation level).

% In this code, the monomial list for the localising matrices for rho_x>=0
% and the localising matrices for sigma>=p_x rho_x are identical.

% If the level of the main moment matrix is insufficient to account for all moments appearing 
% in the localising moment matrices, then the code crashes. This means that
% the main monomial list must be extended apropriately.


% Requires YALMIP

% Requires axuiliary files: "fromSeveralBases.mat", "listReduce.mat",
% "listPadding.mat", "findPosition.mat", "linkedOperators.mat"



global lr lm nX nY nB nOp Vac S L num coef  R2 R3 num2 
global prior GuessingProb maxString 


%% Choose the scenario 


nX=5; % number of preparation inputs
nY=4; % number of measurement inputs
nB=2; % number of measurement outcomes.

%% Choose the witness
%%%%%%%%%%% 2bit-RAC %%%%%%%%%%%%
% For EACC with d=2, use local level 1 and main level {level 1, rho*M, rho*rho} for a
% tight bound.
% For EACC with d=3, use local level {level 1, M*M, M*M*M} and main level
% {level 1, rho*M, rho*rho, M*M, Vac*rho, Vac*M, Vac*Vac, M*M*Vac, M*M*M,
% M*rho*M, M*M*M*M} for a tight bound.
% T=[1 1;1 -1;-1 1;-1 -1];
% for x=1:nX
%     for y=1:nY
%         for b=1:nB
%             coef(b,x,y)=(-1)^(b-1)*T(x,y);
%         end
%     end
% end

%%%%%%%%%%% F-RAC %%%%%%%%%%%%
% For EAQC with d=2, use local level {level 1, M*M, M*M*M, M*M*M*M} and main level
% {level 1, rho*M, rho*rho, M*M, Vac*rho, Vac*M, Vac*Vac, M*M*Vac, M*M*M, M*rho*M,
% M*M*M*M, M*M*M*M*Vac, rho*M*M*M, rho*M*M*M*M}  for a tight bound.
% beta=4;
% T=[1 1 beta;1 -1 beta;-1 1 beta;-1 -1 beta;0 0 -4*beta];
% for x=1:nX
%     for y=1:nY
%         for b=1:nB
%             coef(b,x,y)=(-1)^(b-1)*T(x,y);
%         end
%     end
% end


%%%%%%%%%%% Gallego witness %%%%%%%%%%%%
% For EACC with d=2,3,4, use local level {level 1, M*M, M*M*M} and main level
% {level 1, rho*M, rho*rho, M*M, Vac*rho, Vac*M, Vac*Vac, M*M*Vac, M*M*M, M*rho*M,
% M*M*M*M}  for the bound reported in the paper.
T=[1 1 1 1;1 1 1 -1;1 1 -1 0;1 -1 0 0;-1 0 0 0];
for x=1:nX
    for y=1:nY
        for b=1:nB
            coef(b,x,y)=(-1)^(b-1)*T(x,y);
        end
    end
end





%% Choose information constraint 

prior=1/nX*ones(1,nX);
GuessingProb=4/nX;


%% Choose the relaxation level (primary moment matrix & localising moment matrix)

% {id, Vac, rho, M}
nOp=1+1+nX+nY*(nB-1);
lOp=[1:nOp]; Vac=2; lr=[3:2+nX]; lm=reshape([3+nX:nOp],[nB-1 nY])'; %(y,b) indexed
%%%%%%%%%% Level 1 %%%%%%%%%%%%%
cc=0;
for k=1:nOp
    cc=cc+1;
    S{k}=lOp(k);
end

%%%%%%%%%% Level 2 %%%%%%%%%%%%%
for x=1:nX % rho*M
    for y=1:nY
        for b=1:nB-1
        cc=cc+1;
        S{cc}=[lr(x) lm(y,b)];
        end
    end
end

% for y=1:nY % M*rho
%     for b=1:nB-1
%         for x=1:nX
%         cc=cc+1;
%         S{cc}=[lm(y,b) lr(x)];
%         end
%     end
% end

for x=1:nX % rho*rho
    for xx=1:nX
        cc=cc+1;
        R2(x,xx)=cc;
        S{cc}=[lr(x) lr(xx)];
    end
end
for y=1:nY % M*M
    for b=1:nB-1
        for yy=1:nY
            for bb=1:nB-1
                cc=cc+1;
                S{cc}=[lm(y,b) lm(yy,bb)];
            end
        end
    end
end
for x=1:nX % Vac*rho
    cc=cc+1;
    S{cc}=[Vac lr(x)];
end
for y=1:nY % Vac*M
    for b=1:nB-1
    cc=cc+1;
    S{cc}=[Vac lm(y,b)];
    end
end
cc=cc+1;
S{cc}=[Vac,Vac];

%%%%%%%%%% Level 3 %%%%%%%%%%%%%
% for x=1:nX % Vac*rho*Vac
%     cc=cc+1;
%     S{cc}=[Vac lr(x) Vac];
% end
for y=1:nY % M*M*Vac
    for b=1:nB-1
        for yy=1:nY
            for bb=1:nB-1
                    cc=cc+1;
                    S{cc}=[lm(y,b) lm(yy,bb) Vac];
            end
        end
    end
end
% for x=1:nX % rho*M*Vac
%     for y=1:nY
%         for b=1:nB-1
%             cc=cc+1;
%             S{cc}=[lr(x), lm(y,b), Vac];
%         end
%     end
% end
% for y=1:nY % rho*M*M
%     for b=1:nB-1
%         for yy=1:nY
%             for bb=1:nB-1
%                 for x=1:nX
%                     cc=cc+1;
%                     S{cc}=[lr(x) lm(y,b) lm(yy,bb) ];
%                 end
%             end
%         end
%     end
% end
for y=1:nY %M*M*M
    for b=1:nB-1
        for yy=1:nY
            for bb=1:nB-1
                for yyy=1:nY
                    for bbb=1:nB-1
                        cc=cc+1;
                        S{cc}=[lm(y,b) lm(yy,bb) lm(yyy,bbb)];
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for x=1:nX %rho*rho*rho
%     for xx=1:nX
%         for xxx=1:nX
%             cc=cc+1;
%             R3(x,xx,xxx)=cc;
%             S{cc}=[lr(x) lr(xx) lr(xxx)];
%         end
%     end
% end
for y=1:nY %M*rho*M
    for b=1:nB-1
        for x=1:nX
            for yy=1:nY
                for bb=1:nB-1
                    cc=cc+1;
                    S{cc}=[lm(y,b) lr(x) lm(yy,bb)];
                end
            end
        end
    end
end

% for x=1:nX % rho*rho*M
%     for xx=1:nX
%         for y=1:nY
%             for b=1:nB-1
%                 cc=cc+1;
%                 S{cc}=[lr(x) lr(xx) lm(y,b)];
%             end
%         end
%     end
% end
% for y=1:nY %M*rho*rho
%     for b=1:nB-1
%         for x=1:nX
%             for xx=1:nX
%                 cc=cc+1;
%                 S{cc}=[lm(y,b) lr(x) lr(xx)];
%             end
%         end
%     end
% end

%%%%%%%%%%%%%% All levels above this line are needed if Localising Matrices
%%%%%%%%%%%%%% involved M*M and rho*M 


% for x=1:nX % Vac*rho*rho
%     for xx=1:nX
%         cc=cc+1;
%         S{cc}=[Vac lr(x) lr(xx)];
%     end
% end
%%%%%%%%%%%%% Fourth level %%%%%%%%%%%%%%%%%%%%%
for y=1:nY %M*M*M*M
    for b=1:nB-1
        for yy=1:nY
            for bb=1:nB-1
                for yyy=1:nY
                    for bbb=1:nB-1
                        for yyyy=1:nY
                            for bbbb=1:nB-1
                        cc=cc+1;
                        S{cc}=[lm(y,b) lm(yy,bb) lm(yyy,bbb) lm(yyyy,bbbb)];
                            end
                        end
                    end
                end
            end
        end
    end
end
% for y=1:nY %M*M*M*M*Vac
%     for b=1:nB-1
%         for yy=1:nY
%             for bb=1:nB-1
%                 for yyy=1:nY
%                     for bbb=1:nB-1
%                         for yyyy=1:nY
%                             for bbbb=1:nB-1
%                         cc=cc+1;
%                         S{cc}=[lm(y,b) lm(yy,bb) lm(yyy,bbb) lm(yyyy,bbbb) Vac];
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% for x=1:nX %rho*M*M*M
%     for y=1:nY
%         for b=1:nB-1
%             for yy=1:nY
%                 for bb=1:nB-1
%                     for yyy=1:nY
%                         for bbb=1:nB-1
%                             cc=cc+1;
%                             S{cc}=[lr(x) lm(y,b) lm(yy,bb) lm(yyy,bbb)];
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% for x=1:nX %rho*M*M*M*M
%     for y=1:nY 
%         for b=1:nB-1
%             for yy=1:nY
%                 for bb=1:nB-1
%                     for yyy=1:nY
%                         for bbb=1:nB-1
%                             for yyyy=1:nY
%                                 for bbbb=1:nB-1
%                                     cc=cc+1;
%                                     S{cc}=[lr(x) lm(y,b) lm(yy,bb) lm(yyy,bbb) lm(yyyy,bbbb)];
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

num=cc % size of primary moment matrix




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Localising matrix %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Level 1 %%%%%%%%%%%%%
cc=0;
for k=1:nOp
    cc=cc+1;
    L{k}=lOp(k);
end

%%%%%%%%%% Level 2 %%%%%%%%%%%%%
% for x=1:nX %rho*rho
%     for xx=1:nX
%         cc=cc+1;
%         L{cc}=[lr(x) lr(xx)];
%     end
% end
for y=1:nY %M*M
    for b=1:nB-1
        for yy=1:nY
            for bb=1:nB-1
                cc=cc+1;
                L{cc}=[lm(y,b) lm(yy,bb)];
            end
        end
    end
end
% for x=1:nX %rho*M
%     for y=1:nY
%         for b=1:nB-1
%         cc=cc+1;
%         L{cc}=[lr(x) lm(y,b)];
%         end
%     end
% end
for y=1:nY %M*M*M
    for b=1:nB-1
        for yy=1:nY
            for bb=1:nB-1
                for yyy=1:nY
                    for bbb=1:nB-1
                        cc=cc+1;
                        L{cc}=[lm(y,b) lm(yy,bb) lm(yyy,bbb)];
                    end
                end
            end
        end
    end
end
% for y=1:nY %M*M*M*M
%     for b=1:nB-1
%         for yy=1:nY
%             for bb=1:nB-1
%                 for yyy=1:nY
%                     for bbb=1:nB-1
%                         for yyyy=1:nY
%                             for bbbb=1:nB-1
%                                 cc=cc+1;
%                                 L{cc}=[lm(y,b) lm(yy,bb) lm(yyy,bbb) lm(yyyy,bbbb)];
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
num2=cc % % size of localising (positivty) moment matrix 

%% Call the SDP

out=executeHierarchy; % The main file that runs the hierarchy
bound=out{1};  % the result
Gamma=out{2}; % the moment matrix
GammaX=out{3}; % localising matrix for positivity of states
GammaVac=out{4}; % localising matrix for guessing probability




