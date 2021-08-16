function out=executeHierarchy

% This file executes the hierarchy problem provided in "launch_hierarchy".
% No modifications need to be made here for a standard implementation.

global lr lm nX nY nB nOp Vac S L num coef  R2 R3 num2
global prior GuessingProb maxString  
t=0;
for k=1:num
    r=length(S{k}); if r>t; t=r; end;
end
maxString=2*t; % longest product in the moment matrix


%% Make primary index matrix

GammaRow=zeros(num^2,maxString);
for i1=1:num
    for i2=i1:num
        ind=fromSeveralBases([i1 i2]-1,[num num])+1;
        left=S{i1}; right=flip(S{i2});
        prodRaw=[left, right]; % raw multiplication
        prodReduced=listReduce(prodRaw); % reduction rules
        prodPadded=listPadding(prodReduced); % padd out with 0s
        GammaRow(ind,:)=prodPadded; % put the moment as a row in a matrix
    end
end


[~,~,reduced]=unique(GammaRow,'rows'); % find the indices for identical rows
indexMatrixHalf=reshape(reduced,[num num]); % make a matrix
for i1=1:num
    for i2=i1:num
        indexMatrix(i1,i2)=indexMatrixHalf(i2,i1);% Extend the upper triangle to the lower triangle
        indexMatrix(i2,i1)=indexMatrixHalf(i2,i1);% Extend the upper triangle to the lower triangle
    end
end

clear indexMatrixHalf GammaRow % Kill this to save memory
indexMatrix=reshape(min([indexMatrix(:), reshape(indexMatrix',num^2,1)],[],2),num,num);
[~,~,c]=unique(indexMatrix,'stable');
indexMatrix=reshape(c,[num num]); % indexMatrix with variables ordered from 1,2,...


nVar=length(unique(indexMatrix)); % number of different elements in moment matrix (roughly the number of SDP variables)
disp(['The index matrix has ', num2str(nVar), ' variables']);


%% localising index matrix (positivity) [tr(L*(rho_x)*L')]


clear matrix
for x=1:nX
    matrix=zeros(num2,num2);
    for i1=1:num2
        for i2=i1:num2
            left=L{i1}; right=[lr(x), flip(L{i2})];
            prodRaw=[left, right]; % raw multiplication
            prodReduced=listReduce(prodRaw); % quantum reduction rules
            if prodReduced==0
                pos=[lm(1,1), lm(1,2)]; % vanishing moments are equated with M1*M2=0
            else
                pos=findPosition(prodReduced); % find position of the moment in primary moment matrix
            end
            matrix(i1,i2)=indexMatrix(pos(1),pos(2));  % Extend the upper triangle to the lower triangle
            matrix(i2,i1)=indexMatrix(pos(1),pos(2));  % Extend the upper triangle to the lower triangle
        end
    end
    indexMatrixPositivity{x}=matrix; % indexMatrix for positivity - for each state x
end
%disp('Localising matrices (Positivity) complete...')


%% localising index matrix (guessing) [tr(L*(sig)*L')]

clear matrix
matrix=zeros(num2,num2);
for i1=1:num2
    for i2=i1:num2
        left=L{i1}; right=[Vac, flip(L{i2})];
        prodRaw=[left, right]; % raw multiplication
        prodReduced=listReduce(prodRaw); % quantum reduction rules
        if prodReduced==0
            pos=[lm(1,1), lm(1,2)]; % vanishing moments are equated with M1*M2=0
        else
            pos=findPosition(prodReduced); % find position of the moment in primary moment matrix
        end
        matrix(i1,i2)=indexMatrix(pos(1),pos(2));  % Extend the upper triangle to the lower triangle
        matrix(i2,i1)=indexMatrix(pos(1),pos(2));  % Extend the upper triangle to the lower triangle
    end
end
indexMatrixSigma=matrix; % indexMatrix for guessing probability
%disp('Localising matrices (Guessing) complete...')

disp('All localising matrices complete...')


%% SDP part

clear Gamma GammaX GammaSig

cvx_begin sdp
variable v(nVar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PRIMARY MOMENT MATRIX %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nVar
    pos=find(i==indexMatrix); % find all positions with index i
    Gamma(pos)=v(i); % assign them the SDP variable v(i)
end
Gamma=reshape(Gamma,[num num]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% CONSTRAINTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gamma(1,2)<=GuessingProb;   % Information restriction
Gamma(1,lr(1))==1;              % trace(state)=1
if nB>2;
Gamma(lm(1,1),lm(1,2))==0;      % trace(M*M')=0
end
Gamma>=0;                       % positivity of moment matrix



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% LOCALISING MOMENT MATRIX (POSITIVITY) %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for x=1:nX
    ind=indexMatrixPositivity{x};
    clear matrix
    for k1=1:num2
        for k2=1:num2
            matrix(k1,k2)=v(ind(k1,k2));
        end
    end
    matrix>=0; % positivity of localising matrix
    GammaX{x}=matrix; % tr(L*(rho_x)*L')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% LOCALISING MOMENT MATRIX (Guessing) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for x=1:nX
    ind=indexMatrixPositivity{x};
    clear matrix
    for k1=1:num2
        for k2=1:num2
            matrix(k1,k2)=v(indexMatrixSigma(k1,k2))-prior(x)*v(ind(k1,k2));
        end
    end
    matrix>=0; % positivity of localising matrix
    GammaSig{x}=matrix; % tr(L*(sig-q(x)*rho_x)*L')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% OBJECTIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj=0;
for x=1:nX
    for y=1:nY
        for b=1:nB-1
        obj=obj+coef(b,x,y)*Gamma(lr(x),lm(y,b)); % The witness
        end
    end
end
for x=1:nX
    for y=1:nY
        s=0;
        for b=1:nB-1
            s=s+Gamma(lr(x),lm(y,b)); % The witness
        end
        obj=obj+coef(nB,x,y)*(1-s);
    end
end
maximise(obj)

cvx_end

out{1}=cvx_optval;
out{2}=Gamma;
out{3}=GammaX
out{4}=GammaSig;


end
