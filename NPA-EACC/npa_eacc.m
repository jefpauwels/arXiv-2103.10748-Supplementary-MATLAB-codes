clear all
clc

%%%%%%%%  NPA scenario for EACC with unrestricted entanglement %%%%%%
% Requires YALMIP
% Auxiliary files: "listflip.mat", "symmetryEACC.mat"

% Note that the file "symmetryEACC" must be adapted to the level of the
% hierarchy chosen in this file. Otherwise there is an error message. Up to
% level 2, this is just a matter of inserting/removing commented sections

%% Input the scenario and witness

nX=5; % #inputs Alice
nYa=4; % #inputs Bob
nB=2; % #outputs Bob
nM=3; % message dimension



nY=nM*nYa; % #effective inputs  Bob
nA=nM; % #effective outputs Alice



%%%% Provide witness coefficients coef(b,x,y) %%%%
coef=zeros(nB,nX,nYa);


% The 2bit-RAC witness
% T=[1 1;1 -1;-1 1;-1 -1];
% for x=1:nX
%     X=de2bi(x-1,2)+1;
%     for y=1:nYa
%         for b=1:nB
%             coef(b,x,y)=(-1)^(b-1)*T(x,y);
%         end
%     end
% end


% % The F-RAC witness
% beta=4;
% T=[1 1 beta;1 -1 beta;-1 1 beta;-1 -1 beta; 0 0 -4*beta];
% for x=1:nX
%     X=de2bi(x-1,2)+1;
%     for y=1:nYa
%         for b=1:nB
%             coef(b,x,y)=(-1)^(b-1)*T(x,y);
%         end
%     end
% end


% % The 5-Gallego witness 
% T=[1 1 1 1;1 1 1 -1;1 1 -1 0;1 -1 0 0;-1 0 0 0];
% for x=1:nX
%     for y=1:nYa
%         for b=1:nB
%             coef(b,x,y)=(-1)^(b-1)*T(x,y);
%         end
%     end
% end




%% Input the monomial list - the list can be extended indefinitly


S{1}=[];
cc=1;
for x=1:nX % A
    for a=1:nA
        cc=cc+1;
        S{cc}=[1 x a];
        CA(a,x)=cc;
    end
end


for y=1:nY % B
    for b=1:nB
        cc=cc+1;
        S{cc}=[2 y b];
        CB(b,y)=cc;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for x=1:nX % AB
    for a=1:nA
        for y=1:nY
            for b=1:nB
                cc=cc+1;
                S{cc}=[1 x a; 2 y b];
            end
        end
    end
end

for x=1:nX % AA
    for a=1:nA
        for xx=1:nX
            for aa=1:nA
                cc=cc+1;
                S{cc}=[1 x a; 1 xx aa];
            end
        end
    end
end


for y=1:nY % BB
    for b=1:nB
        for yy=1:nY
            for bb=1:nB
                cc=cc+1;
                S{cc}=[2 y b; 2 yy bb];
            end
        end
    end
end



% for x=1:nX % AAA
%     for a=1:nA
%         for xx=1:nX
%             for aa=1:nA
%                 for xxx=1:nX
%                     for aaa=1:nA
%                         cc=cc+1;
%                         S{cc}=[1 x a; 1 xx aa; 1 xxx aaa];
%                     end
%                 end
%             end
%         end
%     end
% end




num=cc;
sizeOfMomentMatrix=num


%% Build moment matrix


Gamma{1,1}=1;
for kk=2:num^2
    Tot=1;
    ii=de2bi(kk-1,2,num)+1; i=ii(1); j=ii(2);
    
    element=sortrows([S{i}; listflip(S{j})],1); %sort on first column
    pos=find(element(:,1)==1);
    AlicePart=element(pos,2:3);
    pos=find(element(:,1)==2);
    BobPart=element(pos,2:3);

    AliceSingle=AlicePart;
    inds=[]; % Send Alice's projectors squared to projectors
    for k=1:size(AliceSingle,1)-1
        if AliceSingle(k,:)==AliceSingle(k+1,:)
            inds=[inds k+1];
        end
    end
    AliceSingle(inds,:)=[];
    
    BobSingle=BobPart;
    inds=[]; % Send Bob's projectors squared to projectors
    for k=1:size(BobSingle,1)-1
        if BobSingle(k,:)==BobSingle(k+1,:)
            inds=[inds k+1];
        end
    end
    BobSingle(inds,:)=[];
    
    
    % Kill the moment if two of Alice's orthogonal projectors meet.
    for k=1:size(AliceSingle,1)-1
        if AliceSingle(k,1)==AlicePart(k+1,1) && AliceSingle(k,2)~=AlicePart(k+1,2)
            Tot=0;
        end
    end
    Alice=[ones(size(AliceSingle,1),1) AliceSingle];
    
    % Kill the moment if two of Bob's orthogonal projectors meet.
    for k=1:size(BobSingle,1)-1
        if BobSingle(k,1)==BobSingle(k+1,1) && BobSingle(k,2)~=BobSingle(k+1,2)
            Tot=0;
        end
    end
    Bob=[ones(size(BobSingle,1),1)+1 BobSingle];

    

    if Tot==0
        Gamma{i,j}=0;
    else
        Gamma{i,j}=[Alice;Bob];
    end
    
end

disp('Moment matrix lists complete...')




%% Build index matrix

t=0;
for k=1:num
    r=size(S{k},1)*size(S{k},2); if r>t; t=r; end;
end
maxString=2*t; % longest product in the moment matrix

tick=0;
for i=1:num
    for j=1:num
        tick=tick+1;
        s=size(Gamma{i,j});
        element=reshape(Gamma{i,j},[1 s(1)*s(2)]);
        element=[element zeros(1,maxString-s(1)*s(2))];
        GammaList(tick,:)=element;
    end
end
[~,~,reduced]=unique(GammaList,'rows');
indexMatrix=reshape(reduced,[num num]);
indexMatrix=reshape(min([indexMatrix(:), reshape(indexMatrix',num^2,1)],[],2),num,num);
[~,~,c]=unique(indexMatrix,'stable');
indexMatrix=reshape(c,[num num]);

disp(['The matrix currently has ', num2str(length(unique(indexMatrix))), ' variables']);

%% Standard symmetrisation for variable elimination
% One may always symmetrise over the classical communication since it does
% not appear explicitely in the witness. The symmetry group is the
% permutation group of [1, ..., nM].


% Generators
gens=[nM [1:nM-1]; % cycle 
     2 1 [3:nM]]; % swap

 
for K=1:2 % Apply both generators
generator=gens(K,:);
monomialcyclex0=symmetryEacc(generator, nX, nYa, nB, nM);
indexMatrixSwap=indexMatrix(monomialcyclex0,monomialcyclex0);
tick=0;
clear t
for i=1:num
    for j=1:num
        tick=tick+1;
        t(tick,:)=sort([indexMatrix(i,j) indexMatrixSwap(i,j)]);
    end
end
uniqueCrossList=unique(t,'rows');
for k=1:size(uniqueCrossList,1)
    indexMatrix(find(indexMatrix == uniqueCrossList(k,2))) = uniqueCrossList(k,1);
end


[~,~,c]=unique(indexMatrix);
indexMatrix=reshape(c,[num num]);
end
disp(['Finished basic symmetrisation ...'])
disp(['The matrix currently has ', num2str(length(unique(indexMatrix))), ' variables']);



%% Build SDP

num=size(indexMatrix,1);

%%%%%%%%%%%%%% The SDP %%%%%%%%%%%%%%%%%
clear M v
v = sdpvar(max(max(indexMatrix)), 1);

%%%% Moment matrix %%%%


for i=1:max(max(indexMatrix))
    pos=find(i==indexMatrix);
    M(pos)=v(i);
end
id=indexMatrix(1,1);
pos=find(id==indexMatrix);
M(pos)=1; % normalisation


Zero=indexMatrix(CA(1,1),CA(2,1));
pos=find(Zero==indexMatrix);
M(pos)=0; % orthogonalities
M=reshape(M,[num num]);

F=[M>=0];
%%%% probabilities in PM scenario %%%
prob=cell(nB,nX,nYa);
obj=0;
for x=1:nX
    for y=1:nYa
        for b=1:nB
            st=0;
            for a=1:nA
                Y=fromSeveralBases([y,a]-1,[nYa nA])+1;
                st=st+M(CA(a,x),CB(b,Y));
                F=[F];
            end
            prob{b,x,y}=st;
        end
    end
end


%%%%% objective %%%%
obj=0;
for x=1:nX
    for y=1:nYa
        s=0;
        for b=1:nB
            s=s+prob{b,x,y};
            obj=obj+coef(b,x,y)*prob{b,x,y};
        end
        F=[F,s==1];
    end
end




disp('Options')
ops=sdpsettings('solver','mosek', 'cachesolvers', 1);
disp('Solving')
diagnostic=solvesdp(F,-obj,ops)
obj=double(obj)
M=double(M);
