function out=symmetryRAC(omega,type)
% Action of permutation omega on the list of monomials. 
% type=1 means a permutation of positions [x1 x2 x3 ...] and y
% type=2 means a permutation of x1 and b


global nX nY nB lr Vac lm

lomega=length(omega); omega2=[];
for k=1:lomega
    omega2=[omega2 find(k==omega)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Type 1  symmetry %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if type==1

cc=0;
t(1)=1; % id stays
t(2)=2; % Vac stays
cc=2;
%%%%%%%%%%%% Level 1 %%%%%%%%%%%%%%%%%%
tick=cc;
for x=1:nX % rho
    cc=cc+1; X=toSeveralBases(x-1,nB*ones(1,nY))+1;
    X=X(omega); X=fromSeveralBases(X-1,nB*ones(1,nY))+1;
    t(cc)=tick+X;   
end
tick=cc;
for y=1:nY % M
    for b=1:nB
        cc=cc+1; X=fromSeveralBases([omega2(y) b]-1,[nY nB])+1;
        t(cc)=tick+X;
    end
end
%%%%%%%%%%% Level 2 %%%%%%%%%%%%%%%%%%%%
tick=cc;
for x=1:nX % rho*M
    for y=1:nY
        for b=1:nB
        cc=cc+1; X=toSeveralBases(x-1,nB*ones(1,nY))+1; 
        X=X(omega); pos=fromSeveralBases([X omega2(y) b]-1,[nB*ones(1,nY) nY nB])+1;
        t(cc)=tick+pos;
        end
    end
end
tick=cc;
for x=1:nX % rho*rho
    for xx=1:nX
    cc=cc+1; X=toSeveralBases(x-1,nB*ones(1,nY))+1; XX=toSeveralBases(xx-1,nB*ones(1,nY))+1;
    X=X(omega); XX=XX(omega); pos=fromSeveralBases([X XX]-1,nB*ones(1,2*nY))+1; 
    t(cc)=tick+pos;
    end
end
% tick=cc;
% for y=1:nY % M*M
%     for b=1:nB
%         for yy=1:nY
%             for bb=1:nB
%                 cc=cc+1; pos=fromSeveralBases([omega2(y) b omega2(yy) bb]-1,[nY nB nY nB])+1;
%                 t(cc)=tick+pos;
%             end
%         end
%     end
% end
tick=cc;
for x=1:nX % Vac*rho
    cc=cc+1; X=toSeveralBases(x-1,nB*ones(1,nY))+1;
    X=X(omega); X=fromSeveralBases(X-1,nB*ones(1,nY))+1;
    t(cc)=tick+X;   
end
tick=cc;
for y=1:nY % Vac*M
    for b=1:nB
        cc=cc+1; X=fromSeveralBases([omega2(y) b]-1,[nY nB])+1;
        t(cc)=tick+X;
    end
end
cc=cc+1;
t(cc)=cc; %Vac*Vac
tick=cc;
% for x=1:nX %rho*rho*rho
%     for xx=1:nX
%         for xxx=1:nX
%     cc=cc+1; X=toSeveralBases(x-1,nB*ones(1,nY))+1; XX=toSeveralBases(xx-1,nB*ones(1,nY))+1;
%     XXX=toSeveralBases(xxx-1,nB*ones(1,nY))+1;
%     X=X(omega); XX=XX(omega); XXX=XXX(omega);
%     pos=fromSeveralBases([X XX XXX]-1,nB*ones(1,3*nY))+1; 
%     t(cc)=tick+pos;
%         end
%     end
% end


end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Type 2  symmetry %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if type==2
   
cc=0;
t(1)=1; % id stays
t(2)=2; % Vac stays
cc=2;

%%%%%%%%%%%% Level 1 %%%%%%%%%%%%%%%%%%
tick=cc;
for x=1:nX % rho
    cc=cc+1; X=toSeveralBases(x-1,nB*ones(1,nY))+1;
    X(1)=omega(X(1)); X=fromSeveralBases(X-1,nB*ones(1,nY))+1;
    t(cc)=tick+X;   
end
tick=cc;
for y=1:nY % M
    for b=1:nB
        cc=cc+1; X=fromSeveralBases([y (omega(b)*[y==1]+b*[y~=1])]-1,[nY nB])+1;
        t(cc)=tick+X;
    end
end
%%%%%%%%%%% Level 2 %%%%%%%%%%%%%%%%%%%%
tick=cc;
for x=1:nX % rho*M
    for y=1:nY
        for b=1:nB
        cc=cc+1; X=toSeveralBases(x-1,nB*ones(1,nY))+1; 
        X(1)=omega(X(1)); pos=fromSeveralBases([X y (omega(b)*[y==1]+b*[y~=1])]-1,[nB*ones(1,nY) nY nB])+1;
        t(cc)=tick+pos;
        end
    end
end
tick=cc;
for x=1:nX % rho*rho
    for xx=1:nX
    cc=cc+1; X=toSeveralBases(x-1,nB*ones(1,nY))+1; XX=toSeveralBases(xx-1,nB*ones(1,nY))+1;
    X(1)=omega(X(1)); XX(1)=omega(XX(1)); pos=fromSeveralBases([X XX]-1,nB*ones(1,2*nY))+1; 
    t(cc)=tick+pos;
    end
end
tick=cc;
% for y=1:nY % M*M
%     for b=1:nB
%         for yy=1:nY
%             for bb=1:nB
%                 cc=cc+1; pos=fromSeveralBases([y (omega(b)*[y==1]+b*[y~=1]) yy (omega(bb)*[yy==1]+bb*[yy~=1])]-1,[nY nB nY nB])+1;
%                 t(cc)=tick+pos;
%             end
%         end
%     end
% end
tick=cc;
for x=1:nX % Vac*rho
    cc=cc+1; X=toSeveralBases(x-1,nB*ones(1,nY))+1;
    X(1)=omega(X(1)); pos=fromSeveralBases(X-1,nB*ones(1,nY))+1;
    t(cc)=tick+pos;   
end
tick=cc;
for y=1:nY % Vac*M
    for b=1:nB
        cc=cc+1; pos=fromSeveralBases([y (omega(b)*[y==1]+b*[y~=1])]-1,[nY nB])+1;
        t(cc)=tick+pos;
    end
end
cc=cc+1;
t(cc)=cc; %Vac*Vac
tick=cc;
% for x=1:nX %rho*rho*rho
%     for xx=1:nX
%         for xxx=1:nX
%     cc=cc+1; X=toSeveralBases(x-1,nB*ones(1,nY))+1; XX=toSeveralBases(xx-1,nB*ones(1,nY))+1;
%     XXX=toSeveralBases(xxx-1,nB*ones(1,nY))+1;
%     X(1)=omega(X(1)); XX(1)=omega(XX(1));  XXX(1)=omega(XXX(1));pos=fromSeveralBases([X XX XXX]-1,nB*ones(1,3*nY))+1; 
%     t(cc)=tick+pos;
%         end
%     end
% end


end















out=t;
end