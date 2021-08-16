%%%%%%%%%%%%% This is a list of witnesses that can be used in the hierarchy
%%%%%%%%%%%%% packaged %%%%%%%%%%%%%%%%%%%

%%%%% RIC %%%%%%%%%
% coef=zeros(nB,nX,nY);
% for x=1:nX
%     for y=1:nY
%         if x==y
%             coef(1,x,y)=1/nX^2;
%         else
%             coef(2,x,y)=1/nX^2;
%         end
%     end
% end


%%%%%%%%%%%% RAC %%%%%%%%%%%%%%
% for x=1:nX
%     for y=1:nY
%         for b=1:nB
%     xx=toSeveralBases(x-1,nB*ones(1,nY))+1;
%     coef(b,x,y)=[b==xx(y)]/(nX*nY);
%         end
%     end
% end

%%%%%%%% Gallego %%%%%%%%%
% T=[-1 -1;-1 1;1 0];
% for x=1:nX
%     for y=1:nY
%         for b=1:nB
%             coef(b,x,y)=T(x,y)*(-1)^(b+1);
%         end
%     end
% end


%%%%%%% information game %%%%%%%%
% coef=zeros(nX,nX,1);
% q=1/nX*ones(1,nX);
% for x=1:nX
%     coef(x,x,1)=q(x);
% end


%%%%%%%%%%%%%% New witness for Concepcion experiment %%%%%%%%%%%%%%5
% c=zeros(6,10,2);
% c(1,1,1)=1; c(2,1,1)=-1;
% c(1,2,1)=1; c(2,2,1)=-1;
% c(2,3,1)=1; c(1,3,1)=-1;
% c(2,4,1)=1; c(1,4,1)=-1;
% c(3,5,1)=1; c(4,5,1)=-1;
% c(3,6,1)=1; c(4,6,1)=-1;
% c(4,7,1)=1; c(3,7,1)=-1;
% c(4,8,1)=1; c(3,8,1)=-1;
% c(5,9,1)=1; 
% c(6,10,1)=1; 
% 
% c(4,1,2)=1; c(3,1,2)=-1;
% c(3,2,2)=1; c(4,2,2)=-1;
% c(4,3,2)=1; c(3,3,2)=-1;
% c(3,4,2)=1; c(4,4,2)=-1;
% c(1,5,2)=1; c(2,5,2)=-1;
% c(2,6,2)=1; c(1,6,2)=-1;
% c(1,7,2)=1; c(2,7,2)=-1;
% c(2,8,2)=1; c(1,8,2)=-1;
% c(6,9,2)=1; 
% c(5,10,2)=1; 
% coef=c;


%%%%%%% State exclusion %%%%%%%%%
coef=zeros(nB,nX,nY);
for x=1:nX
    for b=1:nB
        if x==b
            coef(b,x,1)=-1/nX;
        end
    end
end



% % The Stefano game
% beta=4;
% gamma=4*beta;
% c=[1 1 beta;1 -1 beta; -1 1 beta; -1 -1 beta; 0 0 -gamma];
% for x=1:nX
%     for y=1:nY
%         for b=1:nB
%             coef(b,x,y)=(-1)^(b-1)*c(x,y);
%         end
%     end
% end



%%% Gallego
% c=[1 1 1 1; 1 1 1 -1; 1 1 -1 0; 1 -1 0 0;-1 0 0 0];
% for x=1:nX
%     for y=1:nY
%         for b=1:nB
%             coef(b,x,y)=(-1)^(b-1)*c(x,y);
%         end
%     end
% end


%%%%%%%%%%%%% Witness for STOCKHOLM experiment %%%%%%%%%%%

% c=zeros(2,5,6);
% c(1,1,1)=1;
% c(1,2,1)=-1;
% c(1,3,1)=-1;
% c(1,1,2)=1;
% c(1,4,2)=-1;
% c(1,1,3)=1;
% c(1,5,3)=-1;
% c(1,2,4)=1;
% c(1,3,4)=-1;
% c(1,4,4)=-1;
% c(1,5,4)=-1;
% c(1,3,5)=1;
% c(1,4,5)=-1;
% c(1,5,5)=-1;
% c(1,4,6)=1;
% c(1,5,6)=-1;
% coef=c;
