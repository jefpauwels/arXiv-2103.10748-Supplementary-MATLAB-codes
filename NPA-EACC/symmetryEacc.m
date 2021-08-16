function out=symmetryEacc(omega, nX, nYa, nB, nM)

% Outputs a reordered monomial list based on the permutation omega.
% For each piece of the monomial list, we implement the permutation
% separately. Therefore the pieces of the monomial list permuted here must match
% the original monomial list in the main file.


clear action out
nY=nYa*nM;

cc=1;
action(1)=1;

L=cc;
for x=1:nX % A
    for a=1:nM
        cc=cc+1;
        pos=fromSeveralBases([x omega(a)]-1,[nX nM])+1;
        action(cc)=L+pos;
    end
end


L=cc;
for y=1:nY % B
    for b=1:nB
        cc=cc+1;
        Y=toSeveralBases(y-1,[nYa, nM])+1;
        pos=fromSeveralBases([Y(1) omega(Y(2)) b]-1,[nYa, nM, nB])+1;
        action(cc)=L+pos;
    end
end



L=cc; % AB
for x=1:nX
    for a=1:nM
        for y=1:nY
            for b=1:nB
                cc=cc+1;
                Y=toSeveralBases(y-1,[nYa, nM])+1;
                pos=fromSeveralBases([x omega(a) Y(1) omega(Y(2)) b]-1,[nX nM, nYa, nM, nB])+1;
                action(cc)=L+pos;
            end
        end
    end
end



L=cc; % AA
for x=1:nX
    for a=1:nM
        for xx=1:nX
            for aa=1:nM
                cc=cc+1;
                pos=fromSeveralBases([x omega(a) xx omega(aa)]-1,[nX nM, nX nM])+1;
                action(cc)=L+pos;
            end
        end
    end
end


L=cc; % BB
for y=1:nY
    for b=1:nB
        for yy=1:nY
            for bb=1:nB
                cc=cc+1;
                Y=toSeveralBases(y-1,[nYa, nM])+1;
                YY=toSeveralBases(yy-1,[nYa, nM])+1;
                pos=fromSeveralBases([Y(1) omega(Y(2)) b YY(1) omega(YY(2)) bb]-1,[nYa, nM, nB, nYa, nM, nB])+1;
                action(cc)=L+pos;
            end
        end
    end
end

% 
% L=cc; % AAA
% for x=1:nX
%     for a=1:nM
%         for xx=1:nX
%             for aa=1:nM
%                 for xxx=1:nX
%                     for aaa=1:nM
%                         cc=cc+1;
%                         pos=fromSeveralBases([x omega(a) xx omega(aa) xxx omega(aaa)]-1,[nX nM, nX nM,nX nM])+1;
%                         action(cc)=L+pos;
%                     end
%                 end
%             end
%         end
%     end
% end


out=action;
end
