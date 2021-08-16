function out=listReduce(list)

% Take the moment and apply the quantum reduction rules (projective measurements, cyclicity and tr(AB)=tr(A'B') in real space)
% Order the final output list systematically.
global lr lm  nOp  Vac 


%% Remove identities
posIdentity=find(1==list); % find all the identities...
if isequal(unique(list),1)==1
    list=[1]; 
else
    list(posIdentity)=[]; % ... and remove them
end




%% sort list such that the first and last operators are disconnected (if possible)

L=length(list);
for k=0:L-1
    cyc=circshift(list,k); % check all cyclic shifts of the moment
    if linkedOperators(cyc(1),cyc(end))==0
        list=cyc;
        break;
    end
end

%% Find elements to be deleted from the moment

kill=0; % =1 if the moment equals zeros
delete=[];
listHold=list;
listHold(delete)=[];
L=length(listHold);
delete=[];
listHold(delete)=[];


L=length(listHold);
delete=[];
for k=1:L-1
    %%%%%%% measurement projectivity %%%%%
    el=listHold(k);
    if ismember(el,lm)==1
        [p1,p2]=find(el==lm);
        elNext=listHold(k+1); [q1,q2]=find(elNext==lm);
        if el==elNext
            delete=[delete k+1]; % M^2=M
        elseif p1==q1 % same measurement, different outcome
            kill=1; % kill the moment
        end
    end
end
listHold(delete)=[];
list=listHold;

%%%%%%%%%% Reduce moment further %%%%%%%%%%%%%%
if kill==1
    list=0; % this means moment equals 0
end

if length(list)==1 && ismember(list,lr)==1 % just a state
    list=lr(1); % assign it to the first state
end


%% Order the output list after smallest element - using the tr(AB)=tr(A'B') relation 

L=length(list);
if L>1
    % find the smallest list element
    minimal=min(list);
    posMinimal=find(minimal==list); % all positions in which the smallest value appears
    val=(nOp)^length(list); % number larger than all possible list values
    for k=1:length(posMinimal)
        candidate=circshift(list,1-posMinimal(k));% order such that the smallest value comes first
        ll=length(candidate); res=fromSeveralBases(candidate-1,nOp*ones(1,ll));
        if res<val % if the list value is smaller ...
            val=res;
            winner=candidate; % ... pick out the list
        end
        % All daggers
        for kk=1:L-1
            candidate2=[flip(candidate(1:kk)) flip(candidate(kk+1:end))];
            ll2=length(candidate2);
            res=fromSeveralBases(candidate2-1,nOp*ones(1,ll2));
            if res<val  % if the list value is smaller ...
                val=res;
                winner=candidate2; % ... pick out the list
            end
        end
    end
    list=winner;
end

out=list;


end

