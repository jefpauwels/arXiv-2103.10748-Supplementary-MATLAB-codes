function out=linkedOperators(A,B)

% out=0 : A and B do not influence each other through the quantum reduction rules
% out=1 : A and B influence each other through the quantum reduction rules

global lr lm Vac


out=0;

if ismember(A,lm)==1 && ismember(B,lm)==1
    [p1,p2]=find(A==lm); [q1,q2]=find(B==lm);
    if p1==q1
        out=1; % measurement operators in the same basis
    end
end

end
