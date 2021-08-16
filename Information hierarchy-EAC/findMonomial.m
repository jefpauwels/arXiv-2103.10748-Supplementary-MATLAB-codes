function out=findMonomial(elem)

% Find the position of elem in the monomial list S

global S num

pos=[];
if isempty(elem)==1
    pos=1;
else
    for k=1:num
        ans=isequal(elem,S{k});
        if ans==1;
            pos=k;
            break;
        end
    end
end

out=pos;

end

