function out=findPosition(elem)

% A row and column index in the moment matrix whose product give elem

global S num

nE=length(elem);

pos=[];
if isempty(elem)==1
    pos=[1 1];
end

for r=1:nE
    elem=[elem(end) elem(1:end-1)];
    for t=1:nE % Loop over all ways of dividing elem into two parts
        left=elem(1:t); right=elem(t+1:end);
        pleft=findMonomial(left); pright=findMonomial(flip(right));
        if isempty(pleft)==0 && isempty(pright)==0
            pos=[pleft, pright];
            break;
        end
    end
end

out=pos;
end