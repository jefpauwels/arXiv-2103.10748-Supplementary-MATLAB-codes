function out=listPadding(list)

% Add zeros to a list so that all moments are of the same size
global maxString
 
L=length(list);

padding=[list zeros(1,maxString-L)];

out=padding;

end

