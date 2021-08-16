function out=listflip(M)

s=size(M,1);
p=flip([1:s]);

out=M(p,:);
end
