function is = toSeveralBases(n, bases)
% is = toSeveralBases(n, bases)
%
% Gives a list describing the number n, where each digit can have a
% different basis.
% 
% Exemple : toSeveralBases(21, [2 3 4]) gives [1 2 1] because 
%   21 = 1*(3*4) + 2*(4) + 1

% Written by Jean-Daniel Bancal, last modified 7 March 2011.

if isempty(bases) && (n == 0)
    is = [];
end

for i = 1:length(bases)
    if bases(i) == 0
        bases(i) = 1;
    end
end

if prod(bases) <= min(n)
    disp('Error in toSeveralBases : not enough digits were described for this input number');
end

if size(n,1) < size(n,2)
    n = n';
end

%if size(bases,1) < size(bases,2)
%    bases = bases';
%end

is = zeros(size(n,1),length(bases));
for j = length(bases):-1:1
    is(:,j) = mod(n, bases(j));
    n = (n - is(:,j))/(bases(j));
end

