function n = fromSeveralBases(is, bases)
% n = fromSeveralBases(is, bases)
%
% Gives a number in basis 10 from a list of digit, each of which having a
% different basis.
% 
% Exemple : 21 = 1*(3*4) + 2*(4) + 1 = fromSeveralBases([1 2 1], [2 3 4])
%
% Also accepts several numbers to be changed simultaneously:
% fromSeveralBases([1 2 1; 1 1 1], [2 3 4]) returns [21; 17]

% Written by Jean-Daniel Bancal <= 2011, last modified 26 May 2014.

bases(bases==0) = 1;

if prod(double(bases <= max(is)))
    disp('Error in fromSeveralBases : not enough digits were described for this input number');
    return;
end

% % One way to do it
%bases2 = cumprod([1 bases]);
%n = sum(bases2(1:length(bases)).*fliplr(is));
%return;

if (size(is,2) == 1) && (size(bases,2) ~= 1)
    is = is';
end

% Uglier but faster (!)
coeff = 1;
n = zeros(size(is,1),1);
for j = length(bases):-1:1
    n = n + coeff*is(:,j);
    coeff = coeff*bases(j);
end

