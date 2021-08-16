function [prob M]=Pguess(rho)

dims=size(rho);
d=dims(1);
NumberOfStates=dims(3);

cvx_begin sdp quiet
variable B(d,d,NumberOfStates) hermitian semidefinite

obj=0;
u=0;
for x=1:NumberOfStates
    u=u+B(:,:,x);
    obj=obj+trace(rho(:,:,x)*B(:,:,x))/NumberOfStates;
end
u==eye(d);

obj=(obj+obj')/2;

maximise(obj)

cvx_end

prob=cvx_optval;
M=B;

end


