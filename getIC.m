function [u0,u0int,locPeak,nPatch,sProp] = getIC(str,epsilon)

if strcmp(str, 'sinusoidal')
    u0 = @(x) -sin(x); %IC
    u0int = @(x) cos(x)-1; %integral of IC - used for quadrature
    locPeak = @(t,x) tay_sinIC(t,x,tay_sinIC(t,x,tay_sinIC(t,x,x)));
    %locPeak improves accuracy of quadrature by locating the mass of the integrand
    nPatch = 34; %minimum total number of patches (of both types) across domain
    sProp = 15; %width of shock patch as multiple of normal patches

elseif strcmp(str, 'tanh')
    delta=epsilon; %width of initial shock
    u0 = @(x) (x/pi -tanh(x/delta))/tanh(pi/delta);
    u0int = @(x) (delLogcosh(pi/delta,delta) - pi/2  + x.^2/2/pi - delLogcosh(x/delta,delta)  )/tanh(pi/delta); %integral of the initial condition
    u0intDiff = (u0int(0)-u0int(pi))/pi^2;
    locPeak = @(t,x) (x/t+2*pi*sign(x)*u0intDiff)/(2*u0intDiff+1/t);
    nPatch = 5; %minimum total number of patches (of both types) across domain
    sProp = 5; %width of shock patch as multiple of normal patches

    
else
    error('Initial conditions not recognised.')
end

%% minor functions
    function out = delLogcosh(x,del) %evaluate delta*log(cosh(x)) at small and large x
        % out = del*log(cosh(x));
        %
        % % out = log( (exp(x/del) + exp(-x/del)).^del ) - del*log(2); %no better than original
        %
        out = zeros(size(x));
        ind = (abs(x)>20);
        
        out(ind) = del*(abs(x(ind))-log(2));
        out(~ind) = del*log(cosh(x(~ind)));
    end

    function ymin = tay_sinIC(t,x,a)
        %This function is intended for use with the initial condition -sin(x).
        %In order to get finite results for two integrals, we need to rescale around the
        %peak of an exponential of (x-y)^2/(2t) + cos(y).
        %This function minimises the taylor expansion (around a) of
        %(x-y)^2/(2t) + cos(y) (below u = y-a)
        
        c=cos(a); s=sin(a);
        umin = fminbnd(@(u) u*(a/t-x/t-s) + u^2*(1/t-c)/2 + u^3*s/6 + u^4*c/24,-2,2);
        ymin = umin+a;
    end
end