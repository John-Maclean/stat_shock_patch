%`exact' solution for Burger's eqn on an infinite domain with an arbitrary
%initial condition.
function [uExact,uSol,nBit,dBit] = fun_Exact(tt,xx,epsilon,u0,u0int,locPeak,sup)

% epsilon = 0.01; %diffusion coefficient

%either get arguments or set them - deprecated
% if ~isempty(u0Opt)
%     u0 = u0Opt;
%     u0int=u0intOpt;
%     locPeak = locPeakOpt;
% else
%     disp('Initial condition not specified: assuming u0 = -sin(x0).')
%     u0 = @(x)-sin(x); %initial condition
%     u0int = @(x) cos(x)-1; %integral from 0 to x of u0
%     locPeak = @(t,x) fsolve(@(y)(y-x)/t-sin(y),x);
% end
% if nargin>6 %user specified a tolerance for the infinite integral
%     tol = tolOpt;
% else
%     tol = Inf;
% end



% u0intDiff = (u0int(0)-u0int(pi))/pi^2; %scaling factor for quadratic
% approximation to u0int - only needed if locPeak is not specified by
% script

partN = @(t,x,y) exp(-(x-y).^2./(4*t)-u0int(y)/(2)); %scaled part of the numerator
nBit = @(t,x,y) (x-y).*exp( (-(x-y).^2./(4*t)-u0int(y)/(2) - log(rescale(t,x)))/epsilon ); %numerator integrand
dBit = @(t,x,y) t*exp( (-(x-y).^2./(4*t)-u0int(y)/(2) - log(rescale(t,x)))/epsilon ); %denominator integrand

%note: uSol(t,x,y) is the integral from -Inf to Inf in y; but the integrand decays
%exponentially as y gets further from x. 
uSol = @(t,x) integral(@(y)nBit(t,x,y),locPeak(t,x)-sup,locPeak(t,x)+sup)/(integral(@(y)dBit(t,x,y),locPeak(t,x)-sup,locPeak(t,x)+sup));

% x = linspace(xl,xr,numX);


% t = 0:0.1:10;

nT = length(tt); nX = length(xx);
uExact=zeros(nT,nX);

uExact(1,:) = u0(xx);
for jj = 2:nT
    for kk=1:nX
        t=tt(jj); x=xx(kk);
        uExact(jj,kk) = uSol(t,x);
    end
end

% [t,x] = meshgrid(tt,xx);
% figure; set(gcf,'color','w','PaperPosition',[0 0 14 10]); hold on
% surf(t,x,uExact')


    function scaler = rescale(t,x)
        if length(t)>1
            error('the rescaling only works with a scalar value of time')
        end
        %Get a factor to rescale the exponent before raising to the power 1/epsilon.
        %The factor needs to be close to the true peak if epsilon is small.
        scaler = (partN(t,x,locPeak(t,x)));
        %if unsure of the peak location, take the max over a range instead
        %scaler = max(abs(partN(t,x,locPeak(t,x)-2:0.001:locPeak(t,x)+2)));
    end

%     function Loc = locPeak(t,x)
%         %this function finds the minimum wrt y of (x-y)^2/(2t) + uint(y) .
%         %It does so by 1) approximating uint(y) by a quadratic of the form
%         %a(|y| - pi)^2, then 2) further approximating uint(y) by
%         %a(y - pi*sign(x))^2.
%         %The minimum Loc gives the maximum value of exp(-(x-y)^2/(2t) -uint(y))
%         %and all the mass in the integral is concentrated close by.  
%         Loc = (x/t+2*pi*sign(x)*u0intDiff)/(2*u0intDiff+1/t);
%     end
end

