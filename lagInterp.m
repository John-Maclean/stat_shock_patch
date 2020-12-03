%get a fast function to evaluate a Lagrange poly at a single point 
%Inputs: 
%xx, fixed x values at which height is observed, should be a column vector
%x is the point at which the poly is evaluated
%Output: function L(y)
%given function heights y **in a column vector**, L(y) interpolates between these heights to give 
%an estimate at x.
function L=fastLagrange(xx,x)
%%%test values
% xx=(1:6)';
% yy=sin(xx);
% xxx= 0:.1:7;
% 
% for jj=1:length(xxx)
%     x=xxx(jj);

nn=length(xx);
xm=zeros(nn,nn-1);

for n=1:length(xx)
xm(n,:) = xx([1:n-1 n+1:end]);
end

lj= prod( (x-xm)./(xx-xm),2);

L = @(y) sum(y.*lj);
% L(jj) = sum(yy.*lj);
% end %end test loop

% figure %plot test results
% hold on
% plot(xx,yy,'bo')
% plot(xxx,L,'k.')
