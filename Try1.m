%% Try to plot  the residual function a_opt

%1. Defining the handle function to solve

%For the parameters, run the file "initialization.m".
%Rest of parameters not defined in the file
theta=0.4;
b=2;
h=1;
Ph=1;

a_opt = @(a, h) (b - a + w_l*(1-h) - P_h*h)^(-par.psi) - ...
                par.beta*(1+r)*( ((h*theta)^par.eta)    *(w_h + (1+r)*a)^(-par.psi) + ...
                             (1-((h*theta)^par.eta))*(w_l + (1+r)*a)^(-par.psi));
                         


%Plot the function for a grid of a from min to max
amin=-w_l/(1+r)+1.0e-2;
amin = -20;
% amax=b+w_l*(1-h)-Ph*h-1.0e-2;
amax = 40;
agrid=amin:0.1:amax;
[~,MM]=size(agrid);
yyyy=zeros(MM,1);
for ii=1:1:MM
    yyyy(ii)=a_opt(agrid(ii),1);
end
plot(agrid,yyyy)



t = fminunc(@(x) general_eq(par, gr, P_h, x(1), x(2), x(3)), [1, 3.7, 0.3], optimset('MaxIter',1))

t = fminunc(@(x) general_eq(par, gr, P_h, 1, x(1), x(2)), [3.7, 0.2], optimset('MaxIter',10))

wlgrid = 0.18:0.05:0.23;
whgrid = 1.9:0.05:2.1;
rgrid = 0.15:0.01:0.19;

mat = zeros(numel(wlgrid),numel(whgrid),numel(rgrid));
for iw_l = 1:numel(wlgrid)
    for iw_h = 1:numel(whgrid)
        for ir = 1:numel(rgrid)
            mat(iw_l, iw_h, ir) = general_eq(par, gr, P_h, wlgrid(iw_l), whgrid(iw_h), rgrid(ir));
        end
    end
end

t = fminunc(@(x) general_eq(par, gr, P_h, x(1), x(2), x(3)), [0.2950, 1.8350, 0.1750], optimset('MaxIter',1))

w_l = t(1);
w_h = t(2);
r = t(3);







% This graphs the region of the people who decide to study

% profile on

mat = zeros(gr.nb, gr.ntheta);
tic
for ib = 1:gr.nb
    for itheta = 1:gr.ntheta
        [ct, at, ht, ut] = utility(par, gr, gr.bgrid(ib), gr.thetagrid(itheta), w_l, w_h, r, P_h);
        mat(ib, itheta) = ht;
    end
end
toc
v=[1,1];

figure(2);
contourf(gr.thetagrid, gr.bgrid, mat,v)
xlabel('\theta'); ylabel('Bequests')
title(['People who study'])

% 


tic
general_eq(par, gr, P_h, w_l, w_h, r)
toc


