

%-------------------------------------------------------------------------%
%                   Housekeeping                                          %
%-------------------------------------------------------------------------%
clear
clc

%-------------------------------------------------------------------------%
%                   Parameters                                            %
%-------------------------------------------------------------------------%

% global par

%--------------------------%
%   Fixed parameters       %
%--------------------------%

% Utility function
par.beta    = 0.971;         % Discount factor
par.psi     = 2;            % Risk aversion

% Government policies
par.tau_l   = 0;        % Tax to unskilled labor
par.tau_h   = 0;        % Tax to skilled labor
par.s       = 0;        % Subsidy to education

% Production technology
par.A       = 1;        % Technological parameter
par.alpha_h = 0.5;      % Weight of skilled labor
par.z       = 1.5;      % Productivity of skilled labor
par.phi     = 2;        % Elasticity of substitution: skilled vs unskilled
par.gamma   = 0.66;     % Share of labor 

% Other parameters
par.eta     = 0.5;      % Concavity of probability of completing college


%----------------------------------------%
%   General equilibrium parameters       %
%----------------------------------------%

% Prices (partial equilibrium)
w_l     = 0.2955;        % Wage rate - unskilled labor
w_h     = 1.8346;        % Wage rate - skilled labor
r       = 0.1754;     % Interest rate
P_h     = 1;      % Price of college

%--------------------------%
%          Grids           %
%--------------------------%

% Bequests grid
gr.nb       = 35;
gr.bgrid    = linspace(0, 5, gr.nb);           % Grid for health

% Bonds
gr.na      = 25;
gr.Abar    = 0;                                 % Exogenous borrowing constraint
gr.Amax    = 20;
% gr.agrid   = linspace(-gr.Abar, gr.Amax, gr.na);      % Grid for savings

% Ability: theta
gr.ntheta       = 35;
gr.thetagrid    = linspace(0, 1, gr.ntheta);          % Grid for savings






