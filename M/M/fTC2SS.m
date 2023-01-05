function [As,Bs,Cs,Ds,Idx] = fTC2SS(A,G,b,C,f,y)
% transforms a model from thermal circuit (TC) to state-space (SS)
%
%Inputs
% A adjancecy (TC connection ) matrix:
%   #cols = #temperature nodes; #rows = #heat flow rates
% G sqaure matrix of conductances
% b vector: 1 for branches with temperature sources, otherwise 0
% C square matrix of capacities
% f vector: 1 for nodes with heat sources, otherwise 0
% y vector: 1 for output nodes, otherwise 0
%
%Outputs
% As state matrix in state equation
% Bs input matrix in state equation
% Cs output matrix in observation equation for thC
% Ds input matrix in obs eq for thC
% Idx{1}    nodes with capacities
%       {2} branches with temp. sources
%       {3} nodes with flow sources
%       {4} nodes output temperatures

%********************
b = b'; f = f'; y = y';
rC = find(diag(C)); r0 = find(not(diag(C))); % find C which is: not-zero; zero
Idx{1} = rC'; Idx{2} = find(b); Idx{3} = find(f); Idx{4} = find(y);
if isempty(rC); error('diag(C) is zero'); end
CC = C(rC,rC); % non-zero C
K = -A'*G*A;
K11 = K(r0,r0);
K12 = K(r0,rC);
K21 = K(rC,r0);
K22 = K(rC,rC);

Kb = A'*G;
Kb1 = Kb(r0,:);
Kb2 = Kb(rC,:);

% state equation
As = inv(CC)*(-K21*inv(K11)*K12 + K22);
Bs = inv(CC)*[-K21*inv(K11)*Kb1+Kb2 -K21*inv(K11) eye(size(CC))];
sb = size(b,2);
Bs(:,[1:sb sb+r0' sb+rC']) = Bs; % rearange in order of f-sources
in = [find(b) size(A,1)+find(f)]; %effective inputs [temp_on_branches flow_in_nodes]
Bs= Bs(:,in);  %extract actual inputs (inputs <> 0)

% observation equation for y: a set of states thC
Ds = zeros(size(y(rC),2),size([b f],2));

% observation equation for y: a set of non-states th0
Cso = -inv(K11)*K12;
Dso = -inv(K11)*[Kb1 eye(size(r0,1)) zeros(size(r0,1),size(CC,1))];

Cx = zeros(size(y,2),size(As,1));
Cs = diag(y(rC));
Cx(rC,:) = Cs;
Cx(r0,:) = Cso;
Cs = Cx(find(y),:);

Dx = zeros(size(y,2),size([b f],2));
Dx(r0,:) = Dso; % feed-through if no capacity
Dx(:,[1:sb sb+r0' sb+rC']) = Dx; % rearange in order of f-sources
Ds = Dx(find(y),in);