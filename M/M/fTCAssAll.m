function [TCa, Idx] = fTCAssAll(TCd, AssX)
% Assembled circuit TCa from dissambled circuit TCd and assembling matrix Ass
% 
% Inputs
% TCd   cell array of cell arrays of [A,G,b,C,f,y] of each TC
% AssX  assembly matrix: 4 elem/row = TC#node TC#node (2nd node collapse in 1st)
% TCa   cell array of [A,G,b,C,f,y] of assembled TC
% Idx   Idx{1}: nodes; Idx{2}: branches; 
%       1st row: local therm. circ.; 
%       2nd row: local indexes; 
%       3rd row: global indexes

% Create assembing matrix Ass from AssX
nf(1) = 0; 
nth = length(TCd{1}{5});    % no of temp nodes
tcth = ones(1,nth);         % index for nodes of thermal circ.
nq = length(TCd{1}{3});     % no of flow branches
tcq = ones(1,nq);           % index for branches of thermal circ.
gq(1) = length(TCd{1}{3});  
gth(1) = length(TCd{1}{5});
for n = 2:size(TCd,2)
    nf(n) = nf(n-1) + length(TCd{n-1}{5});  % 
    nth = nth + length(TCd{n}{5});          % total # nodes
    tcth = [tcth n*ones(1,length(TCd{n}{5}))]; % index for nodes of thermal circ.
    nq = nq + length(TCd{n}{3});            % total # branches
    tcq = [tcq n*ones(1,length(TCd{n}{3}))];         % index for branches of thermal circ.
    gq(n) = gq(n-1) + length(TCd{n}{3});    % global dissambled indexes q
    gth(n) = gth(n-1) + length(TCd{n}{5});  % global dissambled indexes th
end
gth = gq(end)+gth;                          % th counted after q

% Assembly matrix 2 columns 1st <- 2nd: global dissembled temp idx
Ass = [nf(AssX(:,1))'+AssX(:,2) nf(AssX(:,3))'+AssX(:,4)];

% Create dissembling matrix Ad (from Adth themperature + Adq flows)
Adth = eye(nth);    % a) create diagonal matrix
for n = 1:size(Ass,1)   % for each line of Ass
    % add 2 columns of Ass indicated by Ass
    Adth(:,Ass(n,1)) = Adth(:,Ass(n,1)) + Adth(:,Ass(n,2));
end
Adth(:,Ass(:,2)) = [];      % c) delete the columns corresponding to 2nd node in merging

Adq = eye(nq);

Adg = blkdiag(Adq, Adth);

Ad = [Adg(1:gq(1),:); Adg(gq(end)+1:gth(1),:)]; % rearrange Ad
for n = 2:size(TCd,2)
    Ad = [Ad;...
            Adg(gq(n-1)+1:gq(n),:); Adg(gth(n-1)+1:gth(n),:)];
end

% Assemble circuits
Kd = [inv(TCd{1}{2}) TCd{1}{1}; -TCd{1}{1}' TCd{1}{4}]; % Kd = [inv(G1) A1; -A1' C1]
ubf = [TCd{1}{3}; TCd{1}{5}];                           % ubf = [b1; f1]
uby = [TCd{1}{3}; TCd{1}{6}];                           % uby = [b1; y1]
for n = 2:size(TCd,2)
    Kd = blkdiag(Kd,[inv(TCd{n}{2}) TCd{n}{1}; -TCd{n}{1}' TCd{n}{4}]);
    ubf = [ubf; TCd{n}{3}; TCd{n}{5}];
    uby = [uby; TCd{n}{3}; TCd{n}{6}];
end

Ka = Ad'*Kd*Ad;
Ga = inv(Ka(1:nq,1:nq));
Aa = Ka(1:nq,nq+1:end);
Ca = Ka(nq+1:end, nq+1:end);

u = Ad'*ubf; % inputs
ba = u(1:nq);
fa = u(nq+1:end);
fa = fa~=0;     % % 1 <- values >=1
u = Ad'*uby;    %outputs
ya = u(nq+1:end);   
ya = ya~=0;     % 1 <- values >=1

TCa = {Aa,Ga,ba,Ca,fa,ya};  % cell array of assembled TC

% Find indexes local - global
[i,j] = find(Adth);   % global indexes: i -before assembling; j -after assembling
thag = sortrows([i j],1);
qag = 1:nq;

thl = 1:length(TCd{1}{5});  % local index of temperatures
ql = 1:length(TCd{1}{3});   % local indexes of flow    
for n = 2:size(TCd,2)
    thl = [thl 1:length(TCd{n}{5})];
    ql = [ql 1:length(TCd{n}{3})];
end
% temp 1st row: index of local therm. circ; 
% 2nd rox: local indexes; 3nd row global idx of assembled
thidx = [tcth; thl; thag(:,2)'];  
% flow 1st row: index of local therm. circ; 
% 2nd row: local indexes; 3rd row global idx of assembled
qidx = [tcq; ql; qag];           
Idx = {thidx,qidx};
