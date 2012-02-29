function [all_names, AAs ,  Cs , ineq ,cost ,in_vars,all_state_vars,ex_vars] = ExtractModel(n_time_steps,dt,integ_method)



if (nargin < 1)
    n_time_steps = 3;
end

if (nargin < 2)
    discrete_sys = true;
else
    discrete_sys = false;
end

if (nargin < 3)
    integ_method = 'Euler'
end
% First, get structure for one time step
[all_names_single,  AAsingle,  Csingle , state_vars_tmp , ineq_vars_single , cost_vars_single,in_vars_single, ex_vars_single, core_vars,core_functions ] = ExtractOnetimeStep;

%reverse the direction of state vars data : 
state_vars = sparse(1,length(state_vars_tmp),0);
for i= find(state_vars_tmp)
    state_vars(state_vars_tmp(i)) = i;
end

% duplicate for many time steps
idx_state_vars = find(state_vars);
n_state_vars = length(idx_state_vars);
N = size(AAsingle{1},2)*n_time_steps ; %-n_state_vars*(n_time_steps-1) ;
n = size(AAsingle{1},2);
all_names = {};
AAs = {};
Cs = {};

ineq_vars = sparse(1,N,0);

for t=1:n_time_steps
    for i=1:length(AAsingle)
        AAs{end+1} = sparse(size(AAsingle{i},1),N,0); 
        AAs{end}(:,(t-1)*n+[1:n]) = AAsingle{i};
        Cs{end+1} = Csingle{i};
    end
    
    for i=1:length(all_names_single)
        all_names{end+1} = [ all_names_single{i} '.' 't' num2str(t) ];
    end
    ineq_vars((t-1)*n+(1:n)) =  ineq_vars_single;
    cost_vars((t-1)*n+(1:n)) =  cost_vars_single;
    in_vars  ((t-1)*n+(1:n)) =  in_vars_single;
    ex_vars  ((t-1)*n+(1:n)) =  ex_vars_single;
    all_state_vars  ((t-1)*n+(1:n)) =  state_vars;
end


if (discrete_sys)
    % introduce equality constraint at state variables between time steps
    % Calc degree of each state variable - propagate the history backwards for
    % pipe of state blocks
    degree = sparse(1,length(state_vars),0);
    for k=1:n_state_vars
        current = idx_state_vars(k);
        while (state_vars(current)~=0)
            current = state_vars(current);
            degree(idx_state_vars(k)) = degree(idx_state_vars(k))  +1;
        end
    end

    toremove_list = [];

    for d=1:full(max(degree))
        for t=d+1:n_time_steps
            for k=1:n_state_vars
                if degree(idx_state_vars(k)) ~= d
                    continue;
                end
                source = state_vars(idx_state_vars(k));
                for i=2:d
                    source = state_vars(source);
                end
                AAs = MoveEqualVar(AAs,source+(t-1-d)*n,idx_state_vars(k)+(t-1)*n); % just copy data, do not remove the column yet
                toremove_list = [toremove_list idx_state_vars(k)+(t-1)*n];
            end
        end
    end
else % continous system
    % Discretize
    switch (integ_method)
        case 'Euler'
            ButcherTableau.A = [0];
            ButcherTableau.c = 0;
            ButcherTableau.b = 1;
        case 'Trapez'
            ButcherTableau.A = [0 0 ; .5 .5];
            ButcherTableau.c = [0 1]';
            ButcherTableau.b = [0.5 0.5];
        case 'RK4'
            ButcherTableau.A = [ 0 0 0 0 ; 
                                .5 0 0 0 ;
                                 0 .5 0 0 ;
                                 0 0 1 0];
            ButcherTableau.c = [0 0.5 0.5 1]';
            ButcherTableau.b = [1/6 1/3 1/3 1/6 ];
    end

    toremove_list = [];
    for i=2:n_time_steps
        [AAs,Cs,toremove_list_tmp,all_names] = CreateOneIntegrationStep(AAs,Cs, core_vars,core_functions,ButcherTableau,...
            (i-2)*n,(i-1)*n,n,state_vars,dt,in_vars_single,ex_vars_single,all_names,all_names_single);
        toremove_list = [toremove_list toremove_list_tmp];
        Nprev = N;
        N = length(all_names);
        if (N > Nprev)
            ineq_vars = [ineq_vars sparse(1, N - Nprev)];
            cost_vars = [cost_vars sparse(1, N - Nprev)];
            in_vars   = [in_vars sparse(1, N - Nprev)];
            ex_vars   = [ex_vars sparse(1, N - Nprev)];
            all_state_vars = [all_state_vars sparse(1, N - Nprev)];
        end
        
        all_state_vars(n+1:N) = 0; % only initial conditions remains
        
    end

end

% filter out input variables that are fixed for more than one time step
idx = find(in_vars_single(1:n) > 1 );
for i=1:length(idx)
    for tt = 1:ceil(n_time_steps/in_vars_single(idx(i)));
        t = 1 + (tt-1)*in_vars_single(idx(i));
        next_vars = idx(i)+[t:(t+in_vars_single(idx(i))-2)]*n;
        next_vars = next_vars(next_vars <= length(all_names)); % cut all variables after the last variable
        if isempty(next_vars)
            continue;
        end
        
        AAs = MoveEqualVar(AAs,idx(i)+(t-1)*n,next_vars ); % just copy data, do not remove the column yet
        toremove_list = [toremove_list next_vars];
    end
end

Ns = ones(1,N);
Ns(toremove_list) = 0;
idx = find(Ns);
for i=1:length(AAs)
    AAs{i} = AAs{i}(:,idx);
end

all_names = {all_names{idx}};
ineq_vars = ineq_vars(idx);
cost_vars = cost_vars(idx);
in_vars   = in_vars(idx);
ex_vars   = ex_vars(idx);
all_state_vars = all_state_vars(idx);



[all_names, AAs ,  Cs, ineq , cost,in_vars,all_state_vars,ex_vars]=MoveEqToCostAndNeq(all_names, AAs ,  Cs , ineq_vars ,cost_vars,in_vars,all_state_vars,ex_vars);

return

function [AAs, Cs, idx_to_stay , cost_vars, all_names] = EliminateIdentityConstraints(AAs,Cs,cost_vars,all_names)

toremove_list = [];

for i=1:length(AAs)
    for j=1:size(Cs{i},1) % for all output variables
        C = Cs{i}(j,:);
        A = AAs{i}(C~=0,:);
        if (sum(C)==0 && (max(C) == 1) && (min(C) == -1) ...
                && (size(A,1) ==2) && (max(A(:))==1 )&&(sum(A(:))==2 ) ...
                && (length(find(A(1,:))~=0)==1) && (length(find(A(2,:))~=0)==1) )
            to_remove = find(A(2,:));
            origin=find(A(1,:));
            if ~isempty(find(origin==toremove_list,1)) % handle the case that the origin variable is already removed
                continue;
            end
            AAs = MoveEqualVar(AAs,origin,to_remove); % just copy data, do not remove the column yet
%             AAs{i}(C~=0,:) = 0; % reset the identity polyblock
            toremove_list = [toremove_list to_remove ];
            cost_vars(origin) = cost_vars(to_remove);
            
        end
    end
end

to_stay = ones(1,length(cost_vars));
to_stay(toremove_list) = 0;
idx_to_stay = find(to_stay);

[AAs,Cs,all_names] = RemoveVariable(AAs,Cs,all_names,toremove_list);

remove_functons  = [];
for i=1:length(AAs)
    remove_lines = [];
    for j=1:size(Cs{i},1) % for all output variables
        C = Cs{i}(j,:);
        A = AAs{i}(C~=0,:);

        if (sum(C)==0 && (max(C) == 1) && (min(C) == -1) ...
                && (size(A,1) ==2) && (max(A(:))==1 )&&(sum(A(:))==2 ) ...
                && (length(find(A(1,:))~=0)==1) && (length(find(A(2,:))~=0)==1) )
            if  find(A(1,:)) == find(A(2,:))
                remove_lines = [remove_lines j];
            end
        end
    end
    idx = ones(1,size(Cs{i},1));
    idx(remove_lines)=0;
    if length(remove_lines) == length(idx)
        remove_functons = [remove_functons i];
    else
        Cs{i} = Cs{i}(idx==1,:);
    end
end
idx = ones(1,length(AAs));
idx(remove_functons)=0;
AAs = { AAs{idx==1} };
Cs  = { Cs{idx==1} };



function [all_names, AAs ,  Cs, ineq , merged_cost,in_vars,state_vars,ex_vars]=  MoveEqToCostAndNeq(all_names, AAs ,  Cs , ineq_vars ,cost_vars,in_vars,state_vars,ex_vars)

[cost, to_remove_var_cost ,all_names, AAs ,  Cs] =  ExtractVars(all_names, AAs ,  Cs ,cost_vars);
[ineq, to_remove_var_ineq ,all_names, AAs ,  Cs] =  ExtractVars(all_names, AAs ,  Cs ,-ineq_vars); % minus is to handle definition of inequality

[AAs,Cs,all_names_new] = RemoveVariable(AAs,Cs,all_names, [to_remove_var_ineq to_remove_var_cost]);
[cost.AAs,cost.Cs] = RemoveVariable(cost.AAs,cost.Cs,all_names, [to_remove_var_ineq to_remove_var_cost]);
[ineq.AAs,ineq.Cs] = RemoveVariable(ineq.AAs,ineq.Cs,all_names, [to_remove_var_ineq to_remove_var_cost]);

stay_places = ones(1,length(all_names));
all_names = all_names_new;
stay_places( [to_remove_var_ineq to_remove_var_cost]) = 0;
in_vars = in_vars(stay_places==1);
ex_vars = ex_vars(stay_places==1);
state_vars = state_vars(stay_places==1);
% merge all cost functions
merged_cost.A = [];%sparse(0,length(all_names),0);
merged_cost.C = [];%sparse(0,1,0);
for i=1:length(cost.AAs)
    merged_cost.A((end+1):(end+length(cost.Cs{i})),:) = cost.AAs{i};
    merged_cost.C((end+1):(end+length(cost.Cs{i}))) = cost.Cs{i};
end


function [new_func, to_remove_var ,all_names, AAs ,  Cs] =  ExtractVars(all_names, AAs ,  Cs ,vars)

to_remove_var = [];
new_func.AAs = {};
new_func.Cs = {};

idx = find(vars);
for i=1:length(idx)
    usage =0;
    last_used = 0;
    for j=1:length(AAs)
        if sum(abs(AAs{j}(:,idx(i))))>0
          usage  = usage  + 1;  
          last_used = j;
        end
    end
    
    unfolded = 0;
    
    if (usage == 1) % special calculation, only for cost purpose
        % check if can be unfolded - no high powers or exponents
        powers = AAs{last_used}(:,idx(i));
        if max(powers)==1 && min(powers)>-1 && (sum(powers) == 1)
            term = find(powers);
            if (length(term) == 1) && (sum(AAs{last_used}(term,:)) == 1) % no other variables in the term
                new_func.AAs{i} = AAs{last_used}([1:term-1 term+1:end],:);
                % find where this term is used in C
                c_line = find(Cs{last_used}(:,term));
                if (length(c_line) ~= 1)
                    error('Something is wrong,length(c_line) ~= 1');
                end
                new_func.Cs{i} = - vars(idx(i))*Cs{last_used}(c_line,[1:term-1 term+1:end])/Cs{last_used}(c_line,term);
                to_remove_var = [to_remove_var idx(i)];
                % if this is the only row in C, remove function,
                if (size(Cs{last_used},1)==1)
                    AAs = {AAs{1:last_used-1} AAs{last_used+1:end}};
                    Cs =  {Cs{1:last_used-1} Cs{last_used+1:end}};
                else % remove only the row from C
                    Cs{last_used} =  Cs{last_used}( [1:c_line-1 c_line+1:end],:);
                end
                    
                unfolded = 1;
            end
        end
    end
    
    if (unfolded == 0) %  just use the existing variable
        new_func.AAs{i} = sparse(1,length(all_names),0);
        new_func.AAs{i}(1,idx(i)) = 1;
        new_func.Cs{i} = vars(idx(i));
    end
end




%==========================================================================
%==========================================================================
%==========================================================================

function [all_names, AAs ,  Cs  , state_vars , ineq_vars ,cost_vars ,in_vars, ex_vars, core_vars,core_functions]...
    = ExtractOnetimeStep


blks = find_system(gcs, 'Tag', 'PolyBlock');
mem_blks = find_system(gcs, 'Tag', 'OptState');
in_blks = find_system(gcs, 'Tag', 'OptInput');
ex_blks = find_system(gcs, 'Tag', 'OptExternal');
demuxes = find_system(gcs, 'BlockType', 'Demux');
muxes = find_system(gcs, 'BlockType', 'Mux');


cost_blks = find_system(gcs, 'Tag', 'OptCost');
ineq_blks = find_system(gcs, 'Tag', 'InequalBlock');

% group all blocks. muxes and demuxes  will be treated as polyblocks later.
all_blks = {blks{:}, demuxes{:} , muxes{:}, mem_blks{:}, in_blks{:}, ex_blks{:}  };

% this is required for port dimension 
eval([gcs '([],[],[],''compile'');'])

% Get all common block properties
for i= 1:length(all_blks)
    
    connect{i} = get_param(all_blks{i},'PortConnectivity');
    dim{i}     = get_param(all_blks{i}, 'CompiledPortDimensions');
    if (i > length(blks)+length(mem_blks)+length(demuxes)+length(muxes))  % store only output port for other than polyblocks and states
        connect{i} = connect{i}(2);
    end
    names{i} =  GetBlockName(all_blks{i});
    handles(i) = get_param(all_blks{i},'handle');
   
end

%required to unlock the diagram
eval([gcs '([],[],[],''term'');']);

% store all functions from polyblocks
for i= 1:length(blks)
    userdata =  get_param(all_blks{i},'UserData');
    
    As{i}  = userdata{1};
    Cs{i}  = userdata{2};
end

% represent demux as a identity function.
for i= length(blks)+1:length(blks)+length(demuxes)
    
    As{i}  = eye(dim{i}.Inport(2)*2);
    Cs{i}  = [eye(dim{i}.Inport(2)) , -eye(dim{i}.Inport(2))];
end

% represent mux as a identity function.
for i= length(blks)+1+length(demuxes):length(blks)+length(demuxes)+length(muxes)
    
    As{i}  = eye(dim{i}.Outport(2)*2);
    Cs{i}  = [eye(dim{i}.Outport(2)) , -eye(dim{i}.Outport(2))];
end

% from here muxes and demuxes are treated as polyblocks.
blks = {blks{:} ,demuxes{:} , muxes{:}};

mem_handles = [];
cost_handles = [];
ineq_handles = [];

for i= 1:length(mem_blks)
    mem_names{i} =  GetBlockName(mem_blks{i});
    mem_handles(i) = get_param(mem_blks{i},'handle');
end

for i= 1:length(ineq_blks)
    ineq_handles(i) = get_param(ineq_blks{i},'handle');
    if (strcmp(get_param(ineq_blks{i}, 'CompareSign'),'> 0'))
        ineq_sign(i)  =  1;
    else
        ineq_sign(i)  =  -1;
    end
end

for i= 1:length(cost_blks)
    cost_handles(i) = get_param(cost_blks{i},'handle');
end





% first phase : combine all matrices

N = 0;

% handle all polyblocks
for i=1:length(blks)
    IIs{i} = N+1:(N + size(As{i},2));
    nInVars = size(As{i},2)- size(Cs{i},1);
    nOutVars = size(Cs{i},1);
    for k=1:nInVars;
        all_names{N+k} = sprintf('%s.In%d',names{i},k);
    end
   
    for k=nInVars+1:nInVars+nOutVars;
        all_names{N+k} = sprintf('%s.Out%d',names{i},k-nInVars);
    end

    k=0;
    for in = 1:length(dim{i}.Inport)/2
        for j=1:dim{i}.Inport(in*2)
            k=k+1;
            VarConnect(N+k) = connect{i}(in);
            VarNum(N+k) = j;
            InPort(N+k) = in;
            BlockHandle(N+k) = handles(i);
            BlockType(N+k)=1;
        end
    end
    for out = 1:length(dim{i}.Outport)/2
        for j=1:dim{i}.Outport(out*2)
            k=k+1;
            VarConnect(N+k) = connect{i}(length(dim{i}.Inport)/2+out);
            VarNum(N+k) = j;
            InPort(N+k) = -out;
            BlockHandle(N+k) = handles(i);
            BlockType(N+k)=1;
        end
    end
    %     all_names{N+size(As{i},2)} = sprintf('%s.Out%d',names{i},1);
    N = N + size(As{i},2);
end

% handle all state blocks
for i=(length(blks)+1):(length(blks)+length(mem_blks))
    for k=1:dim{i}.Inport(2)
        all_names{N+k} = sprintf('%s.In%d',names{i},k);
    end
    for k=1:dim{i}.Outport(2)
        all_names{N+dim{i}.Inport(2)+k} = sprintf('%s.Out%d',names{i},k);
    end
    k=0;
    for in = 1:length(dim{i}.Inport)/2
        for j=1:dim{i}.Inport(in*2)
            k=k+1;
            VarConnect(N+k) = connect{i}(in);
            VarNum(N+k) = j;
            InPort(N+k) = in;
            BlockHandle(N+k) = handles(i);
            BlockType(N+k)=2;
        end
    end
    for out = 1:length(dim{i}.Outport)/2
        for j=1:dim{i}.Outport(out*2)
            k=k+1;
            VarConnect(N+k) = connect{i}(length(dim{i}.Inport)/2+out);
            VarNum(N+k) = j;
            InPort(N+k) = -out;
            BlockHandle(N+k) = handles(i);
            BlockType(N+k)=2;
        end
    end
    IIs{end+1} = N + [ 1:dim{i}.Inport(2) , (dim{i}.Inport(2)+[1:dim{i}.Outport(2)])];
    N = N + dim{i}.Outport(2) + dim{i}.Inport(2);
end

% handle input and external blocks
for i=(length(blks)+length(mem_blks)+1):length(all_blks)
    IIs{end+1} = N+[ 1:dim{i}.Outport(2)];
    for k=1:dim{i}.Outport(2)
        all_names{N+k} = sprintf('%s.Out%d',names{i},k);
    end
    k=0;
    for out = 1:length(dim{i}.Outport)/2
        for j=1:dim{i}.Outport(out*2)
            k=k+1;
            VarConnect(N+k) = connect{i}(out);
            VarNum(N+k) = j;
            InPort(N+k) = -out;
            BlockHandle(N+k) = handles(i);
            BlockType(N+k)=3;
        end
    end
    N = N + dim{i}.Outport(2);
end


% exapnd functions to functions of all variables
for i=1:length(blks)
    AAs{i} = sparse(size(As{i},1),N,0);
    AAs{i}(:,IIs{i}) = As{i};
end


% Store step ratio of input vars
in_vars = sparse(1,N,0);
for i= 1:length(in_blks)
    in_vars(IIs{length(blks) + length(mem_blks)+i}) = evalin('base',get_param(in_blks{i},'step_ratio'));
end

% Mark external vars
ex_vars = sparse(1,N,0);
for i= 1:length(ex_blks)
    ex_vars(IIs{length(blks) + length(mem_blks)+length(in_blks) + i}) = 1;
end



ineq_vars = sparse(1,N,0);
cost_vars = sparse(1,N,0);

% second phase : eliminate equal variables
toremove_list = [];

for i=1:length(VarConnect) % all variables
    if ~isempty(VarConnect(i).DstBlock) % all outgoing variables
        for k = 1:length(VarConnect(i).DstBlock) % all destination blocks
           idx      = find(BlockHandle == VarConnect(i).DstBlock(k) ...
                               & InPort == VarConnect(i).DstPort(k)+1  ...
                               & VarNum == VarNum(i) & BlockType ~= 3); % destination variable
           ineq_idx = find(ineq_handles == VarConnect(i).DstBlock(k)); % destination ineq block
           cost_idx = find(cost_handles == VarConnect(i).DstBlock(k)); % destination cost block

            if ~isempty(idx) % destination is a polyblock or state
                AAs = MoveEqualVar(AAs,i,idx); % just copy data, do not remove the column yet
                toremove_list = [toremove_list idx];
            end
            if ~isempty(ineq_idx)
                ineq_vars(i) =   ineq_idx;
            elseif ~isempty(cost_idx)
                cost_vars(i) =   cost_idx;
            end
        end
    end
end



Ns = ones(1,N);
Ns(toremove_list) = 0;
idx = find(Ns);
for i=1:length(AAs)
    AAs{i} = AAs{i}(:,idx);
end

all_names = {all_names{idx}};
BlockHandle = BlockHandle(idx);
InPort = InPort(idx);
BlockType = BlockType(idx);
VarNum = VarNum(idx);
VarConnect = VarConnect(idx);

ineq_vars = ineq_vars(idx);
cost_vars = cost_vars(idx);
in_vars   = in_vars(idx);
ex_vars   = ex_vars(idx);

prev_idx = [];
while length(idx) ~= length(prev_idx)
    prev_idx = idx;
    [AAs, Cs, idx , cost_vars,all_names] = EliminateIdentityConstraints(AAs,Cs,cost_vars,all_names);
    
    BlockHandle = BlockHandle(idx);
    InPort = InPort(idx);
    BlockType = BlockType(idx);
    VarNum = VarNum(idx);
    VarConnect = VarConnect(idx);
    ineq_vars = ineq_vars(idx);
    in_vars   = in_vars(idx);
    ex_vars   = ex_vars(idx);
    cost_vars = cost_vars(idx);
end

state_vars = sparse(1,length(all_names),0);

% another loop to handle all state vars.
for i=1:length(VarConnect)
    if ~isempty(VarConnect(i).DstBlock)
        for k = 1:length(VarConnect(i).DstBlock)
           mem_idx  = find(BlockHandle == VarConnect(i).DstBlock(k) ...
                               & -InPort == VarConnect(i).DstPort(k)+1  ...
                               & VarNum == VarNum(i) & BlockType == 2);
            if ~isempty(mem_idx)
                state_vars(i) =  mem_idx;
            end
        end
    end
end



idx = find(ineq_vars);
for i=1:length(idx)
    ineq_vars(idx(i)) = ineq_sign(ineq_vars(idx(i)));
end

[core_vars,core_functions]  = FindCoreModel(AAs,Cs,state_vars,ex_vars,in_vars,ineq_vars,cost_vars);


function [AAs,Cs,names] = RemoveVariable(AAs,Cs,names,vars_to_remove)


stay_places = ones(1,length(names));
stay_places(vars_to_remove) = 0;
stay_idx= find(stay_places);

for i=1:length(AAs)
    AAs{i} = AAs{i}(:,stay_idx);
%     Cs{i} = Cs{i}; 
end

names = {names{stay_idx} };



function[core_vars,core_functions]  = FindCoreModel(AAs,Cs,state_vars,ex_vars,in_vars,ineq_vars,cost_vars)

core_vars = -1*ones(1,size(AAs{1},2));
% core_vars(ineq_vars ~=0 ) = -1;

core_vars(state_vars ~= 0) = 1e6;  % never remove state vars and derivatives
core_vars(state_vars(state_vars ~= 0) ) = 1e6;  % never remove state vars and derivatives
core_vars(ex_vars~=0 ) = core_vars(ex_vars~=0 ) + 1;
core_vars(in_vars~=0 ) = core_vars(in_vars~=0 ) + 1;


% count all variables appearance in functions
for i=1:length(AAs)
    for j=1:size(Cs{i},1)
        core_vars = core_vars + ( sum(abs(AAs{i}(Cs{i}(j,:)~=0,:)),1) > 0) ;
    end
end

removed = true;
while removed
    removed = false;
    % check if there are functions that uses 0-grade variables
    for i=1:length(AAs)
        for j=1:size(Cs{i},1)
            used_vars =  ( sum(abs(AAs{i}(Cs{i}(j,:)~=0,:)),1) > 0) ;
            idx =  find(core_vars(used_vars~=0)<=0,1); % zero grade vars
            if ~isempty(idx)
                % decrease use counter of those variables
                core_vars(used_vars~=0) = core_vars(used_vars~=0) -1;
                % remove this function
                Cs{i}(j,:) = 0;
                removed = true;                
            end
        end
    end
end

% filter all zero functions - copy just the non zeros

k = 1;
for i=1:length(AAs)
    m=1;
    for j=1:size(Cs{i},1)
        if (~isempty(find(Cs{i}(j,:),1)))
            core_functions.Cs{k}(m,:) = Cs{i}(j,:);
            m = m+1;
        end
    end
    if m > 1
        core_functions.AAs{k} = AAs{i};
        k = k + 1;
    end
end

core_vars = core_vars > 0;


function  [AAs,Cs,remove_var,all_names] = CreateOneIntegrationStep(AAs,Cs, core_vars,core_functions,ButcherTableau,...
        prev_step,new_step,n,state_vars,dt,in_vars,ex_vars,all_names,all_names_single)
    
states = find(state_vars );
derivs = state_vars(state_vars~=0);

remove_var = [];

k_start = [];

if  (   (size( ButcherTableau.A,1) == 2) ...
        && (sum(abs(ButcherTableau.A(1,:)))==0) ...
        && sum(ButcherTableau.A(2,:) ~= ButcherTableau.b)==0 ...
        && (ButcherTableau.c(1) == 0)) % trapezoidal

    k_start = [prev_step new_step ];
% elseif isempty(ButcherTableau.A) % euler
%     k_start = prev_step;
else

    
    for i=1:size(ButcherTableau.A,1) % for each k_i
        if sum(abs(ButcherTableau.A(i,:)))==0 % use time step data
            if (ButcherTableau.c(i) < 1)
                 k_start(i) = prev_step;
            else
                k_start(i) = new_step;
            end
        else
            k_start(i)  = size(AAs{1},2); % remember the previous size
            [AAs,Cs] = AppendVarsFuncs(AAs,Cs, core_functions); % Add core_functions and expand all matrices
            remove_var = [remove_var (k_start(i)+find(core_vars==0)) ]; % remember what to remove later
        end
    end

    % create equality constraint for k_i input
    for i=1:length(k_start)
        if (k_start(i) == prev_step) || (k_start(i) == new_step)
            continue;
        end
        % y_k = y_n + h*sum(a_ij * k_j)
        A = sparse(length(states),size(AAs{1},2));
        C = sparse(length(states),length(states));
        for s=1:length(states) % y_n
            A(s,states(s)+prev_step)=1;
            C(s,s) = 1;
        end

        for s=1:length(states) % -y_k
            A(length(states)+s,states(s)+k_start(i))=1;
            C(s,length(states)+s) = -1;
        end

        for j=1:length(k_start) % h*sum(a_ij * k_j)
            for s=1:length(derivs) % k_j
                A(end+1,derivs(s)+k_start(j))=1;
                C(s,size(A,1)) = ButcherTableau.A(i,j)*dt;
            end
        end
        AAs = {AAs{:} ,A};
        Cs = {Cs{:} ,C};
        
        % expand names
        for j=1:length(all_names_single)
            all_names{end+1} = [ all_names_single{j} '.' 't' num2str(new_step/n+1) 'rk' num2str(i) ];
        end

        % connect all input and external vars to the original time step
        % values
        for v=[find(in_vars),find(ex_vars)]
            if (ButcherTableau.c(i) < 1) % connect to the prev time step if c < 1
                source = prev_step + v;
            else
                source = new_step + v;
            end
            
            AAs = MoveEqualVar(AAs,source,k_start(i)+v); % just copy data, do not remove the column yet
            remove_var = [remove_var  (k_start(i)+v)];
        end


    end
end

% create equality constraint for y_(n+1)
% y_(n+1) = y_n + h*sum(b_i * k_i)
A = sparse(length(states),size(AAs{1},2));
C = sparse(length(states),length(states));
for s=1:length(states) % y_n
    A(s,states(s)+prev_step)=1;
    C(s,s) = 1;
end

for s=1:length(states) % -y_(n+1)
    A(length(states)+s,states(s)+new_step)=1;
    C(s,length(states)+s) = -1;
end

for j=1:length(k_start) % h*sum(b_i * k_i)
    for s=1:length(derivs) % k_j
        A(end+1,derivs(s)+k_start(j))=1;
        C(s,end+1) = ButcherTableau.b(j)*dt;
    end
end

AAs = {AAs{:} ,A};
Cs = {Cs{:} ,C};

function [AAs,Cs] = AppendVarsFuncs(AAs,Cs, core_functions)
% Add core_functions and expand all matrices
prev = size(AAs{1},2);
N = size(AAs{1},2)+size(core_functions.AAs{1},2);
for i=1:length(AAs)
       [I,J,s] = find(AAs{i});
       [m,n] = size(AAs{i});
       AAs{i} = sparse(I,J,s,m,N);
end

for i=1:length(core_functions.AAs)
       [I,J,s] = find(core_functions.AAs{i});
       [m,n] = size(core_functions.AAs{i});
       AAs{end+1} = sparse(I,J+prev,s,m,N);
       Cs{end+1} = core_functions.Cs{i};
end