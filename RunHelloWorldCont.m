clear
clc
mode = 'fmincon' ;


open('HelloWorldCont');

%% Convert model to optimization problem
[all_names, AAs ,  Cs , ineq_vars ,cost_vars,in_vars,state_vars, ex_vars] = ExtractModel(10,1,'RK4');
sim('HelloWorldCont');
%% Print problem - just for display
code = PrintProblem(all_names, AAs ,  Cs ,ineq_vars,cost_vars);
for i=1:length(code)
    disp(code{i});
end

%% create matlab optimization format
if (strcmp(mode,'fmincon'))
    disp('Creating fmincon callbacks...');
        
        [cost_name costGrad_name eqconstr_name eqconstrGrad_name neqconstr_name neqconstrGrad_name ...
            Aeq Beq A B ] = CreateProblemCBs('HelloWorld',all_names, AAs ,  Cs , ineq_vars,cost_vars,true);
        disp('Done.');
        rehash;
        
        pr.objective = @(x)GenericCostFun(x,str2func(cost_name),str2func(costGrad_name)) ;
        pr.nonlcon = @(x)GenericConstrFun(x,str2func(neqconstr_name),str2func(eqconstr_name),...
            str2func(neqconstrGrad_name),str2func(eqconstrGrad_name)) ;
        pr.Aineq = A;
        pr.bineq = B;
        pr.Aeq = Aeq;
        pr.beq = Beq;
end


%% Take initial guess from simulation
disp('Filling the initial guess...');
pr.x0 = zeros(length(all_names),1);
for i=1:length(all_names)
    idx = strfind(all_names{i},'.');
    if length(idx)~=2
        continue;
    end
    name = all_names{i}(1:idx(1)-1);
    if (isempty(who(name)))
        warning(['Var ' name ' not found in workspace']);
        continue;
    end
    time = sscanf(all_names{i}(idx(2)+2:end),'%d');
%     pr.x0(i) =  eval([ name '.signals.values(' num2str(time) ')']) ;
      port = str2double(all_names{i}(idx(1)+4:idx(2)-1));
    pr.x0(i) =  eval([ name '.signals.values(' num2str(time) ',' num2str(port) ')']) ; 
  
end
disp('Done.');


%% Take initial external variables and intial state from simulation

disp('Fixing external and initial state variables');
k=0;
% fix external vars
idx = find(ex_vars);
for i=1:length(idx)
    if (strcmp(mode,'fmincon'))
        pr.Aeq(end+1,idx(i))=1;
        pr.beq(end+1) = pr.x0(idx(i));
    end
    k = k+1;
    fixed.AAs{k} = sparse(2,length(all_names),0);
    fixed.AAs{k}(1,idx(i))=1;
    fixed.Cs{k}(1) = -1;
    fixed.Cs{k}(2) = pr.x0(idx(i));
end

% fix initial state vars
idx = find(state_vars);
for i=1:length(idx)
    if (strcmp(mode,'fmincon'))
        pr.Aeq(end+1,idx(i))=1;
        pr.beq(end+1) = pr.x0(idx(i));
    end
    k = k+1;
    fixed.AAs{k} = sparse(2,length(all_names),0);
    fixed.AAs{k}(1,idx(i))=1;
    fixed.Cs{k}(1) = -1;
    fixed.Cs{k}(2) = pr.x0(idx(i));
end
disp('Done.');



%% 

switch (mode)
    case 'IPOPT_C'
        
        % profile on
        CreateIpoptCPP('test',all_names, AAs ,  Cs , ineq_vars,fixed,cost_vars);
        % profile off; profile report
        CreateIpoptDAT('test',fixed,pr.x0);

        prev = pwd;
        
         cd BLOM_Ipopt
        ! make clean
        ! make all
        cd(prev);
        ! ./BLOM_Ipopt/BLOM_NLP
        
        x = load('result.dat');
        
        
    case 'fmincon'
 
          pr.lb = [];
        pr.ub = [];
        pr.solver = 'fmincon';
        
        pr.options = optimset('GradObj','on','GradConstr','on','MaxIter',1000);
        
        disp('Running fmincon.');
        
        x = fmincon(pr);
        disp('Done.');
      
        
end
%% plot the results
subplot(211);
plot(x(strmatch('BL_System_Continuous_state.',all_names)))
ylabel('X');
subplot(212);
plot(x(strmatch('BL_System_u',all_names)))
ylabel('U');

