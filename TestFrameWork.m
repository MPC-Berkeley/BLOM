clear
clc

AR_A = [ 0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1];
AR_C =[ 0.001630889584137   0.000960687393672   0.000526525284885   0.000097391141316   0.000179393992191  -0.000171463047318 ...
  -0.000107550015753  -0.000181196558137   0.000181946649478   3.267693499582606   1.183669483246657  -0.197967751006581   -1.000000000000000 ];

% Tr = k*Toa + (1-k)Tz
% Toa - 1 ; k - 2 ; Tz - 3; Tr - 4
mix_A = [1 1 0 0 ; 
         0 0 1 0 ;
         0 1 1 0 ;
         0 0 0 1 ];
mix_C = [1 1 -1 -1];     

% Out = x1 - x2 
% x1 - 1 ; x2 - 2 ; Out - 3;
delta_A = [1 0 0 ; 
         0 1 0 ;
         0 0 1];
delta_C = [1 -1 -1 ];     

open('BLOM_ex');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[all_names, AAs ,  Cs , ineq_vars ,cost_vars,in_vars,state_vars ] = ExtractModel(10);
sim('BLOM_ex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
code = PrintProblem(all_names, AAs ,  Cs ,ineq_vars,cost_vars);
for i=1:length(code)
    disp(code{i});
end

[cost_name costGrad_name eqconstr_name eqconstrGrad_name neqconstr_name neqconstrGrad_name ...
          Aeq Beq A B ] = CreateProblemCBs('test',all_names, AAs ,  Cs , ineq_vars,cost_vars,true);

rehash;

pr.objective = @(x)GenericCostFun(x,str2func(cost_name),str2func(costGrad_name)) ;
pr.nonlcon = @(x)GenericConstrFun(x,str2func(neqconstr_name),str2func(eqconstr_name),...
                                    str2func(neqconstrGrad_name),str2func(eqconstrGrad_name)) ;
pr.Aineq = A;
pr.bineq = B;
pr.Aeq = Aeq;
pr.beq = Beq;

pr.lb = [];
pr.ub = [];
pr.solver = 'fmincon';

pr.options = optimset('GradObj','on','GradConstr','on','MaxIter',1000);

pr.x0 = zeros(length(all_names),1);
for i=1:length(all_names)
    idx = strfind(all_names{i},'.');
    if length(idx)~=2
        continue;
    end
    name = all_names{i}(1:idx(1)-1);
    if (isempty(who(name)))
        name
        continue;
    end
    time = str2double(all_names{i}(idx(2)+2:end));
    pr.x0(i) =  eval([ name '.signals.values(' num2str(time) ')']) ;
    
end

% fix TOA and Tref
k=0;
idx = strmatch('AHU_VAV_TOA.',all_names);
for i=1:length(idx)
    pr.Aeq(end+1,idx(i))=1;
    pr.beq(end+1) = pr.x0(idx(i));
    k = k+1;
    fixed.AAs{k} = sparse(2,length(all_names),0);
    fixed.AAs{k}(1,idx(i))=1;
    fixed.Cs{k}(1) = -1;
    fixed.Cs{k}(2) = pr.x0(idx(i));
    
end
idx = strmatch('AHU_VAV_Tref.',all_names);
for i=1:length(idx)
    pr.Aeq(end+1,idx(i))=1;
    pr.beq(end+1) = pr.x0(idx(i));
    k = k+1;
    fixed.AAs{k} = sparse(2,length(all_names),0);
    fixed.AAs{k}(1,idx(i))=1;
    fixed.Cs{k}(1) = -1;
    fixed.Cs{k}(2) = pr.x0(idx(i));
end

% fix initial state vars
idx = find(state_vars);
for i=1:length(idx)
    pr.Aeq(end+1,idx(i))=1;
    pr.beq(end+1) = pr.x0(idx(i));
    k = k+1;
    fixed.AAs{k} = sparse(2,length(all_names),0);
    fixed.AAs{k}(1,idx(i))=1;
    fixed.Cs{k}(1) = -1;
    fixed.Cs{k}(2) = pr.x0(idx(i));
end

mode = 'fmincon' ;

switch (mode)
    case 'IPOPT_C'
        
        % profile on
        CreateIpoptCPP('test',all_names, AAs ,  Cs , ineq_vars,fixed,cost_vars);
        % profile off; profile report
        CreateIpoptDAT('test',fixed,pr.x0);

        prev = pwd;
        
         cd ~/workspace/BLOM_Ipopt
        ! make clean
        ! make all
        
        ! ./BLOM_NLP
        
        x = load('result.dat');
        
        cd(prev);
    case 'fmincon'
 
        
        x = fmincon(pr);
        
end
%% plot the results
subplot(211);
plot(x(strmatch('AHU_VAV_Tin.',all_names))-273)
ylabel('Tin [^oC]');
subplot(212);
plot(x(strmatch('AHU_VAV_ARmodel.Out1',all_names))-273)
ylabel('T zone [^oC]');

