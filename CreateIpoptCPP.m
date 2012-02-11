function [cost_name costGrad_name eqconstr_name eqconstrGrad_name neqconstr_name neqconstrGrad_name ...
          Aeq Beq A B ]...
    =CreateIpoptCPP(name,all_names, AAs ,  Cs , ineq,fixed,cost)

[AAs , Cs] = ConvertToSingleValue(AAs,Cs);
[iAAs , iCs] = ConvertToSingleValue(ineq.AAs,ineq.Cs);
% [tAAs , tCs] = ConvertToSingleValue(fixed.AAs,fixed.Cs);

constr.AAs = {iAAs{:} AAs{:} };
constr.Cs = {iCs{:} Cs{:} };

n_eq = length(AAs);
n_ineq = length(ineq.AAs);
i_fixed = n_eq+n_ineq+1;
clear iAAs  iCs  tAAs tCs AAs Cs

for i = 1:length(all_names)
    idx_names{i} = sprintf('x[%d]',i-1);
end


cost_grad = CreateGradient(cost.A,cost.C);
cost_hessian = {};
for j=1:length(cost_grad.AAs)
    cost_hessian{j} = CreateGradient(cost_grad.AAs{j},cost_grad.Cs{j});
end


n = 100;
nnz_jac_g = 0;
Hessian.Idx = sparse(length(idx_names),length(idx_names));
Hessian.str = {};

textprogressbar('Creating Jacobian,Hessian - ');

for nn=1:n:length(constr.AAs)
    textprogressbar(nn/length(constr.AAs)*100);
    n_range = nn:min(length(constr.AAs),nn+n-1);
    constr_grad = cell(1,length(n_range));
    for i=1:length(n_range);
        constr_grad{i} = CreateGradient(constr.AAs{n_range(i)},constr.Cs{n_range(i)});
    end
    
    nnz_jac_g = CreateJacobian(name,idx_names,constr_grad,nnz_jac_g,n_range);
    
    % for i=1:length(AAs)
    %     eq_grad{i} = CreateGradient(AAs{i},Cs{i});
    % end
    %
    % for i=1:length(fixed.AAs)
    %     eq_grad{end+1} = CreateGradient(fixed.AAs{i},fixed.Cs{i});
    % end
    %
    %
    % for i=1:length(ineq.AAs)
    %     ineq_grad{i} = CreateGradient(ineq.AAs{i},ineq.Cs{i});
    % end
    
    
    hessian_valid = sparse(length(constr_grad),length(idx_names));
    hessian = cell(length(constr_grad),length(idx_names));
    for i=1:length(constr_grad)
        idx = find(constr_grad{i}.valid);
        for j=idx % 1:length(idx_names)
            if (~isempty(constr_grad{i}.Cs{j}))
                hessian{i,j} = CreateGradient(constr_grad{i}.AAs{j},constr_grad{i}.Cs{j});
                hessian_valid(i,j) = 1;
            else
                
            end
            
            
        end
    end
    
    Hessian = CreateHessian(idx_names,cost_hessian,hessian,hessian_valid,n_range,Hessian);
    
    % nnz_h_lag = CreateHessian(name,idx_names,cost_hessian,ineq_hessian,eq_hessian,ineq_hessian_valid,eq_hessian_valid);
    
    % ineq_hessian_valid = sparse(length(ineq_grad),length(idx_names),0);
    % ineq_hessian = cell(length(ineq_grad),length(idx_names));
    % for i=1:length(ineq_grad)
    %     idx = find(ineq_grad{i}.valid);
    %     for j=idx % 1:length(idx_names)
    %         if (~isempty(ineq_grad{i}.Cs{j}))
    %             ineq_hessian{i,j} = CreateGradient(ineq_grad{i}.AAs{j},ineq_grad{i}.Cs{j});
    %             ineq_hessian_valid(i,j) = 1;
    %         else
    %
    %         end
    %
    %
    %     end
    % end
    %
    % eq_hessian_valid = sparse(length(eq_grad),length(idx_names),0);
    % eq_hessian = cell(length(eq_grad),length(idx_names));
    %
    % for i=1:length(eq_grad)
    %     idx = find(eq_grad{i}.valid);
    %     for j=idx %1:length(idx_names)
    %         if (~isempty(eq_grad{i}.Cs{j}))
    %             eq_hessian{i,j} = CreateGradient(eq_grad{i}.AAs{j},eq_grad{i}.Cs{j});
    %             eq_hessian_valid(i,j) = 1;
    %         else
    %
    %         end
    %     end
    % end
    
    %%%%%%%%%%%%%%%%%%%%%%


end

textprogressbar(100);
nnz_h_lag = SaveHessian(name,Hessian);



CreateConstraints(name,idx_names,constr,i_fixed);

CreateCost(name,idx_names,cost);

CreateCostGrad(name,idx_names,cost_grad);

Create_get_nlp_info(name,length(idx_names), length(constr.AAs) ,nnz_jac_g,nnz_h_lag);
Create_get_bounds_info(name,length(idx_names),length(constr.AAs) ,n_ineq,fixed);
CreateConstructor(name,length(idx_names),length(fixed.AAs) );
 Createget_starting_point(name,length(idx_names));

function Createget_starting_point(name,n)

get_starting_point = fopen([name 'get_starting_point.cpp'],'wt');

fprintf(get_starting_point,'for (int i=0; i < %d ; i ++) { ;\n', n);
fprintf(get_starting_point,'x[i] = x0[i]; }\n', n);


fclose(get_starting_point);


function CreateConstructor(name,n,n_fixed)

Constructor = fopen([name 'Constructor.cpp'],'wt');

fprintf(Constructor,'fixed = new Number[%d] ;\n', n_fixed);
fprintf(Constructor,'x0 = new Number[%d] ;\n \n', n);

fprintf(Constructor,'FILE * f = fopen("%s","rt");\n', [name 'Fixed.dat']);
fprintf(Constructor,'if (!f) { printf("Fixed file not found"); exit(1);} \n  ');
fprintf(Constructor,'for (int i=0; i < %d ; i ++) { ;\n', n_fixed);
fprintf(Constructor,'fscanf(f,"%%lf ",&fixed[i]); }\n\n');
fprintf(Constructor,'fclose(f);\n \n');

fprintf(Constructor,' f = fopen("%s","rt");\n', [name 'X0.dat']);
fprintf(Constructor,'if (!f) { printf("X0 file not found"); exit(1);} \n  ');
fprintf(Constructor,'for (int i=0; i < %d ; i ++) { ;\n', n);
fprintf(Constructor,'fscanf(f,"%%lf ",&x0[i]); }\n\n');
fprintf(Constructor,'fclose(f);\n \n');


fclose(Constructor);

 

function Create_get_nlp_info(name,n,m,nnz_jac_g,nnz_h_lag)

get_nlp_info = fopen([name 'get_nlp_info.cpp'],'wt');

%   // The problem described in MyNLP.hpp has 2 variables, x1, & x2,
%   n = 2;
% 
%   // one equality constraint,
%   m = 1;
% 
%   // 2 nonzeros in the jacobian (one for x1, and one for x2),
%   nnz_jac_g = 2;
% 
%   // and 2 nonzeros in the hessian of the lagrangian
%   // (one in the hessian of the objective for x2,
%   //  and one in the hessian of the constraints for x1)
%   nnz_h_lag = 2;
% 
%   // We use the standard fortran index style for row/col entries
%   index_style = FORTRAN_STYLE;
  
  

fprintf(get_nlp_info,'n=%d ;\n', n);
fprintf(get_nlp_info,'m=%d ;\n', m);
fprintf(get_nlp_info,'nnz_jac_g =%d ;\n', nnz_jac_g);
fprintf(get_nlp_info,'nnz_h_lag =%d ;\n', nnz_h_lag);
fprintf(get_nlp_info,' index_style = FORTRAN_STYLE;\n');

fclose(get_nlp_info);

function Create_get_bounds_info(name,n,m,ineq_length,fixed)


get_bounds_info  = fopen([name 'get_bounds_info.cpp'],'wt');

fixed_vars = sparse(1,n,length(fixed.AAs));
for i=1:length(fixed.AAs)
    fixed_vars(find(fixed.AAs{i}(1,:)))=i;
end

for i=1:n
    if (fixed_vars(i) ~= 0)
         fprintf(get_bounds_info,'x_l[%d] = fixed[%d];\n', i-1,fixed_vars(i)-1);
         fprintf(get_bounds_info,'x_u[%d] = fixed[%d];\n', i-1,fixed_vars(i)-1);
    else
        fprintf(get_bounds_info,'x_l[%d] = -1.0e19;\n', i-1);
        fprintf(get_bounds_info,'x_u[%d] = +1.0e19;\n', i-1);
    end
end

for i=1:ineq_length
    fprintf(get_bounds_info,'g_l[%d] = -1.0e19;\n', i-1);
    fprintf(get_bounds_info,'g_u[%d] = 0;\n', i-1);
end


for i=ineq_length+1:m
    fprintf(get_bounds_info,'g_l[%d] = 0;\n', i-1);
    fprintf(get_bounds_info,'g_u[%d] = 0;\n', i-1);
end

fclose(get_bounds_info);

function CreateCost(name,names,cost)

feval_f = fopen([name 'eval_f.cpp'],'wt');



str = CreatePolyFromMatrix(cost.A,cost.C,names,'high','C');


fprintf(feval_f,'obj_value=%s ;\n', str);

fclose(feval_f);

function CreateCostGrad(name,names,cost_grad)

feval_grad_f = fopen([name 'eval_grad_f.cpp'],'wt');

for i=1:length(cost_grad.AAs) 
        str = '0';
        
            if  ~isempty(cost_grad.AAs) && ~isempty(cost_grad.AAs{i})

                str = CreatePolyFromMatrix(cost_grad.AAs{i},cost_grad.Cs{i},names,'high','C');
            end
        
        
          fprintf(feval_grad_f,'grad_f[%d]=%s ;\n', i-1,str);

end

fclose(feval_grad_f);


function CreateConstraints(name,names,constr,i_fixed)

feval_g = fopen([name 'eval_g.cpp'],'wt');

for i=1:length(constr.AAs) 
        str = '0';
        
         if (i < i_fixed)
            if  ~isempty(constr.AAs) && ~isempty(constr.AAs{i})

                str = CreatePolyFromMatrix(constr.AAs{i},constr.Cs{i},names,'high','C');
            end
         else
            k = i - i_fixed+1; 
            if ~isempty(constr.AAs) && ~isempty(constr.AAs{i})
                [jA ] = find(constr.AAs{i}(1,:));
                
                str = [ num2str(constr.Cs{i}(1)) '*'  names{jA}  ' + fixed[' num2str(k-1) ']'];
            end
            
        end
        
        
          fprintf(feval_g,'g[%d]=%s ;\n', i-1,str);

end

fclose(feval_g);


function Hessian = CreateHessian(names,cost_hessian,hessian,hessian_valid,n_range,Hessian)

% nnz_h_lag=  0 ;

% if nnz_h_lag == 0
%     feval_h_val = fopen([name 'eval_h_val.cpp'],'wt');
%     feval_h_idx = fopen([name 'eval_h_idx.cpp'],'wt');
% else
%     feval_h_val = fopen([name 'eval_h_val.cpp'],'at');
%     feval_h_idx = fopen([name 'eval_h_idx.cpp'],'at');
% end

if (~isempty(cost_hessian)&& n_range(1) == 1) % first time only
    for i= 1:length(names)
        if(isempty(cost_hessian{i}.valid))
            continue;
        end
        idx = find(cost_hessian{i}.valid(1:i));
        if (~isempty(idx))
            for j=idx %1:i
                tmp = CreatePolyFromMatrix(cost_hessian{i}.AAs{j},cost_hessian{i}.Cs{j},names,'high','C');
                str = [  'obj_factor*(' tmp ')' ];
                Hessian = AddToHessian(Hessian,i,j,str);
            end
        end
    end   
end

for k=1:length(n_range)
        idx_i = find(hessian_valid(k,:));    
    for i=idx_i % 1:length(names)
         idx = find(hessian{k,i}.valid(1:i));
        if (~isempty(idx))
            for j=idx %1:i
                tmp = CreatePolyFromMatrix(hessian{k,i}.AAs{j},hessian{k,i}.Cs{j},names,'high','C');
                if (~strcmp(tmp,'0'))
                    str = [ '+lambda[' num2str(n_range(k)-1) ']*(' tmp ')' ];
                    Hessian = AddToHessian(Hessian,i,j,str);
                end
            end
        end
    end
end 
% 
% for i=1:length(n_range) % 1:length(names)
%         ineqK = zeros(1,size(hessian,1));
%         idx = find(hessian_valid(:,i)');
%         for k = idx %1:size(ineq_hessian,1)
%             if  ~isempty(hessian{k,i}.AAs) 
%                     ineqK(k) = 1;
%             end
%         end
%         ineqK_idx{i} = find(ineqK);
% end
% 
% % for i=n_range % 1:length(names)
% %         eqK = zeros(1,size(eq_hessian,1));
% %         idx = find(eq_hessian_valid(:,i)');
% %         for k = idx %1:size(eq_hessian,1)
% %             if  ~isempty(eq_hessian{k,i}.AAs) 
% %                     eqK(k) = 1;
% %             end
% %         end
% %         eqK_idx{i} = find(eqK);
% % end
% 
% % if  ~isempty(cost_hessian)
%     for i=1:length(names)
%         valid = sparse(1,length(names));
%         if (~isempty(cost_hessian{i}.valid))
%             valid = cost_hessian{i}.valid;
%         end
%         range = ineqK_idx{i};
%         if ~isempty(range)
%             for k = range
%                 valid = valid |  hessian{k,i}.valid;
%             end
%         end
%         
%         
%         idx = find(valid(1:i));
%         if (~isempty(idx))
%             for j=idx %1:i
%                 nz = 0;
%                 str = '';
%                 if  ~isempty(cost_hessian{i}.AAs) &&  ~isempty(cost_hessian{i}.Cs{j})
%                     nz = 1;
%                     str = CreatePolyFromMatrix(cost_hessian{i}.AAs{j},cost_hessian{i}.Cs{j},names);
%                     str = [ 'obj_factor*(' str ')'];
%                 end
%                 
%                 range = ineqK_idx{i};
%                 if ~isempty(range)
%                     for k = range
%                         if    ~isempty(ineq_hessian{k,i}.Cs{j})
%                             nz = 1;
%                             tmp = CreatePolyFromMatrix(ineq_hessian{k,i}.AAs{j},ineq_hessian{k,i}.Cs{j},names);
%                             str = [str  '+lambda[' num2str(k-1) ']*(' tmp ')' ];
%                         end
%                     end
%                 end
%                 
%                 range = eqK_idx{i};
%                 if ~isempty(range);
%                     for k =  range
%                         if   ~isempty(eq_hessian{k,i}.Cs{j})
%                             nz = 1;
%                             tmp = CreatePolyFromMatrix(eq_hessian{k,i}.AAs{j},eq_hessian{k,i}.Cs{j},names);
%                             str = [str  '+lambda[' num2str(k-1+size(ineq_hessian,1)) ']*(' tmp ')' ];
%                         end
%                     end
%                 end
%                 if (nz == 1)
%                     Hessian = AddToHessian(Hessian,n_range(i),j,str);
% %                     nnz_h_lag = nnz_h_lag + 1;
% %                     fprintf(feval_h_val,'values[%d]=%s ;\n', nnz_h_lag-1,str);
% %                     fprintf(feval_h_idx,'iRow[%d]=%d ;\n', nnz_h_lag-1,i);
% %                     fprintf(feval_h_idx,'jCol[%d]=%d ;\n', nnz_h_lag-1,j);
%                 end
%             end
%         end
%     end
% end
% fclose(feval_h_val);
% fclose(feval_h_idx);

function nnz_jac_g = CreateJacobian(name,names,constr_grad,nnz_jac_g,n_range)

% nnz_jac_g=  0 ;

if (nnz_jac_g ==  0) % first time reset the file
    feval_jac_g_val = fopen([name 'eval_jac_g_val.cpp'],'wt');
    feval_jac_g_idx = fopen([name 'eval_jac_g_idx.cpp'],'wt');
else
    feval_jac_g_val = fopen([name 'eval_jac_g_val.cpp'],'at');
    feval_jac_g_idx = fopen([name 'eval_jac_g_idx.cpp'],'at');
end

for i=1:length(n_range)
        idx = find(constr_grad{i}.valid); 
    
    for j=idx %1:length(names)
        nz = 0;
%         str = '';
        
            if  ~isempty(constr_grad{i}.AAs) && ~isempty(constr_grad{i}.Cs{j})
                
                str = CreatePolyFromMatrix(constr_grad{i}.AAs{j},constr_grad{i}.Cs{j},names,'high','C');
                if (~strcmp(str,'0'))
                    nz = 1;
                end
            end
        
        
        if (nz == 1)
          nnz_jac_g = nnz_jac_g + 1;
          fprintf(feval_jac_g_val,'values[%d]=%s ;\n', nnz_jac_g-1,str);
          fprintf(feval_jac_g_idx,'iRow[%d]=%d ;\n', nnz_jac_g-1,n_range(i));
          fprintf(feval_jac_g_idx,'jCol[%d]=%d ;\n', nnz_jac_g-1,j);
        end
    end
end

fclose(feval_jac_g_val);
fclose(feval_jac_g_idx);



function [ Matrix , B,   nlAAs , nlCs ] = ExtractLinearPart(AAs,Cs)

l=1;
nl=1;
Matrix = [];
B = [];
nlAAs = {};
nlCs  = {};
for i = 1:length(AAs)
    if IsLinear(AAs{i})
        for j=1:size(Cs{i},1)
            Matrix(l,:) = Cs{i}(j,:)*AAs{i};
            idx = find(sum(AAs{i}')==0);% find constant terms
            if isempty(idx)
                B(l) = 0;
            else
                B(l) = -sum(Cs{i}(j,idx));
            end

            l=l+1;
        end
       else
           nlAAs{nl} = AAs{i};
           nlCs{nl} = Cs{i};
           nl = nl+1;
       end
           
    
end

Matrix = sparse(Matrix);

function is_linear = IsLinear(A)
if (max(sum(A'))<=1)&& (min(sum(A')) > -1 )
    is_linear = true;
else
    is_linear = false;
end
   
function [AAs , Cs] = ConvertToSingleValue(inAAs,inCs)
% convert to single value functions
k=1;
for i=1:length(inAAs)
    for j=1:size(inCs{i},1)
        AAs{k} = inAAs{i};
        [I,J,s] = find(abs(inCs{i}(j,:))>1e-10);
        [m,n] = size(inCs{i}(j,:));
        Cs{k} = sparse(I,J,inCs{i}(j,J),m,n); %inCs{i}(j,:);
        k=k+1;
    end
end

function Hessian = AddToHessian(Hessian,i,j,str)

if (Hessian.Idx(i,j) == 0) % new entry
    Hessian.Idx(i,j) = length(Hessian.str)+1;
    Hessian.str{Hessian.Idx(i,j)} = '';
end
Hessian.str{Hessian.Idx(i,j)} = [Hessian.str{Hessian.Idx(i,j)} str];
    
    
    
    
function nnz_h_lag = SaveHessian(name,Hessian)


feval_h_val = fopen([name 'eval_h_val.cpp'],'wt');
feval_h_idx = fopen([name 'eval_h_idx.cpp'],'wt');


[I,J,s] = find(Hessian.Idx);
nnz_h_lag = length(I);
for i=1:length(I);
    fprintf(feval_h_val,'values[%d]=%s ;\n', i-1,Hessian.str{s(i)});
    fprintf(feval_h_idx,'iRow[%d]=%d ;\n', i-1,I(i));
    fprintf(feval_h_idx,'jCol[%d]=%d ;\n', i-1,J(i));
end
fclose(feval_h_val);
fclose(feval_h_idx);
