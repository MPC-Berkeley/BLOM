// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: MyNLP.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "BLOM_NLP.hpp"

#include <cassert>



using namespace Ipopt;

/* Constructor. */
MyNLP::MyNLP()
{
	fixed = 0;
	x0 = 0;
	
    FILE * fdat = fopen("params.dat", "rt");
    if (!fdat) { printf("Error: params.dat file not found\n"); exit(1);}
    fscanf(fdat," %d %d %d %d %d %d ",&m_n,&m_n_fixed,&m_m,&m_nnz_jac,&m_nnz_hessian,&m_m_ineq_constrs);
    fclose(fdat);

    
    fixed = new Number[m_n_fixed] ;
    x0 = new Number[m_n] ;
    
    FILE * f = fopen("testFixed.dat", "rt");
    if (!f) { printf("Error: Fixed file not found\n"); exit(1);}
    for (int i=0; i < m_n_fixed ; i ++) { ;
    fscanf(f, "%lf ", &fixed[i]); }
    
    fclose(f);
    
    f = fopen("testX0.dat", "rt");
    if (!f) { printf("Error: X0 file not found\n"); exit(1);}
    for (int i=0; i < m_n ; i ++) { ;
    fscanf(f, "%lf ", &x0[i]); }
    
    fclose(f);
    

    
	ReadAandC();
//	printf("fixed[0] = %lf \n fixed[51] = %lf \n",fixed[0],fixed[51]);  	
//	printf("x0[0] = %lf \n x0[191] = %lf \n",x0[0],x0[191]);  
//  printf("\n Constructor ");
}


MyNLP::~MyNLP()
{
	delete [] fixed;
	delete [] x0;
	

}

bool MyNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	
	
// n=28400 ;
// m=29086 ;
// nnz_jac_g =120364 ;
// nnz_h_lag =3640 ;
// m_m_ineq_constrs=5940 ;
//  // The problem described in MyNLP.hpp has 2 variables, x1, & x2,
//  n = 2;
//
//  // one equality constraint,
//  m = 1;
//
//  // 2 nonzeros in the jacobian (one for x1, and one for x2),
//  nnz_jac_g = 2;
//
//  // and 2 nonzeros in the hessian of the lagrangian
//  // (one in the hessian of the objective for x2,
//  //  and one in the hessian of the constraints for x1)
//  nnz_h_lag = 2;
//
//  // We use the standard fortran index style for row/col entries
//  index_style = FORTRAN_STYLE;
//  printf("\n get_nlp_info, %d,%d ",m,n);
            
    nnz_jac_g = m_nnz_jac;
    nnz_h_lag = m_nnz_hessian;
    n = m_n;
    m = m_m;
     index_style = FORTRAN_STYLE;
     
  return true;
}

bool MyNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
 
    GetBounds_info(n,  x_l,  x_u,m,g_l, g_u);
//  printf("\n get_bounds_info, %d,%d ",m,n);
  return true; 
}

bool MyNLP::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
    // your own NLP, you can provide starting values for the others if
    // you wish.
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);
//  printf("\n get_starting_point");
    
    for (int i=0; i < m_n ; i ++) {
        x[i] = x0[i];
    }
    
//
//Number g[292];
//	eval_g(n,  x, true, 292,  g);
//
//	for (int j =0 ; j < 292; j ++)
//		printf(" %f \n", g[j]);
		

  return true;
}

bool MyNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
//  printf("\n eval_f");
  // return the value of the objective function
  obj_value =  CalcValue(m_A,m_C,0,x);
  return true;
}

bool MyNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // return the gradient of the objective function grad_{x} f(x)
//  printf("\n eval_grad_f");
  // grad_{x1} f(x): x1 is not in the objective
    for (int i=0; i < m_n ; i ++)
    {
       grad_f[i] =  CalcDerValue(m_A,m_C,0 ,x,i);
    }
  return true;
}

bool MyNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
//	 printf("\n eval_g");
  // return the value of the constraints: g(x)
    CalcConstraintsValues(g,x) ;

  return true;
}

bool MyNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian of the constraints

      GetJacobianStruct(iRow,jCol) ;
//       eval_jac_g_idx(n, x,  new_x, m, nele_jac, iRow, jCol, values);
  }
  else {
    // return the values of the jacobian of the constraints
//#include "testeval_jac_g_val.cpp"
      CalcJacobianValues(values,x) ;
//       eval_jac_g_val(n,x,  new_x,m, nele_jac,iRow, jCol,values);

  }
//  printf("\n eval_jac_g");

  return true;
}

bool MyNLP::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
  if (values == NULL) {

	GetHessianStruct(iRow,jCol) ;


    // Note: off-diagonal elements are zero for this problem
  }
  else {
    // return the values

	CalcHessianValues(values,x,obj_factor,lambda) ;

    // Note: off-diagonal elements are zero for this problem
  }

//  printf("\n eval_h");
  return true;
}

void MyNLP::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
			      const IpoptData* ip_data,
			      IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution. Since the solution is displayed to the console,
  // we currently do nothing here.
  
  FILE * f = fopen("result.dat","wt");
  for (int i = 0 ; i < n ; i ++ )
  	fprintf(f,"%16.15e \n",x[i]);
  
  fclose(f);
  
}
