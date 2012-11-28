#include "BLOM_NLP.hpp"
#include <string.h>
#include <math.h>

__inline double EvalSpecial(double tmp,const double p, const double * x,const int idx )
{
    if (p == BLOM_TYPE_EXP) //exp
    {
        tmp *= exp(x[idx]);
        if (std::isinf(tmp))
            printf("BLOM_NLP_Sparse:exp overflow,j=%d, %f %f \n",idx,x[idx],exp(x[idx]));
        if (std::isnan(tmp))
            printf("BLOM_NLP_Sparse:exp nan, j=%d, %f %f \n",idx,x[idx],exp(x[idx]));
    }
    else if (p == BLOM_TYPE_LOG) //log
    {
        tmp *= log(x[idx]);
    }
    else if (p == BLOM_TYPE_SIN) //sin
    {
        tmp *= sin(x[idx]);
    }
    else if (p == BLOM_TYPE_COS) //cos
    {
        tmp *= cos(x[idx]);
    }
    else if (p == BLOM_TYPE_TANH) //tanh
    {
        tmp *= tanh(x[idx]);
    }
    else if (p == BLOM_TYPE_ATAN) //atan
    {
        tmp *= atan(x[idx]);
    }
    else
    {
        printf("BLOM_NLP_Sparse:EvalSpecial: unknown special code %g ! \n",p);
    }
    
    return tmp;
}


double MyNLP::  CalcValue(CompRow_Mat_double& A,CompRow_Mat_double&  C,int f,const double x[])
{
	int row = C.row_ptr(f);
	int row_end = C.row_ptr(f+1);;
	
	double val = 0;
	
	for (int i =row ; i < row_end ; i ++)
	{
			int A_row = A.row_ptr(C.col_ind(i));
			int A_row_end = A.row_ptr(C.col_ind(i)+1);;
			
			double tmp = C.val(i);
			 for (int j=A_row ; j < A_row_end ; j ++ )
			 {
			 	double p = A.val(j);
                if (p < BLOM_TYPE_EXP ) // The main case - powers: integer and non-integer
                {
                    if (floor(p)==p && p >= 0 ) // integer and non negative
                    {
                        for (; p--;  tmp*=x[A.col_ind(j)]);
                    }
                    else  // non-integer power
                    {
                        tmp *=pow(x[A.col_ind(j)],p);
                        if (std::isnan(tmp))
                            printf("BLOM_NLP_Sparse:tmp is nan, j=%d, var  %d, %f %f %f\n", j,A.col_ind(j), x[A.col_ind(j)], p,pow(x[A.col_ind(j)],p));
                    }
                }
			 	else 
                {
                    tmp = EvalSpecial(tmp,p,x,A.col_ind(j));
                }
                if (std::isnan(tmp))
                    printf("BLOM_NLP_Sparse:tmp is nan, j=%d,var  %d, %f %f \n", j,A.col_ind(j), x[A.col_ind(j)], p);

			 }
			 
			 val += tmp;
		
	}
	
    if (std::isinf(val))
        printf("BLOM_NLP_Sparse:val overflow, f=%d\n", f);
    
    if (std::isnan(val))
        printf("BLOM_NLP_Sparse:val is nan number, f=%d\n", f);    
    
	return val;
}

double MyNLP::  CalcDerValue(CompRow_Mat_double& A,CompRow_Mat_double&  C,int f,const Number x[],int dvar)
{
	int row = C.row_ptr(f);
	int row_end = C.row_ptr(f+1);;
	
	double val = 0;
	
	for (int i =row ; i < row_end ; i ++) // Loop on sparse indices of the row - C coefficients
	{
            // fetch the associated row from A
			int A_row = A.row_ptr(C.col_ind(i)); // begining of the current A row
			int A_row_end = A.row_ptr(C.col_ind(i)+1); // end of the current A row
			
			double tmp = C.val(i); // initialize with coefficient value
            // Check if dvar is present in this term
			bool found=false;
			 for (int j=A_row ; j < A_row_end ; j ++ )
			 {
			 	if (A.col_ind(j) == dvar && A.val(j) !=0)
			 	{
			 		found = true;
			 		break;
			 	}
			 }
			 if (!found) // if this term does not include the dvar, the derivative is zero.
			 {
                continue;
             }
            
            // evaluate all term multipliers
            for (int j=A_row ; j < A_row_end ; j ++ )
            {
                double p = A.val(j);
                if (A.col_ind(j) == dvar) // special treatment for dvar
                {
                    if(p < BLOM_TYPE_EXP)
                    {
                        tmp *= p --;
                    }
                    else if (p == BLOM_TYPE_EXP) //exp
                    {
                        // do nothing for exp
                        
                    }
                    else if (p ==  BLOM_TYPE_LOG) //log
                    {
                        p=-1;
                        
                    }
                    else if (p ==  BLOM_TYPE_SIN) //sin
                    {
                        p = BLOM_TYPE_COS;
                    }
                    else if (p == BLOM_TYPE_COS) //cos
                    {
                        tmp *= -1;
                        p = BLOM_TYPE_SIN;
                    }
                    else if (p == BLOM_TYPE_TANH) //tanh
                    { // special treatment: evaluate it here:
                        double th = tanh(x[A.col_ind(j)]);
                        tmp *= 1 - th*th;
                        continue; // continue to the next variable
                    }
                    else if (p == BLOM_TYPE_ATAN) //atan
                    { // special treatment: evaluate it here:
                        tmp *= 1/(x[A.col_ind(j)]*x[A.col_ind(j)] + 1);
                        continue; // continue to the next variable
                    }
                    else
                    {
                        printf("BLOM_NLP_Sparse:EvalSpecial: unknown special code %g ! \n",p);
                    }
                    
                    
                }

                
                
                if (p < BLOM_TYPE_EXP ) // The main case - powers: integer and non-integer
                {
                    if (floor(p)==p && p >= 0 ) // integer and non negative
                    {
                        for (; p--;  tmp*=x[A.col_ind(j)]);
                    }
                    else  // non-integer power
                    {
                        tmp *=pow(x[A.col_ind(j)],p);
                        if (std::isnan(tmp))
                            printf("BLOM_NLP_Sparse:CalcDerValue: tmp is nan, j=%d, var  %d, %f %f %f\n", j,A.col_ind(j), x[A.col_ind(j)], p,pow(x[A.col_ind(j)],p));
                    }
                }
                else
                {
                    tmp = EvalSpecial(tmp,p,x,A.col_ind(j));
                }
			 	
			 }
			 
			 val += tmp;
		
	}
    
    if (std::isinf(val))
        printf("BLOM_NLP_Sparse:der val overflow, f=%d\n", f);

    if (std::isnan(val))
        printf("BLOM_NLP_Sparse:der val is not normal number, f=%d\n", f);
    
	return val;
}

double MyNLP::  CalcDoubleDerValue(CompRow_Mat_double& A,CompRow_Mat_double&  C,int f,const Number x[],int dvar1,int dvar2)
{
	int row = C.row_ptr(f);
	int row_end = C.row_ptr(f+1);;
	
	double val = 0;
	
	for (int i =row ; i < row_end ; i ++)
	{
			int A_row = A.row_ptr(C.col_ind(i));
			int A_row_end = A.row_ptr(C.col_ind(i)+1);;
			
			double tmp = C.val(i);
			bool found=false;
			 for (int j=A_row ; j < A_row_end ; j ++ )
			 {
			 	if (A.col_ind(j) == dvar1 && A.val(j) !=0)
			 	{
			 		found = true;
			 		break;
			 	}
			 }
			 if (!found)
			 {
			 		continue;
			 }
			 
			 found=false;		 
			 for (int j=A_row ; j < A_row_end ; j ++ )
			 {
			 	if (A.col_ind(j) == dvar2 && A.val(j) !=0)
			 	{
			 		found = true;
			 		break;
			 	}
			 }
			 if (!found)
			 {
			 		continue;
             }
             
             for (int j=A_row ; j < A_row_end ; j ++ )
             {
                 double p = A.val(j);
                 bool   eval = true;
                 
                 for (int k=0, dvar = dvar1; k<2; k++, dvar = dvar2) // check both derivatives
                 { // Apply the derivative each time dvar is found
                     if (A.col_ind(j) == dvar)
                     {
                         if(p < BLOM_TYPE_EXP)
                         {
                             tmp *= p --;
                         }
                         else if (p == BLOM_TYPE_EXP) //exp
                         {
                             // do nothing for exp
                         }
                         else if (p == BLOM_TYPE_LOG) //log
                         {
                             p=-1;
                         }
                         else if (p == BLOM_TYPE_SIN) //sin
                         {
                             p = BLOM_TYPE_COS;
                         }
                         else if (p == BLOM_TYPE_COS) //cos
                         {
                             tmp *= -1;
                             p = BLOM_TYPE_SIN;
                         }
                         else if (p == BLOM_TYPE_TANH) //tanh
                         { // special treatment: evaluate it here:
                             double th = tanh(x[A.col_ind(j)]);
                             if (dvar1==dvar2) // 2nd derivative of tanh
                             {
                                 // the first time we are in dvar1
                                 tmp *= -2*th*(1 - th*th);
                                 eval = false; // continue to the next variable
                                 break; // do not check the dvar2
                             }
                             else
                             {
                                 tmp *= 1 - th*th;
                                 eval = false; // continue to the next variable
                             }
                         }
                         else if (p == BLOM_TYPE_ATAN) //atan
                         { // special treatment: evaluate it here:
                             double xj = x[A.col_ind(j)];
                             if (dvar1==dvar2) // 2nd derivative of atan
                             {
                                 // the first time we are in dvar1
                                 tmp *= -2*xj/((xj*xj + 1)*(xj*xj + 1));
                                 eval = false; // continue to the next variable
                                 break; // do not check the dvar2
                             }
                             else
                             {
                                 tmp *= 1/(xj*xj + 1);
                                 eval = false; // continue to the next variable
                             }
                         }
                         else
                         {
                             printf("BLOM_NLP_Sparse:EvalSpecial: unknown special code %g ! \n",p);
                         }
                     }
                 }
                 
                 
                 if (eval) { // evalutate if was not evaluted before, such as tanh.
                     if (p < BLOM_TYPE_EXP ) // The main case - powers: integer and non-integer
                     {
                         if (floor(p)==p && p >= 0 ) // integer and non negative
                         {
                             for (; p--;  tmp*=x[A.col_ind(j)]);
                         }
                         else  // non-integer power
                         {
                             tmp *=pow(x[A.col_ind(j)],p);
                             if (std::isnan(tmp))
                                 printf("BLOM_NLP_Sparse:CalcDoubleDerValue: tmp is nan, j=%d, var %d, %f %f %f\n", j,A.col_ind(j), x[A.col_ind(j)], p,pow(x[A.col_ind(j)],p));
                         }
                     }
                     else
                     {
                         tmp = EvalSpecial(tmp,p,x,A.col_ind(j));
                     }
                 }
                 
             }
             
             val += tmp;
             
    }
    
    if (std::isinf(val))
        printf("BLOM_NLP_Sparse:dder val overflow, f=%d\n", f);
    
    return val;
}


void MyNLP::ReadAandC()
{
//	CompRow_Mat_double A;
//	CompRow_Mat_double C;
//	double x[] = { 1, 2 ,7, 4, 5, 6};
	
	 readtxtfile_mat("A.txt", &m_A);
	 readtxtfile_mat("C.txt", &m_C);
	 
	 readtxtfile_mat("JacobianStruct.txt", &m_JacobianStruct);
	 readtxtfile_mat("HessianStruct.txt", &m_HessianStruct);
     readtxtfile_mat("LambdaStruct.txt", &m_LambdaStruct);
     readtxtfile_mat("FixedStruct.txt", &m_FixedStruct);
     
	 // printf("\n\n HERE \n\n");
	// // printf("\n\n %f %f %f \n\n", CalcValue( m_A, m_C,0, x),CalcDerValue( A, C,0, x,1),CalcDoubleDerValue( A, C,0, x,2,1));
}


void MyNLP::  CalcConstraintsValues(Number values[],const Number x[]) 
{
	for (int i=1; i < m_C.dim(0); i ++ )
	{
        // printf("%d\n",i);
		values[i-1] = CalcValue(m_A,m_C,i,x);
	}
}
 
void MyNLP::CalcJacobianValues(Number values[],const Number x[]) 
{
	for (int i=0; i<m_nnz_jac ; i ++)
	{
		values[i] = CalcDerValue(m_A,m_C,m_JacobianStruct.row_ind(i)+1 ,x,m_JacobianStruct.col_ind(i));
	}
	
}


void MyNLP::GetJacobianStruct(int iRow[],int jCol[]) 
{
            // printf("GetJacobianStruct  \n");

	for (int i=0; i<m_nnz_jac ; i ++)
	{
                          // printf("GetJacobianStruct %d\n",i);

		iRow[i] = m_JacobianStruct.row_ind(i)+1;
		jCol[i] = m_JacobianStruct.col_ind(i)+1;
	}
                    // printf("GetJacobianStruct end \n");

}



void MyNLP::CalcHessianValues(Number values[],const Number x[],Number obj_factor,const Number* lambda)
 {
        // printf("CalcHessianValues \n");

    
    memset(values, 0, sizeof(double)*m_nnz_hessian);

    if (obj_factor) {
        for (int i=m_LambdaStruct.row_ptr(0); i<m_LambdaStruct.row_ptr(1) ; i ++) {
            // printf("o %d %d\n",m_LambdaStruct.col_ind(i),i);
            values[m_LambdaStruct.col_ind(i)] = obj_factor*CalcDoubleDerValue(m_A, m_C, 0, x, m_HessianStruct.row_ind(m_LambdaStruct.col_ind(i)), m_HessianStruct.col_ind(m_LambdaStruct.col_ind(i)));
        }
    }
    
    for (int j =0 ; j<m_JacobianStruct.dim(0) ; j ++) // for all constraints
    {
                // printf(" j %d\n",j);

        if (lambda[j]) {
            
            for (int i=m_LambdaStruct.row_ptr(j+1); i<m_LambdaStruct.row_ptr(j+2) ; i ++) 
            {
            // printf("%d %d %d\n",j,m_LambdaStruct.col_ind(i),i);
                values[m_LambdaStruct.col_ind(i)] += lambda[j]*CalcDoubleDerValue(m_A, m_C, j+1, x, m_HessianStruct.row_ind(m_LambdaStruct.col_ind(i)), m_HessianStruct.col_ind(m_LambdaStruct.col_ind(i)));
            }
        }
    }
        // printf("CalcHessianValues end \n");

}

void MyNLP::GetHessianStruct(int iRow[],int jCol[]) 
{
            // printf("GetHessianStruct \n");

	for (int i=0; i<m_nnz_hessian ; i ++)
	{
          // printf("GetHessianStruct %d\n",i);
		iRow[i] = m_HessianStruct.row_ind(i)+1;
		jCol[i] = m_HessianStruct.col_ind(i)+1;
	}

        // printf("GetHessianStruct end \n");

}

void MyNLP::GetBounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
    for(int i=0; i < m_n ; i ++)
    {
        x_l[i] =-1.0e19;
        x_u[i] =+1.0e19;
    }
   
    // Add variable boundaries for trivial linear constraints
    for (int i=1; i <= m_m  ; i ++) // Do inequalities first, because if the same variable appears in equality and inequality, the equality shoud prevail.
    {
//        printf("i=%d ,%d \n ",i,m_C.row_ptr(i+1));
        if  ((m_C.row_ptr(i+1) - m_C.row_ptr(i))<=2)   // 2 or less terms in constraint
        {
            int linear_var = -1;
            int const_term = -1;
            double linear_coeff = 0;
            double const_coeff  = 0;
            bool skip = false;
            
            
            for (int j=m_C.row_ptr(i) ; j < m_C.row_ptr(i+1) ; j ++)
            {
                int A_row = m_A.row_ptr(m_C.col_ind(j)); // begining of the current A row
                int A_row_end = m_A.row_ptr(m_C.col_ind(j)+1); // end of the current A row
                
//                 printf("%d , ",A_row_end-A_row );
                if (A_row_end-A_row > 1) // more than one variable in this term
                {
                    skip = true;
                    break;
                }
                
                
                if (A_row_end-A_row == 1)// The row of A has only one term.
                {
                     if (m_A.val(A_row) == 1 && linear_var == -1) // linear term for the first time
                     {
                         linear_var = m_A.col_ind(A_row);
                         linear_coeff = m_C.val(j);
                     }
                     else if (m_A.val(A_row) == 0 && const_term == -1) { // constant term for the first time
                         const_term = m_A.col_ind(A_row);
                         const_coeff = m_C.val(j);
                     }
                     else {
                         skip = true;
                         //printf("skipped, const %d, term %d\n", i, j);
                         break;
                     }
                }
                else {// empty row - zeros
                     if ( const_term == -1) { // constant term for the first time
                         const_term = 0;
                         const_coeff = m_C.val(j);
                     }
                     else {
                         skip = true;
                         //printf("skipped - two constants, const %d, term %d\n", i, j);
                         break;
                     }
                }
            } // end of for on j
            
            if (skip)
            {    continue; }
            
            
            if (linear_var >= 0 && const_term == -1) // no const term
            {
                const_term = 0;
                const_coeff = 0;
            }
            
            if (const_term >= 0 && linear_var >= 0 && linear_coeff != 0)
            {// If we are here the constraint is linear with single constant and single linear term.
            
                double val = -const_coeff/linear_coeff;
//                printf("val = %f ",val);
                if (i > m_m_ineq_constrs)   
                { // Equality constraint
//                   printf("Replaced equality, var %d\n",linear_var);                    
                   x_l[linear_var] = val;
                   x_u[linear_var] = val;
                   
                }
                else 
                {// Inequality constraint
//                    printf("Replaced inequality, var %d\n",linear_var);
                    if (linear_coeff < 0 ) {
                        x_l[linear_var] = val;
                    }
                    else {
                        x_u[linear_var] = val;
                    }
                   
                    
                }
            } // end of if (const_term >= 0 && linear_var >= 0 && linear_coeff != 0)
        } // end of   "if  ((m_C.row_ptr(i+1) - m_C.row_ptr(i))<=2)"
            
    }    // end of for loop on i
  
    for(int i= 0 ; i <  m_n_fixed ;  i ++)
    {
        int idx = m_FixedStruct.col_ind(i);
        x_l[idx] = fixed[(int)m_FixedStruct.val(i)-1];
        x_u[idx] = fixed[(int)m_FixedStruct.val(i)-1];
     //   printf("%d %d\n",idx,(int)m_FixedStruct.val(i)-1);
    }
    
    for(int i=0; i < m_m_ineq_constrs; i ++)
    {
        g_l[i] = -1.0e19;
        g_u[i] = 0;
    }
    for(int i= m_m_ineq_constrs; i < m_m; i ++)
    {
        g_l[i] = 0;
        g_u[i] = 0;
    }
    
//    printf("\n end func \n");

}

