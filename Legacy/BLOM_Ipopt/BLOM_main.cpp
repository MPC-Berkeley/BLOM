// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: cpp_example.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "BLOM_NLP.hpp"

#include <iostream>
#include <time.h>
#include <sys/stat.h>

using namespace Ipopt;

int main(int argv, char* argc[])
{
  // Create an instance of your nlp...
  SmartPtr<TNLP> mynlp = new MyNLP();

  // Create an instance of the IpoptApplication
  //
  // We are using the factory, since this allows us to compile this
  // example with an Ipopt Windows DLL
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  // Initialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
//  app->Options()->SetStringValue("derivative_test", "second-order");
//  app->Options()->SetStringValue("print_user_options", "yes");
//  app->Options()->SetStringValue("print_timing_statistics", "yes");
//  app->Options()->SetStringValue("print_options_documentation", "yes");
//  app->Options()->SetStringValue("replace_bounds", "yes");
//  app->Options()->SetStringValue("inexact_algorithm", "yes");

  struct stat st;
  //printf("stat(./results/) = %d \n", stat("./results/",&st));
  stat("./results",&st);
  //printf("./results st.st_mode = %d \n", st.st_mode);
  //stat("./A.txt",&st);
  //printf("./A.txt st.st_mode = %d \n", st.st_mode);
  if (S_ISDIR(st.st_mode)) // if results folder exists, put output there
  {
    // get current time and set output_file option accordingly
    time_t timenow = time(NULL);
    struct tm * timestruct = localtime(&timenow);
    char timechar[35];
    strftime(timechar, 35, "./results/output_%y%m%d_%H%M%S.txt", timestruct);
    app->Options()->SetStringValue("output_file", timechar);
  }
  
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    return (int) status;
  }


//for (int i =1; i < 1000; i ++)
  status = app->OptimizeTNLP(mynlp);

  if (status == Solve_Succeeded) {
    // Retrieve some statistics about the solve
    Index iter_count = app->Statistics()->IterationCount();
    std::cout << std::endl << std::endl << "*** The problem solved in " << iter_count << " iterations!" << std::endl;

    Number final_obj = app->Statistics()->FinalObjective();
    std::cout << std::endl << std::endl << "*** The final value of the objective function is " << final_obj << '.' << std::endl;
  }

  return (int) status;
}
