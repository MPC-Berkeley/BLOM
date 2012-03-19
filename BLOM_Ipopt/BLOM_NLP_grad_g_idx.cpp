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
bool MyNLP::eval_jac_g_idx(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  //if (values == NULL) {
    // return the structure of the jacobian of the constraints

#include "testeval_jac_g_idx.cpp"
        return true;
}
