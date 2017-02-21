// File: sparseMatrix.cc -- implements sparse matrix functionality 

// Author: Suvrit Sra <suvrit@tuebingen.mpg.de>
// (c) Copyright 2010   Suvrit Sra
// Max-Planck-Institute for Biological Cybernetics

// sparseMatrix.cc - implements sparse matrix functionality
// Copyright (C) 2010 Suvrit Sra (suvrit@tuebingen.mpg.de)

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#include "sparseMatrix.h"
#include <cstdio>
#include <string>

using namespace nsNNLS;

/// Get the (i,j) entry of the matrix
double sparseMatrix::operator()   (size_t i, size_t j)
{
  // size_t sz = cols[j+1]-cols[j];
  // int idx = binary_search(m_rowindx + m_colptrs[j], i, sz);
  // if (idx != -1) 
  //   return (m_val + m_colptrs[j])[idx];
  return 0.0;
}

/// Get the (i,j) entry of the matrix
double sparseMatrix::get (size_t i, size_t j)
{
  // size_t sz = cols[j+1]-cols[j];
  // int idx = binary_search(m_rowindx + m_colptrs[j], i, sz);
  // if (idx != -1) 
  //   return (m_val + m_colptrs[j])[idx];
  return 0.0;
}

/// Set the (i,j) entry of the matrix. If entry does not exist, function bombs.
int sparseMatrix::set (size_t i, size_t j, double val)
{
  return -1;
}
    
/// Returns 'r'-th row into pre-alloced vector
int sparseMatrix::get_row (size_t i, vector*& r)
{
  return -1;
}

/// Returns 'c'-th col as a vector
int sparseMatrix::get_col (size_t j, vector*& c)
{
  return -1;
}

/// Returns main or second diagonal (if p == true)
int sparseMatrix::get_diag(bool p, vector*& d) 
{
  return -1;
}

/// Vector l_p norms for this matrix, p > 0
double sparseMatrix::norm (double p)
{
  return -1;
}

/// Vector l_p norms, p is 'l1', 'l2', 'fro', 'inf'
double sparseMatrix::norm (const char*  p)
{
  return -1;
}

/// r = a*row(i) + r
int  sparseMatrix::row_daxpy(size_t i, double a, vector* r)
{
  return -1;
}

/// c = a*col(j) + c
int  sparseMatrix::col_daxpy(size_t j, double a, vector* c)
{
  return -1;
}
