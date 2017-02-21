// File: nnls.cc 
// Author: Suvrit Sra
// Time-stamp: <08 March 2011 01:40:04 PM CET --  suvrit>
// nnls solver

// Copyright (C) 2009, 2010 Suvrit Sra (suvrit@tuebingen.mpg.de)
// Copyright Max-Planck-Institute for Biological Cybernetics

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


#include "nnls.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <algorithm>

using namespace nsNNLS;

