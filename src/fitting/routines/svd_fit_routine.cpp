/***
Copyright (c) 2016, UChicago Argonne, LLC. All rights reserved.

Copyright 2016. UChicago Argonne, LLC. This software was produced
under U.S. Government contract DE-AC02-06CH11357 for Argonne National
Laboratory (ANL), which is operated by UChicago Argonne, LLC for the
U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is
modified to produce derivative works, such modified software should
be clearly marked, so as not to confuse it with the version available
from ANL.

Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.

    * Neither the name of UChicago Argonne, LLC, Argonne National
      Laboratory, ANL, the U.S. Government, nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago
Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
***/

/// Initial Author <2016>: Arthur Glowacki



#include "svd_fit_routine.h"

#include <Eigen/SVD>

//debug
#include <iostream>

namespace fitting
{
namespace routines
{

SVD_Fit_Routine::SVD_Fit_Routine() : Matrix_Optimized_Fit_Routine()
{

}

// ----------------------------------------------------------------------------

SVD_Fit_Routine::~SVD_Fit_Routine()
{

}



// ----------------------------------------------------------------------------

void SVD_Fit_Routine::_generate_fitmatrix(const unordered_map<string, Spectra> * const element_models,
                                          const struct Range energy_range)
{

    _element_row_index.clear();

    _fitmatrix.resize(energy_range.count(), element_models->size());
    int i = 0;
    for(const auto& itr : *element_models)
    {
        //Spectra element_model = itr.second;
        for (int j=0; j<itr.second.size(); j++)
        {
            _fitmatrix(j,i) = itr.second[j];
        }
        //save element index for later
        _element_row_index[itr.first] = i;
        //(*fit_params)[itr.first].opt_array_index = i;
        i++;
    }

}

// ----------------------------------------------------------------------------

void SVD_Fit_Routine::fit_spectra(const models::Base_Model * const model,
                                  const Spectra * const spectra,
                                  const Detector * const detector,
                                  const Fit_Element_Map_Dict * const elements_to_fit,
                                  Fit_Count_Dict *out_counts_dic,
                                  size_t row_idx,
                                  size_t col_idx)
{

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(_fitmatrix, Eigen::ComputeThinU | Eigen::ComputeThinV );
    Eigen::VectorXd rhs;
    rhs.resize(spectra->size());
    //assert _fitmatrix.rows() == spectra->size()

    for(size_t i=0; i<spectra->size(); i++)
    {
        rhs[i] = (*spectra)[i];
    }

    Eigen::VectorXd result = svd.solve(rhs);
    //std::cout << "SVD Result : "<<std::endl<< result <<std::endl;

    for(const auto& itr : *elements_to_fit)
    {
        //[itr.first];
        (*out_counts_dic)[itr.first][row_idx][col_idx] = result[_element_row_index[itr.first]];
    }

}

// ----------------------------------------------------------------------------

void SVD_Fit_Routine::initialize(models::Base_Model * const model,
                                 const Detector * const detector,
                                 const Fit_Element_Map_Dict * const elements_to_fit,
                                 const struct Range energy_range)
{

    unordered_map<string, Spectra> element_models = _generate_element_models(model, detector, elements_to_fit, energy_range);

    _generate_fitmatrix(&element_models, energy_range);

}

// ----------------------------------------------------------------------------

} //namespace routines
} //namespace fitting