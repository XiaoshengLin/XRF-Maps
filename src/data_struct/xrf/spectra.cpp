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



#include "spectra.h"
#include <algorithm>
#include <math.h>
#include <iostream>
#include <vector>

namespace data_struct
{
namespace xrf
{

template<typename T>
std::vector<T> conv_valid(std::vector<T> const &f, std::vector<T> const &g)
{
	size_t const nf = f.size();
	size_t const ng = g.size();
	std::vector<T> const &min_v = (nf < ng)? f : g;
	std::vector<T> const &max_v = (nf < ng)? g : f;
	size_t const n  = std::max(nf, ng) - std::min(nf, ng) + 1;
	std::vector<T> out(n, T());
	for(auto i(0); i < n; ++i)
	{
        for(int j(min_v.size() - 1), k(i); j >= 0; --j)
		{
			out[i] += min_v[j] * max_v[k];
			++k;
		}
	}
	return out;
}

EArrayXr convolve1d( EArrayXr &arr, size_t boxcar_size)
{
    EArrayXr boxcar;
    boxcar.resize(boxcar_size);
    boxcar.setOnes(boxcar_size);
    return convolve1d(arr, boxcar);
}

EArrayXr convolve1d( EArrayXr &arr, EArrayXr &boxcar)
{
	EArrayXr new_background;
	new_background.resize(arr.size());
	new_background.setZero(arr.size());
    //convolve 1d

    size_t const nf = arr.size();
	size_t const ng = boxcar.size();
    EArrayXr const &min_v = (nf < ng)? arr : boxcar;
    EArrayXr const &max_v = (nf < ng)? boxcar : arr;
	size_t const n  = std::max(nf, ng) - std::min(nf, ng) + 1;
    EArrayXr out;
    out.resize(n);
	out.setZero(n);
    for(auto i(0); i < n; ++i)
    {
        for(int j(min_v.size() - 1), k(i); j >= 0; --j)
        {
            out(i) += min_v(j) * max_v(k);
            ++k;
        }
    }
    for(size_t i=0; i< n; i++)
    {
        if( out(i) != (real_t)0.0)
        {
            new_background(i) = out(i) / real_t(boxcar.size());
        }
    }

    return new_background;
}

EArrayXr snip_background(const Spectra* const spectra,
									  real_t energy_offset,
									  real_t energy_linear,
									  real_t energy_quadratic,
								      real_t spectral_binning,
									  real_t width,
									  real_t xmin,
									  real_t xmax)
{
	EArrayXr energy;
	EArrayXr index;
	EArrayXr background;
	size_t buffer_size = spectra->size();
	energy.resize(buffer_size);
	index.resize(buffer_size);
	background.resize(buffer_size);
	for (size_t i = 0; i < buffer_size; i++)
	{
		background(i) = (*spectra)(i);
		energy(i) = real_t(i);
		index(i) = real_t(i);
	}

	if (spectral_binning > 0)
	{
		energy = energy * spectral_binning;
	}

	energy = energy_offset + energy * energy_linear + energy.pow((real_t)2.0) * energy_quadratic;

	EArrayXr tmp = std::pow((energy_offset / (real_t)2.3548), (real_t)2.0) + energy * (real_t)2.96 * energy_linear;
	for (size_t i = 0; i<tmp.size(); i++)
	{
		if (tmp(i) < 0.0)
		{
			tmp(i) = 0.0;
		}
	}
	//std::valarray<real_t> fwhm = 2.35 * std::sqrt(tmp);
	EArrayXr current_width = (real_t)2.35 * tmp.sqrt();


	EArrayXr boxcar;
	EArrayXr new_background;
	new_background.resize(background.size());
	// smooth the background
	if (spectral_binning > 0)
	{
		boxcar.resize(3);
	}
	else
	{
		boxcar.resize(5);
	}

	for (size_t i = 0; i<boxcar.size(); i++)
	{
		boxcar(i) = 1.0;
	}
	//convolve 1d
	for (size_t i = 0; i< background.size(); i++)
	{
		new_background(i) = 0.0;
		for (size_t j = 0; j<boxcar.size(); j++)
		{
			if ((i - j) >= 0)
				new_background(i) += background(i - j) * boxcar(j);
		}
	}
	background = new_background / real_t(boxcar.size());
	//clear out
	new_background.resize(1);
	boxcar.resize(1);
	//fwhm
	current_width = width * current_width / energy_linear;  // in channels
	if (spectral_binning > 0)
	{
		current_width = current_width / (real_t)2.0;
	}

	background += 1.0;
	background = background.log();
	background += 1.0;
	background = background.log();

	// FIRST SNIPPING
	int no_iterations = 2;
	if (spectral_binning > 0)
	{
		no_iterations = 3;
	}

	real_t max_of_xmin = (std::max)(xmin, (real_t)0.0);
	real_t min_of_xmax = (std::min)(xmax, real_t(buffer_size - 1));
	for (int j = 0; j<no_iterations; j++)
	{
		for (size_t k = 0; k<background.size(); k++)
		{
			real_t lo_index = k - current_width(k);
			real_t hi_index = k + current_width(k);
			if (lo_index < max_of_xmin)
			{
				lo_index = max_of_xmin;
			}
			if (hi_index > min_of_xmax)
			{
				hi_index = min_of_xmax;
			}
			real_t temp = (background(lo_index) + background(hi_index)) / 2.0;
			if (background(k) > temp)
			{
				background(k) = temp;
			}
		}
	}

	while (current_width.maxCoeff() >= 0.5)
	{
		for (size_t k = 0; k<background.size(); k++)
		{
			real_t lo_index = k - current_width(k);
			real_t hi_index = k + current_width(k);
			if (lo_index < max_of_xmin)
			{
				lo_index = max_of_xmin;
			}
			if (hi_index > min_of_xmax)
			{
				hi_index = min_of_xmax;
			}
			real_t temp = (background(lo_index) + background(hi_index)) / 2.0;
			if (background(k) > temp)
			{
				background(k) = temp;
			}
		}

		current_width = current_width / real_t(M_SQRT2); // window_rf
	}

	background -= (real_t) 1.0;
	background = background.exp();
	background -= (real_t) 1.0;
	background = background.exp();
	//background = std::exp(std::exp(background) - (real_t)1.0) - (real_t)1.0;

	for (size_t i = 0; i<background.size(); i++)
	{
		if (std::isnan(background(i)))
		{
			background(i) = 0.0;
		}
	}

	return background;

}

} //namespace data_struct
} //namespace xrf
