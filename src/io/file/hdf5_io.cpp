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


#include "hdf5_io.h"

#include <iostream>
#include <string>

#include <chrono>
#include <ctime>
#include <thread>
#include <mutex>

#include "data_struct/xrf/element_info.h"

#define HDF5_SAVE_VERSION 10.0

const std::vector<std::string> hdf5_copy_dset_names = {"Element_Weights",
                                                       "Element_Weights_Names",
                                                       "DS_IC",
                                                       "US_IC",
                                                       "Standard_Name",
                                                       "Channel_Names",
                                                       "Names",
                                                       "vresion"
                                                      };


const std::vector<std::string> hdf5_copy_grp_names = {"Scan"
                                                      };

const std::vector<std::string> xrf_analysis_save_names = {"ROI",
                                                          "Params",
                                                          "Fitted",
                                                          "SVD",
                                                          "NNLS"
                                                         };





namespace io
{
namespace file
{

std::mutex HDF5_IO::_mutex;

HDF5_IO* HDF5_IO::_this_inst(0);


struct Detector_HDF5_Struct
{
    hid_t    dset_id;
    hid_t    dataspace_id;
    real_t * buffer;
};


/*
int Find_Index(hid_t file_id, loc)
{

}
*/

H5_Spectra_Layout Generate_Layout(H5_SPECTRA_LAYOUTS layout_def, size_t detector_num)
{
    H5_Spectra_Layout layout;

    layout.detector_num = detector_num;
    if (layout_def == MAPS_RAW)
    {
        switch(detector_num)
        {
        case 0:
            layout.spectra.path = "/MAPS_RAW/data_a";
            break;
        case 1:
            layout.spectra.path = "/MAPS_RAW/data_b";
            break;
        case 2:
            layout.spectra.path = "/MAPS_RAW/data_c";
            break;
        case 3:
            layout.spectra.path = "/MAPS_RAW/data_d";
            break;
        default:
            layout.spectra.path = "";
            break;
        }
        layout.spectra.set_order(SAMPLE, ROW, COL);

        layout.elt.path = "/MAPS_RAW/livetime";
        layout.elt.set_order(DETECTOR, ROW, COL);
        layout.ert.path = "/MAPS_RAW/realtime";
        layout.ert.set_order(DETECTOR, ROW, COL);
        layout.incnt.path = "/MAPS_RAW/inputcounts";
        layout.incnt.set_order(DETECTOR, ROW, COL);
        layout.outcnt.path = "/MAPS_RAW/ouputcounts";
        layout.outcnt.set_order(DETECTOR, ROW, COL);

    }
    else if (layout_def == MAPS_V9)
    {
        layout.spectra.path = "/MAPS/mca_arr";
        layout.spectra.set_order(SAMPLE, ROW, COL);
        layout.elt.path = "/MAPS/scalers"; //index from attributes /MAPS/scalers_names "ELT1"
        layout.elt.set_order(ROW, COL);
        layout.ert.path = "/MAPS/scalers"; //index from attributes /MAPS/scalers_names "ERT1"
        layout.ert.set_order(ROW, COL);
        layout.incnt.path = "/MAPS/scalers"; //index from attributes /MAPS/scalers_names "ICR1"
        layout.incnt.set_order(ROW, COL);
        layout.outcnt.path = "/MAPS/scalers"; //index from attributes /MAPS/scalers_names "OCR1"
        layout.outcnt.set_order(ROW, COL);
    }
    else if (layout_def == MAPS_V10)
    {
        layout.spectra.path = "/MAPS/Spectra/mca_arr";
        layout.spectra.set_order(SAMPLE, ROW, COL);
        layout.elt.path = "/MAPS/Spectra/Elapsed_Livetime";
        layout.elt.set_order(ROW, COL);
        layout.ert.path = "/MAPS/Spectra/Elapsed_Realtime";
        layout.ert.set_order(ROW, COL);
        layout.incnt.path = "/MAPS/Spectra/Input_Counts";
        layout.incnt.set_order(ROW, COL);
        layout.outcnt.path = "/MAPS/Spectra/Output_Counts";
        layout.outcnt.set_order(ROW, COL);
    }
    else if (layout_def == XSPRESS)
    {
        std::string live_time_dataset_name = "CHAN" + std::to_string(detector_num+1) + "SCA0";

        layout.spectra.path = "/entry/data/data";
        layout.spectra.set_order(COL, DETECTOR, SAMPLE);
        layout.elt.path = "/entry/instrument/NDAttributes/"+live_time_dataset_name;
        layout.elt.set_order(COL);
        //layout.ert.path = "";
        //layout.incnt.path = "";
        //layout.outcnt.path = "";
    }
    else if (layout_def == APS_SEC20)
    {
        layout.spectra.path = "/2D Scan/MCA "+ std::to_string(detector_num+1);
        layout.spectra.set_order(ROW, COL, SAMPLE);
        layout.elt.path = "/2D Scan/Detectors"; //index from attributes "Timer"
        layout.spectra.set_order(ROW, COL);
        //layout.ert.path = "";
        layout.incnt.path = "/2D Scan/Detectors"; //index from attributes  "ICR Ch 1"
        layout.incnt.set_order(ROW, COL);
        layout.outcnt.path = "/2D Scan/Detectors"; //index from attributes "OCR Ch 1"
        layout.outcnt.set_order(ROW, COL);
    }
    else
    {
        //throw exception
    }

    return layout;
}


//-----------------------------------------------------------------------------

HDF5_IO::HDF5_IO() : Base_File_IO()
{
	//disable hdf print to std err
	hid_t status;
	status = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    _cur_file_id = -1;
}

//-----------------------------------------------------------------------------

HDF5_IO* HDF5_IO::inst()
{
    std::lock_guard<std::mutex> lock(_mutex);

    if (_this_inst == nullptr)
    {
        _this_inst = new HDF5_IO();
    }
    return _this_inst;
}

//-----------------------------------------------------------------------------

void HDF5_IO::lazy_load()
{
   std::lock_guard<std::mutex> lock(_mutex);

   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   //logit<<"lazy_load "<< _filename <<std::endl;

   hid_t    file_id = -1;
   hid_t    dset_id = -1;
   hid_t    dataspace_id = -1;
   hid_t    memoryspace = -1;
   hid_t    datatype = -1;
   herr_t   error = -1;

    H5T_class_t dtype_class;
   //H5T_order_t order;
   //size_t      size;

    //logit<<"pre open file"<< std::endl;
    file_id = H5Fopen(_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    //logit<<"pre open dset "<<std::endl;
    dset_id = H5Dopen2(file_id, "/MAPS_RAW/data_a", H5P_DEFAULT);


    datatype = H5Dget_type(dset_id);
    dtype_class = H5Tget_class(datatype);
    //if (dtype_class == H5T_INTEGER)
    //   printf("Data set has INTEGER type \n");
    //order = H5Tget_order(datatype);
    //if (order == H5T_ORDER_LE)
    //   printf("Little endian order \n");
    dataspace_id = H5Dget_space(dset_id);
    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    hsize_t* dims_out = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    hsize_t* sel_dims = new hsize_t[rank];
    //logit<<"rank = "<<rank<< std::endl;
    unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_out[0], NULL);

    unsigned long total = 1;
    for (int i=0; i < rank; i++)
    {
       logit<<"dims ["<<i<<"] ="<<dims_out[i]<< std::endl;
       total *= dims_out[i];
       offset[i] = 0;
       sel_dims[i] = count[i] = dims_out[i];
    }
    //logit<<"total = "<<total<<std::endl;

    //unsigned short *buffer = new unsigned short[total];
    float *buffer = new float[total];
    //logit<<"allocated "<<std::endl;
    memoryspace = H5Screate_simple(rank, sel_dims, NULL);

    H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    //error = H5Dread(dset_id, H5T_NATIVE_USHORT, memoryspace, dataspace_id, H5P_DEFAULT, buffer);
    error = H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace, dataspace_id, H5P_DEFAULT, buffer);
    printf("read in: error = %d \n", error);
/*
    for(unsigned long i=0; i<total; i++)
    {
       logit<<buffer[i]<<" ";
       if( i%2048 == 1)
          logit<<std::endl;
    }
*/
    delete [] dims_out;
    delete [] offset;
    delete [] count;
    delete [] sel_dims;
    delete [] buffer;

    H5Dclose(dset_id);
    H5Sclose(memoryspace);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logit << "finished computation at " << std::ctime(&end_time)
                 << "elapsed time: " << elapsed_seconds.count() << "s\n";

    _is_loaded = LAZY_LOAD;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_dataset(std::string path, Base_Dataset *dset)
{
    return false;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::_open_h5_object(hid_t &id, H5_OBJECTS obj, std::stack<std::pair<hid_t, H5_OBJECTS> > &close_map, std::string s1, hid_t id2)
{
    if (obj == H5O_FILE)
    {
        id = H5Fopen(s1.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if(id < 0)
        {
            _close_h5_objects(close_map);
            logit<<"Error opening file "<<s1<<std::endl;
            return false;
        }
    }
    else if (obj == H5O_GROUP)
    {
        id = H5Gopen(id2, s1.c_str(), H5P_DEFAULT);
        if(id < 0)
        {
            _close_h5_objects(close_map);
            logit<<"Error opening group "<<s1<<std::endl;
            return false;
        }
    }
    else if (obj ==  H5O_DATASET)
    {
        id = H5Dopen2(id2, s1.c_str(), H5P_DEFAULT);
        if(id < 0)
        {
            _close_h5_objects(close_map);
            logit<<"Error opening dataset "<<s1<<std::endl;
            return false;
        }
    }
    else
    {
        _close_h5_objects(close_map);
        logit<<"Error unknown H5_OBJECT type "<<obj<<std::endl;
        return false;
    }

    close_map.push({id, obj});
    return true;
}

void HDF5_IO::_close_h5_objects(std::stack<std::pair<hid_t, H5_OBJECTS> >  &close_map)
{

    herr_t err;
    while (false == close_map.empty())
    {
        auto& obj = close_map.top();
        close_map.pop();

        if (obj.second == H5O_FILE)
        {
            err = H5Fclose(obj.first);
            if (err > 0)
            {
                logit<<"Warning: could not close h5 file id "<<obj.first<<std::endl;
            }
        }
        else if (obj.second == H5O_GROUP)
        {
            err = H5Gclose(obj.first);
            if (err > 0)
            {
                logit<<"Warning: could not close h5 group id "<<obj.first<<std::endl;
            }
        }
        else if (obj.second == H5O_DATASPACE)
        {
            err = H5Sclose(obj.first);
            if (err > 0)
            {
                logit<<"Warning: could not close h5 dataspace id "<<obj.first<<std::endl;
            }
        }
        else if (obj.second == H5O_DATASET)
        {
            err = H5Dclose(obj.first);
            if (err > 0)
            {
                logit<<"Warning: could not close h5 dataset id "<<obj.first<<std::endl;
            }
        }
        else if(obj.second == H5O_ATTRIBUTE)
        {
            err = H5Aclose(obj.first);
            if (err > 0)
            {
                logit<<"Warning: could not close h5 dataset id "<<obj.first<<std::endl;
            }
        }
    }
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_volume(std::string path, size_t detector_num, data_struct::xrf::Spectra_Volume* spec_vol)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

   logit<< path <<" detector : "<<detector_num<<std::endl;

   hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
   hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
   herr_t   error;
   std::string detector_path;
   real_t * buffer;
   hsize_t offset_row[2] = {0,0};
   hsize_t count_row[2] = {0,0};
   hsize_t offset_meta[3] = {0,0,0};
   hsize_t count_meta[3] = {1,1,1};


   switch(detector_num)
   {
   case 0:
       detector_path = "data_a";
       break;
   case 1:
       detector_path = "data_b";
       break;
   case 2:
       detector_path = "data_c";
       break;
   case 3:
       detector_path = "data_d";
       break;
   default:
       detector_path = "";
       break;
   }

    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
        return false;

    if ( false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "MAPS_RAW", file_id) )
       return false;

    if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id) )
       return false;
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "livetime", maps_grp_id) )
       return false;
    dataspace_lt_id = H5Dget_space(dset_lt_id);
    close_map.push({dataspace_lt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "realtime", maps_grp_id) )
       return false;
    dataspace_rt_id = H5Dget_space(dset_rt_id);
    close_map.push({dataspace_rt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "inputcounts", maps_grp_id) )
       return false;
    dataspace_inct_id = H5Dget_space(dset_incnt_id);
    close_map.push({dataspace_inct_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "ouputcounts", maps_grp_id) )
       return false;
    dataspace_outct_id = H5Dget_space(dset_outcnt_id);
    close_map.push({dataspace_outct_id, H5O_DATASPACE});

    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logit<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<std::endl;
        return false;
       //throw exception ("Dataset is not a volume");
    }
    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    logit<<"rank = "<<rank<< std::endl;
    unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logit<<"Error getting dataset rank for MAPS_RAW/"<< detector_path<<std::endl;
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       logit<<"dims ["<<i<<"] ="<<dims_in[i]<< std::endl;
       offset[i] = 0;
       count[i] = dims_in[i];
    }

    buffer = new real_t [dims_in[0] * dims_in[2]]; // spectra_size x cols
    count_row[0] = dims_in[0];
    count_row[1] = dims_in[2];

/* TODO: maybe use greatest size (like xpress) becaue x_axis and y_axis will be diff than images size
    size_t greater_rows = std::max(spec_vol->rows() , dims_in[1]);
    size_t greater_cols = std::max(spec_vol->cols() , dims_in[2]);
    size_t greater_channels = std::max(spec_vol->samples_size() , dims_in[0]);
*/

    if(spec_vol->rows() < dims_in[1] || spec_vol->cols() < dims_in[2] || spec_vol->samples_size() < dims_in[0])
    {
        spec_vol->resize(dims_in[1], dims_in[2], dims_in[0]);
    }

    count[1] = 1; //1 row

    memoryspace_id = H5Screate_simple(2, count_row, NULL);
    close_map.push({memoryspace_id, H5O_DATASPACE});
    memoryspace_meta_id = H5Screate_simple(1, count_meta, NULL);
    close_map.push({memoryspace_meta_id, H5O_DATASPACE});
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, NULL, count_row, NULL);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;

    offset_meta[0] = detector_num;
    for (size_t row=0; row < dims_in[1]; row++)
    {
         offset[1] = row;
         offset_meta[1] = row;

         H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
         error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

         if (error > -1 )
         {
             for(size_t col=0; col<count_row[1]; col++)
             {
                 offset_meta[2] = col;
                 data_struct::xrf::Spectra *spectra = &((*spec_vol)[row][col]);

                 H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                 spectra->elapsed_lifetime(live_time);

                 H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                 spectra->elapsed_realtime(real_time);

                 H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                 spectra->input_counts(in_cnt);

                 H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                 spectra->output_counts(out_cnt);

                 spectra->recalc_elapsed_lifetime();

                 for(size_t s=0; s<count_row[0]; s++)
                 {
                     (*spectra)[s] = buffer[(count_row[1] * s) + col];
                 }
                 //logit<<"saved col "<<col<<std::endl;
             }

            //logit<<"read row "<<row<<std::endl;
         }
         else
         {
            logit<<"Error: reading row "<<row<<std::endl;
         }
    }

    delete [] dims_in;
    delete [] offset;
    delete [] count;
    delete [] buffer;

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;
    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_line_xspress3(std::string path, size_t detector_num, data_struct::xrf::Spectra_Line* spec_row)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   logit<< path <<" detector : "<<detector_num<<std::endl;

   std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
   hid_t    file_id, dset_id, dataspace_id, maps_grp_id, scaler_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
   hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
   herr_t   error;
   std::string detector_path;
   real_t * buffer;
   hsize_t offset_row[3] = {0,0,0};
   hsize_t count_row[3] = {0,0,0};
   hsize_t offset_meta[1] = {0};
   hsize_t count_meta[1] = {1};

   std::string live_time_dataset_name = "CHAN" + std::to_string(detector_num+1) + "SCA0";
   //std::string real_time_dataset_name = "CHAN" + std::to_string(detector_num+1) + "EventWidth";


    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
       return false;

    if ( false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "/entry/data", file_id) )
      return false;

    if ( false == _open_h5_object(scaler_grp_id, H5O_GROUP, close_map, "/entry/instrument/NDAttributes", file_id) )
      return false;

    if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, "data", maps_grp_id) )
      return false;
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, live_time_dataset_name, scaler_grp_id) )
      return false;
    dataspace_lt_id = H5Dget_space(dset_lt_id);
    close_map.push({dataspace_lt_id, H5O_DATASPACE});

    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logit<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<std::endl;
        return false;
       //throw exception ("Dataset is not a volume");
    }
    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    logit<<"rank = "<<rank<< std::endl;
    unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logit<<"Error getting dataset rank for MAPS_RAW/"<< detector_path<<std::endl;
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       logit<<"dims ["<<i<<"] ="<<dims_in[i]<< std::endl;
       offset[i] = 0;
       count[i] = dims_in[i];
    }

    buffer = new real_t [dims_in[0] * dims_in[2]]; // cols x spectra_size
    count_row[0] = dims_in[0];
    count_row[1] = 1;
    count_row[2] = dims_in[2];

    size_t greater_cols = std::max(spec_row->size() , (size_t)dims_in[0]);
    size_t greater_channels = std::max(spec_row[0].size() , (size_t)dims_in[2]);

    if( spec_row->size() < dims_in[0] || spec_row[0].size() < dims_in[2])
    {
        spec_row->resize(greater_cols, greater_channels);
    }

    memoryspace_id = H5Screate_simple(3, count_row, NULL);
    close_map.push({memoryspace_id, H5O_DATASPACE});
    memoryspace_meta_id = H5Screate_simple(1, count_meta, NULL);
    close_map.push({memoryspace_meta_id, H5O_DATASPACE});
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, NULL, count_row, NULL);
    //H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;

    //offset[1] = row;

    offset_row[1] = detector_num;

    H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset_row, NULL, count_row, NULL);
    error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);


    if (error > -1 )
    {

        for(size_t col=0; col < dims_in[0]; col++)
        {
            offset_meta[0] = col;
            data_struct::xrf::Spectra *spectra = &((*spec_row)[col]);

            H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
            error = H5Dread(dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
            spectra->elapsed_lifetime(live_time * 0.000000125 );
            /*
            H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
            error = H5Dread(dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
            spectra->elapsed_realtime(real_time);

            H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
            error = H5Dread(dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
            spectra->input_counts(in_cnt);

            H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
            error = H5Dread(dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
            spectra->output_counts(out_cnt);

            spectra->recalc_elapsed_lifetime();
            */


            for(size_t s=0; s<count_row[2]; s++) // num samples
            {
                //(*spectra)[s] = buffer[(count_row[1] * s) + col];
                (*spectra)[s] = buffer[(count_row[2] * col) + s];
            }
            //logit<<"saved col "<<col<<std::endl;
        }
    }

    delete [] dims_in;
    delete [] offset;
    delete [] count;
    delete [] buffer;

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;
    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_volume_confocal(std::string path, size_t detector_num, data_struct::xrf::Spectra_Volume* spec_vol)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

   logit<< path <<" detector : "<<detector_num<<std::endl;

   hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_xpos_id, dset_ypos_id, dset_detectors_id;
   hid_t    dataspace_detectors_id, dataspace_xpos_id, dataspace_ypos_id;
   hid_t    attr_detector_names_id, attr_timebase_id;
   herr_t   error;
   std::string detector_path;
   char* detector_names[256];
   real_t * buffer;
   hsize_t offset_row[2] = {0,0};
   hsize_t count_row[2] = {0,0};
   hsize_t offset_meta[3] = {0,0,0};
   hsize_t count_meta[3] = {1,1,1};
   std::unordered_map<std::string, int> detector_lookup;
    std::string incnt_str = "ICR Ch ";
    std::string outcnt_str = "OCR Ch ";

   switch(detector_num)
   {
   case 0:
       detector_path = "MCA 1";
       incnt_str += "1";
       outcnt_str += "1";
       break;
   case 1:
       detector_path = "MCA 2";
       incnt_str += "2";
       outcnt_str += "2";
       break;
   case 2:
       detector_path = "MCA 3";
       incnt_str += "3";
       outcnt_str += "3";
       break;
   case 3:
       detector_path = "MCA 4";
       incnt_str += "4";
       outcnt_str += "4";
       break;
   default:
       detector_path = "";
       break;
   }

    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
        return false;

    if ( false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "2D Scan", file_id) )
       return false;

    if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id) )
       return false;
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_detectors_id, H5O_DATASET, close_map, "Detectors", maps_grp_id) )
       return false;
    dataspace_detectors_id = H5Dget_space(dset_detectors_id);
    close_map.push({dataspace_detectors_id, H5O_DATASPACE});
    attr_detector_names_id=H5Aopen(dset_detectors_id, "Detector Names", H5P_DEFAULT);
    close_map.push({dataspace_detectors_id, H5O_ATTRIBUTE});
    attr_timebase_id=H5Aopen(dset_detectors_id, "TIMEBASE", H5P_DEFAULT);
    close_map.push({dataspace_detectors_id, H5O_ATTRIBUTE});

    if ( false == _open_h5_object(dset_xpos_id, H5O_DATASET, close_map, "X Positions", maps_grp_id) )
       return false;
    dataspace_xpos_id = H5Dget_space(dset_xpos_id);
    close_map.push({dataspace_xpos_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_ypos_id, H5O_DATASET, close_map, "Y Positions", maps_grp_id) )
       return false;
    dataspace_ypos_id = H5Dget_space(dset_ypos_id);
    close_map.push({dataspace_ypos_id, H5O_DATASPACE});

    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logit<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<std::endl;
        return false;
       //throw exception ("Dataset is not a volume");
    }
    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    logit<<"rank = "<<rank<< std::endl;
    unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logit<<"Error getting dataset rank for MAPS_RAW/"<< detector_path<<std::endl;
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       logit<<"dims ["<<i<<"] ="<<dims_in[i]<< std::endl;
       offset[i] = 0;
       count[i] = dims_in[i];
    }


    //chunking is 1 x col x samples
    buffer = new real_t [dims_in[1] * dims_in[2]]; //  cols x spectra_size
    count_row[0] = dims_in[1];
    count_row[1] = dims_in[2];

// TODO: maybe use greatest size (like xpress) becaue x_axis and y_axis will be diff than images size
//    size_t greater_rows = std::max(spec_vol->rows() , dims_in[1]);
//    size_t greater_cols = std::max(spec_vol->cols() , dims_in[2]);
//    size_t greater_channels = std::max(spec_vol->samples_size() , dims_in[0]);


    if(spec_vol->rows() < dims_in[0] || spec_vol->cols() < dims_in[1] || spec_vol->samples_size() < dims_in[2])
    {
        spec_vol->resize(dims_in[0], dims_in[1], dims_in[2]);
    }



    int det_rank = H5Sget_simple_extent_ndims(dataspace_detectors_id);
    hsize_t* det_dims_in = new hsize_t[det_rank];
    H5Sget_simple_extent_dims(dataspace_detectors_id, &det_dims_in[0], NULL);

    hid_t ftype = H5Aget_type(attr_detector_names_id);
    hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
    error = H5Aread(attr_detector_names_id, type, detector_names);
    if(error == 0)
    {
        //generate lookup map
        for(int z = 0; z < det_dims_in[2]; z++)
        {
            detector_lookup[std::string(detector_names[z])] = z;
            free(detector_names[z]);
        }
    }
    delete [] det_dims_in;


    count[0] = 1; //1 row

    memoryspace_id = H5Screate_simple(2, count_row, NULL);
    close_map.push({memoryspace_id, H5O_DATASPACE});
    memoryspace_meta_id = H5Screate_simple(1, count_meta, NULL);
    close_map.push({memoryspace_meta_id, H5O_DATASPACE});
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, NULL, count_row, NULL);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;


    for (size_t row=0; row < dims_in[0]; row++)
    {
         offset[0] = row;
         offset_meta[0] = row;

         H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
         error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

         if (error > -1 )
         {
             for(size_t col=0; col<dims_in[1]; col++)
             {
                 offset_meta[1] = col;

                 data_struct::xrf::Spectra *spectra = &((*spec_vol)[row][col]);

                 offset_meta[2] = detector_lookup["Timer"];
                 H5Sselect_hyperslab (dataspace_detectors_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_detectors_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &live_time);
                 spectra->elapsed_lifetime(live_time);

//                 H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
//                 error = H5Dread(dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
//                 spectra->elapsed_realtime(real_time);

                 offset_meta[2] = detector_lookup[incnt_str];
                 H5Sselect_hyperslab (dataspace_detectors_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_detectors_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &in_cnt);
                 spectra->input_counts(in_cnt);

                 offset_meta[2] = detector_lookup[outcnt_str];
                 H5Sselect_hyperslab (dataspace_detectors_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_detectors_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &out_cnt);
                 spectra->output_counts(out_cnt);

//                 spectra->recalc_elapsed_lifetime();

                 for(size_t s=0; s<dims_in[2]; s++)
                 {
                     (*spectra)[s] = buffer[(col * dims_in[2] ) + s];
                 }
             }
         }
         else
         {
            logit<<"Error: reading row "<<row<<std::endl;
         }
    }

    delete [] dims_in;
    delete [] offset;
    delete [] count;
    delete [] buffer;

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;
    return true;
}

/*
bool HDF5_IO::_load_spectra_volume(std::string path, size_t detector_num, H5_Spectra_Layout layout, data_struct::xrf::Spectra_Volume* spec_vol)
{
    std::lock_guard<std::mutex> lock(_mutex);

   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

   logit<< path <<" detector : "<<detector_num<<std::endl;

   hid_t    file_id, dset_id, dataspace_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
   hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
   herr_t   error;
   real_t * buffer;
   hsize_t offset_row[2] = {0,0};
   hsize_t count_row[2] = {0,0};
   hsize_t offset_meta[3] = {0,0,0};
   hsize_t count_meta[3] = {1,1,1};

   //layout order
   int spec_row_idx = layout.spectra.get_order_index(ROW);
   int spec_col_idx = layout.spectra.get_order_index(COL);
   int spec_sample_idx = layout.spectra.get_order_index(SAMPLE);
   int spec_detector_idx = layout.spectra.get_order_index(DETECTOR);

   if(spec_sample_idx < 0)
   {
       logit<<"Error loading, malformed layout. Could not find Spectra samples index"<<std::endl;
       return false;
   }


    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
        return false;

    if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, layout.spectra.path, file_id) )
       return false;
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, layout.elt.path, file_id) )
       return false;
    dataspace_lt_id = H5Dget_space(dset_lt_id);
    close_map.push({dataspace_lt_id, H5O_DATASPACE});

    if(layout.ert.path.length() > 0)
    {
        if ( false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, layout.ert.path, file_id) )
           return false;
        dataspace_rt_id = H5Dget_space(dset_rt_id);
        close_map.push({dataspace_rt_id, H5O_DATASPACE});
    }

    if(layout.incnt.path.length() > 0)
    {
        if ( false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, layout.incnt.path, file_id) )
           return false;
        dataspace_inct_id = H5Dget_space(dset_incnt_id);
        close_map.push({dataspace_inct_id, H5O_DATASPACE});
    }

    if(layout.outcnt.path.length() > 0)
    {
        if ( false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, layout.outcnt.path, file_id) )
           return false;
        dataspace_outct_id = H5Dget_space(dset_outcnt_id);
        close_map.push({dataspace_outct_id, H5O_DATASPACE});
    }

    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != layout.spectra.order_size)
    {
        _close_h5_objects(close_map);
        logit<<"Dataset "<<layout.spectra.path<<" rank != "<< layout.spectra.order_size<<". rank = "<<rank<<". Can't load dataset. returning"<<std::endl;
        return false;
    }
    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    logit<<"rank = "<<rank<< std::endl;
    unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logit<<"Error getting dataset rank for "<< layout.spectra.path<<std::endl;
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       logit<<"dims ["<<i<<"] ="<<dims_in[i]<< std::endl;
       offset[i] = 0;
       count[i] = dims_in[i];
    }



    if(spec_row_idx <0)
    {
        //load_single_row();
    }
    else
    {
        count[spec_row_idx] = 1;
    }

    if(spec_col_idx > -1)
    {
        count[spec_col_idx] = 1;
    }



//    buffer = new real_t [dims_in[0] * dims_in[2]]; // spectra_size x cols
//    count_row[0] = dims_in[0];
//    count_row[1] = dims_in[2];



    if(spec_vol->rows() < dims_in[spec_row_idx] || spec_vol->cols() < dims_in[spec_col_idx] || spec_vol->samples_size() < dims_in[spec_sample_idx])
    {
        spec_vol->resize(dims_in[spec_row_idx], dims_in[spec_col_idx], dims_in[spec_sample_idx]);
    }



    memoryspace_id = H5Screate_simple(2, count_row, NULL);
    close_map.push({memoryspace_id, H5O_DATASPACE});
    memoryspace_meta_id = H5Screate_simple(1, count_meta, NULL);
    close_map.push({memoryspace_meta_id, H5O_DATASPACE});
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, NULL, count_row, NULL);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);

    offset_meta[meta_detector_idx] = detector_num;
    for (size_t row=0; row < dims_in[spec_row_idx]; row++)
    {
         offset[spec_row_idx] = row;
         offset_meta[spec_row_idx] = row;

         H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
         error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

         if (error > -1 )
         {
             for(size_t col=0; col<dims_in[spec_col_idx]; col++)
             {
                 real_t live_time = 1.0;
                 real_t real_time = 1.0;
                 real_t in_cnt = 1.0;
                 real_t out_cnt = 1.0;

                 offset_meta[2] = col;
                 data_struct::xrf::Spectra *spectra = &((*spec_vol)[row][col]);

                 H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                 spectra->elapsed_lifetime(live_time);

                 H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                 spectra->elapsed_realtime(real_time);

                 H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                 spectra->input_counts(in_cnt);

                 H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                 spectra->output_counts(out_cnt);

                 spectra->recalc_elapsed_lifetime();

                 for(size_t s=0; s<count_row[0]; s++)
                 {
                     (*spectra)[s] = buffer[(count_row[1] * s) + col];
                 }
             }
         }
         else
         {
            logit<<"Error: reading row "<<row<<std::endl;
         }
    }

    delete [] dims_in;
    delete [] offset;
    delete [] count;
    delete [] buffer;

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;

}
*/
//-----------------------------------------------------------------------------


bool HDF5_IO::load_spectra_volume_with_callback(std::string path,
                                                size_t detector_num_start,
                                                size_t detector_num_end,
                                                IO_Callback_Func_Def callback_func,
                                                void* user_data)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();
   std::map<size_t, struct Detector_HDF5_Struct> detector_hid_map;

   std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

   logit<< path <<" detectors : "<<detector_num_start<<":"<<detector_num_end<<std::endl;

   hid_t    file_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
   hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
   herr_t   error;
   std::string detector_path;
   hsize_t offset_row[2] = {0,0};
   hsize_t count_row[2] = {0,0};
   hsize_t offset_meta[3] = {0,0,0};
   hsize_t count_meta[3] = {1,1,1};

   if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
       return false;

   if ( false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "MAPS_RAW", file_id) )
      return false;

    logit<<"pre open dset "<<std::endl;
    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {
        detector_hid_map.insert( {detector_num, Detector_HDF5_Struct()} );

        switch(detector_num)
        {
        case 0:
            detector_path = "data_a";
            break;
        case 1:
            detector_path = "data_b";
            break;
        case 2:
            detector_path = "data_c";
            break;
        case 3:
            detector_path = "data_d";
            break;
        default:
            detector_path = "";
            break;
        }


        if ( false == _open_h5_object(detector_hid_map[detector_num].dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id) )
           return false;
        detector_hid_map[detector_num].dataspace_id = H5Dget_space(detector_hid_map[detector_num].dset_id);
        close_map.push({detector_hid_map[detector_num].dataspace_id, H5O_DATASPACE});
    }

    if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "livetime", maps_grp_id) )
       return false;
    dataspace_lt_id = H5Dget_space(dset_lt_id);
    close_map.push({dataspace_lt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "realtime", maps_grp_id) )
       return false;
    dataspace_rt_id = H5Dget_space(dset_rt_id);
    close_map.push({dataspace_rt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "inputcounts", maps_grp_id) )
       return false;
    dataspace_inct_id = H5Dget_space(dset_incnt_id);
    close_map.push({dataspace_inct_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "ouputcounts", maps_grp_id) )
       return false;
    dataspace_outct_id = H5Dget_space(dset_outcnt_id);
    close_map.push({dataspace_outct_id, H5O_DATASPACE});


    int rank = H5Sget_simple_extent_ndims(detector_hid_map[detector_num_start].dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logit<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<std::endl;
        return false;
       //throw exception ("Dataset is not a volume");
    }
    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    logit<<"rank = "<<rank<< std::endl;
    unsigned int status_n = H5Sget_simple_extent_dims(detector_hid_map[detector_num_start].dataspace_id, &dims_in[0], NULL);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logit<<"Error getting dataset rank for MAPS_RAW/"<< detector_path<<std::endl;
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       logit<<"dims ["<<i<<"] ="<<dims_in[i]<< std::endl;
       offset[i] = 0;
       count[i] = dims_in[i];
    }

    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {
        detector_hid_map[detector_num].buffer = new real_t [dims_in[0] * dims_in[2]]; // spectra_size x cols
    }
    count_row[0] = dims_in[0];
    count_row[1] = dims_in[2];

    //if(spec_vol->rows() < dims_in[1] || spec_vol->cols() < dims_in[2] || spec_vol->samples_size() < dims_in[0])
    //{
    //    spec_vol->resize(dims_in[1], dims_in[2], dims_in[0]);
    //}

    count[1] = 1; //1 row

    memoryspace_id = H5Screate_simple(2, count_row, NULL);
    memoryspace_meta_id = H5Screate_simple(1, count_meta, NULL);
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, NULL, count_row, NULL);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;

    for (size_t row=0; row < dims_in[1]; row++)
    {
         offset[1] = row;
         offset_meta[1] = row;

         for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
         {
            H5Sselect_hyperslab (detector_hid_map[detector_num].dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            error = H5Dread(detector_hid_map[detector_num].dset_id, H5T_NATIVE_REAL, memoryspace_id, detector_hid_map[detector_num].dataspace_id, H5P_DEFAULT, detector_hid_map[detector_num].buffer);
         }

         if (error > -1 )
         {
             for(size_t col=0; col<count_row[1]; col++)
             {
                 offset_meta[2] = col;

                 for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
                 {
                     offset_meta[0] = detector_num;
                     data_struct::xrf::Spectra * spectra = new data_struct::xrf::Spectra(dims_in[0]);

                     H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                     error = H5Dread(dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                     spectra->elapsed_lifetime(live_time);

                     H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                     error = H5Dread(dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                     spectra->elapsed_realtime(real_time);

                     H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                     error = H5Dread(dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                     spectra->input_counts(in_cnt);

                     H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                     error = H5Dread(dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                     spectra->output_counts(out_cnt);

                     for(size_t s=0; s<count_row[0]; s++)
                     {
                         (*spectra)[s] = detector_hid_map[detector_num].buffer[(count_row[1] * s) + col];
                     }
                     callback_func(row, col, dims_in[1], count_row[1], detector_num, spectra, user_data);
                }
             }

         }
         else
         {
            logit<<"Error: reading row "<<row<<std::endl;
         }
    }

    delete [] dims_in;
    delete [] offset;
    delete [] count;

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_and_integrate_spectra_volume(std::string path, size_t detector_num, data_struct::xrf::Spectra* spectra)
{
    std::lock_guard<std::mutex> lock(_mutex);


    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    logit<< path <<" detector : "<<detector_num<<std::endl;

    hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
    hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
    herr_t   error;
    std::string detector_path;
    real_t * buffer;
    hsize_t offset_row[2] = {0,0};
    hsize_t count_row[2] = {0,0};
    hsize_t offset_meta[3] = {0,0,0};
    hsize_t count_meta[3] = {1,1,1};


    switch(detector_num)
    {
    case 0:
        detector_path = "data_a";
        break;
    case 1:
        detector_path = "data_b";
        break;
    case 2:
        detector_path = "data_c";
        break;
    case 3:
        detector_path = "data_d";
        break;
    default:
        detector_path = "";
        break;
    }

    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
        return false;

    if ( false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "MAPS_RAW", file_id) )
       return false;

    if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id) )
       return false;
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "livetime", maps_grp_id) )
       return false;
    dataspace_lt_id = H5Dget_space(dset_lt_id);
    close_map.push({dataspace_lt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "realtime", maps_grp_id) )
       return false;
    dataspace_rt_id = H5Dget_space(dset_rt_id);
    close_map.push({dataspace_rt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "inputcounts", maps_grp_id) )
       return false;
    dataspace_inct_id = H5Dget_space(dset_incnt_id);
    close_map.push({dataspace_inct_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "ouputcounts", maps_grp_id) )
       return false;
    dataspace_outct_id = H5Dget_space(dset_outcnt_id);
    close_map.push({dataspace_outct_id, H5O_DATASPACE});


     int rank = H5Sget_simple_extent_ndims(dataspace_id);
     if (rank != 3)
     {
         _close_h5_objects(close_map);
         logit<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<std::endl;
         return false;
        //throw exception ("Dataset is not a volume");
     }
     hsize_t* dims_in = new hsize_t[rank];
     hsize_t* offset = new hsize_t[rank];
     hsize_t* count = new hsize_t[rank];
     logit<<"rank = "<<rank<< std::endl;
     unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
     if(status_n < 0)
     {
         _close_h5_objects(close_map);
          logit<<"Error getting dataset rank for MAPS_RAW/"<< detector_path<<std::endl;
          return false;
     }

     for (int i=0; i < rank; i++)
     {
        logit<<"dims ["<<i<<"] ="<<dims_in[i]<< std::endl;
        offset[i] = 0;
        count[i] = dims_in[i];
     }

     buffer = new real_t [dims_in[0] * dims_in[2]]; // spectra_size x cols
     count_row[0] = dims_in[0];
     count_row[1] = dims_in[2];

     count[1] = 1; //1 row

     if(spectra->size() != dims_in[0])
     {
        spectra->resize(dims_in[0]);
     }

     memoryspace_id = H5Screate_simple(2, count_row, NULL);
     memoryspace_meta_id = H5Screate_simple(1, count_meta, NULL);
     H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, NULL, count_row, NULL);
     H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);

     real_t live_time = 1.0;
     real_t real_time = 1.0;
     real_t in_cnt = 1.0;
     real_t out_cnt = 1.0;

     real_t live_time_total = 0.0;
     real_t real_time_total = 0.0;
     real_t in_cnt_total = 0.0;
     real_t out_cnt_total = 0.0;

     offset_meta[0] = detector_num;
     for (size_t row=0; row < dims_in[1]; row++)
     {
          offset[1] = row;
          offset_meta[1] = row;

          H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
          error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

          if (error > -1 )
          {
              for(size_t col=0; col<count_row[1]; col++)
              {
                  offset_meta[2] = col;

                  H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                  error = H5Dread(dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                  live_time_total += live_time;

                  H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                  error = H5Dread(dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                  real_time_total += real_time;

                  H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                  error = H5Dread(dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                  in_cnt_total += in_cnt;

                  H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                  error = H5Dread(dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                  out_cnt_total += out_cnt;

                  for(size_t s=0; s<count_row[0]; s++)
                  {
                      (*spectra)[s] += buffer[(count_row[1] * s) + col];
                  }
                  //logit<<"saved col "<<col<<std::endl;
              }

             logit<<"read row "<<row<<std::endl;
          }
          else
          {
             logit<<"Error: reading row "<<row<<std::endl;
          }
     }

     spectra->elapsed_lifetime(live_time_total);
     spectra->elapsed_realtime(real_time_total);
     spectra->input_counts(in_cnt_total);
     spectra->output_counts(out_cnt_total);

     delete [] dims_in;
     delete [] offset;
     delete [] count;
     delete [] buffer;

     _close_h5_objects(close_map);

     end = std::chrono::system_clock::now();
     std::chrono::duration<double> elapsed_seconds = end-start;
     //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

     logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;

     return true;
}

//-----------------------------------------------------------------------------


bool HDF5_IO::load_spectra_vol_analyzed_h5(std::string path,
                                           size_t detector_num,
                                           data_struct::xrf::Spectra_Volume* spectra_volume,
                                           int row_idx_start,
                                           int row_idx_end,
                                           int col_idx_start,
                                           int col_idx_end)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

   logit<< path <<" detector : "<<detector_num<<std::endl;

   path += ".h5" + std::to_string(detector_num);

   hid_t    file_id, dset_id, dataspace_id, spec_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
   hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
   herr_t   error;
   hsize_t dims_in[3] = {0,0,0};
   hsize_t offset[3] = {0,0,0};
   hsize_t count[3] = {1,1,1};
   hsize_t offset_time[2] = {0,0};
   hsize_t count_time[2] = {1,1};

   if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
       return false;

   if ( false == _open_h5_object(spec_grp_id, H5O_GROUP, close_map, "/MAPS/Spectra", file_id) )
      return false;

   if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, "mca_arr", spec_grp_id) )
      return false;
   dataspace_id = H5Dget_space(dset_id);
   close_map.push({dataspace_id, H5O_DATASPACE});

   if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "Elapsed_Livetime", spec_grp_id) )
      return false;
   dataspace_lt_id = H5Dget_space(dset_lt_id);
   close_map.push({dataspace_lt_id, H5O_DATASPACE});

   if ( false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "Elapsed_Realtime", spec_grp_id) )
      return false;
   dataspace_rt_id = H5Dget_space(dset_rt_id);
   close_map.push({dataspace_rt_id, H5O_DATASPACE});

   if ( false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "Input_Counts", spec_grp_id) )
      return false;
   dataspace_inct_id = H5Dget_space(dset_incnt_id);
   close_map.push({dataspace_inct_id, H5O_DATASPACE});

   if ( false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "Output_Counts", spec_grp_id) )
      return false;
   dataspace_outct_id = H5Dget_space(dset_outcnt_id);
   close_map.push({dataspace_outct_id, H5O_DATASPACE});


    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logit<<"Dataset /MAPS/Spectra/mca_arr  rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<std::endl;
        return false;
       //throw exception ("Dataset is not a volume");
    }
    logit<<"rank = "<<rank<< std::endl;
    unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logit<<"Error getting dataset dims for /MAPS/Spectra/mca_arr"<<std::endl;
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       logit<<"dims ["<<i<<"] ="<<dims_in[i]<< std::endl;
       offset[i] = 0;
       count[i] = dims_in[i];
    }

    if(row_idx_end < row_idx_start)
    {
        row_idx_end = dims_in[1];
    }

    if(col_idx_end < col_idx_start)
    {
        col_idx_end = dims_in[2];
    }

    spectra_volume->resize(dims_in[1], dims_in[2], dims_in[0]);

    //buffer = new real_t [dims_in[0] * dims_in[2]]; // cols x spectra_size
    count[0] = dims_in[0];
    count[1] = 1;
    count[2] = 1;

    memoryspace_id = H5Screate_simple(3, count, NULL);
    memoryspace_meta_id = H5Screate_simple(2, count_time, NULL);
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;

    //offset[1] = row;

    for(size_t row=row_idx_start; row < row_idx_end; row++)
    {
        offset[1] = row;
        offset_time[0] = row;
        for(size_t col=col_idx_start; col < col_idx_end; col++)
        {
            data_struct::xrf::Spectra *spectra = &((*spectra_volume)[row][col]);
            offset[2] = col;
            offset_time[1] = col;
            H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

            error = H5Dread (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(*spectra)[0]);

            H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);
            H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);
            H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);
            H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);

            error = H5Dread (dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, (void*)&real_time);
            error = H5Dread (dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, (void*)&live_time);
            error = H5Dread (dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, (void*)&in_cnt);
            error = H5Dread (dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, (void*)&out_cnt);

            spectra->elapsed_lifetime(live_time);
            spectra->elapsed_realtime(real_time);
            spectra->input_counts(in_cnt);
            spectra->output_counts(out_cnt);
        }
    }

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_integrated_spectra_analyzed_h5(std::string path,
                                           size_t detector_num,
                                           data_struct::xrf::Spectra* spectra)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   //std::chrono::time_point<std::chrono::system_clock> start, end;
   //start = std::chrono::system_clock::now();

   logit<< path <<" detector : "<<detector_num<<std::endl;

   path += ".h5" + std::to_string(detector_num);

   hid_t    file_id;

    file_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if(file_id < 0)
    {

        logit<<"Error opening file "<<path<<std::endl;
        return false;
    }

    return load_integrated_spectra_analyzed_h5(file_id, spectra);
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_integrated_spectra_analyzed_h5(hid_t file_id, data_struct::xrf::Spectra* spectra)
{

    hid_t    dset_id, dataspace_id, spec_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
    hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
    herr_t   error;
    hsize_t dims_in[1] = {0};
    hsize_t offset[1] = {0};
    hsize_t count[1] = {1};
    hsize_t offset_time[1] = {0};
    hsize_t count_time[1] = {1};

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    close_map.push({file_id, H5O_FILE});

    if ( false == _open_h5_object(spec_grp_id, H5O_GROUP, close_map, "/MAPS/Spectra/Integrated_Spectra", file_id) )
       return false;

    if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, "Spectra", spec_grp_id) )
       return false;
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "Elapsed_Livetime", spec_grp_id) )
       return false;
    dataspace_lt_id = H5Dget_space(dset_lt_id);
    close_map.push({dataspace_lt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "Elapsed_Realtime", spec_grp_id) )
       return false;
    dataspace_rt_id = H5Dget_space(dset_rt_id);
    close_map.push({dataspace_rt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "Input_Counts", spec_grp_id) )
       return false;
    dataspace_inct_id = H5Dget_space(dset_incnt_id);
    close_map.push({dataspace_inct_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "Output_Counts", spec_grp_id) )
       return false;
    dataspace_outct_id = H5Dget_space(dset_outcnt_id);
    close_map.push({dataspace_outct_id, H5O_DATASPACE});


    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 1)
    {
        _close_h5_objects(close_map);
        logit<<"Dataset /MAPS/Spectra/mca_arr  rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<std::endl;
        return false;
       //throw exception ("Dataset is not a volume");
    }
    logit<<"rank = "<<rank<< std::endl;
    unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logit<<"Error getting dataset dims for /MAPS/Spectra/mca_arr"<<std::endl;
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       logit<<"dims ["<<i<<"] ="<<dims_in[i]<< std::endl;
       offset[i] = 0;
       count[i] = dims_in[i];
    }

    spectra->resize(dims_in[0]);

    count[0] = dims_in[0];

    memoryspace_id = H5Screate_simple(1, count, NULL);
    memoryspace_meta_id = H5Screate_simple(1, count_time, NULL);
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;

    H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    error = H5Dread (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(*spectra)[0]);

    H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);
    H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);
    H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);
    H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);

    error = H5Dread (dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, (void*)&real_time);
    error = H5Dread (dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, (void*)&live_time);
    error = H5Dread (dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, (void*)&in_cnt);
    error = H5Dread (dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, (void*)&out_cnt);

    spectra->elapsed_lifetime(live_time);
    spectra->elapsed_realtime(real_time);
    spectra->input_counts(in_cnt);
    spectra->output_counts(out_cnt);

    _close_h5_objects(close_map);

    //end = std::chrono::system_clock::now();
    //std::chrono::duration<double> elapsed_seconds = end-start;
    ////std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    //logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::start_save_seq(const std::string filename, bool force_new_file)
{

    if (_cur_file_id > -1)
    {
        logit<<"Warning: file already open, calling close() before opening new file. "<<std::endl;
        end_save_seq();
    }

    if(false == force_new_file)
        _cur_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (_cur_file_id < 1)
    {
        logit<<"Creating file "<<filename<<std::endl;
        _cur_file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    if(_cur_file_id < 0)
    {
        logit<<"Error opening file "<<filename<<std::endl;
        return false;
    }

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::end_save_seq()
{

    if(_cur_file_id > 0)
    {
        logit<<"closing file\n";
        logit<<"=========================================================\n\n";
        H5Fflush(_cur_file_id, H5F_SCOPE_LOCAL);

        ssize_t obj_cnt = H5Fget_obj_count( _cur_file_id, H5F_OBJ_DATASET | H5F_OBJ_LOCAL );
        if(obj_cnt > 0)
        {
            logit<<" closing forgotten datasets: "<<obj_cnt<<std::endl;
            hid_t* objs = new hid_t[obj_cnt];
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_DATASET, -1, objs ) > -1)
            {
                for(int i=0; i<obj_cnt; i++)
                {
                    H5Dclose(objs[i]);
                }
            }
            delete [] objs;
        }
        obj_cnt = H5Fget_obj_count( _cur_file_id, H5F_OBJ_GROUP | H5F_OBJ_LOCAL );
        if(obj_cnt > 0)
        {
            logit<<" closing forgotten groups: "<<obj_cnt<<std::endl;
            hid_t* objs = new hid_t[obj_cnt];
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_GROUP, -1, objs ) > -1)
            {
                for(int i=0; i<obj_cnt; i++)
                {
                    H5Gclose(objs[i]);
                }
            }
            delete [] objs;
        }
        obj_cnt = H5Fget_obj_count( _cur_file_id, H5F_OBJ_DATATYPE | H5F_OBJ_LOCAL );
        if(obj_cnt > 0)
        {
            logit<<" closing forgotten datatypes: "<<obj_cnt<<std::endl;
            hid_t* objs = new hid_t[obj_cnt];
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_DATATYPE, -1, objs ) > -1)
            {
                for(int i=0; i<obj_cnt; i++)
                {
                    H5Tclose(objs[i]);
                }
            }
            delete [] objs;
        }
        obj_cnt = H5Fget_obj_count( _cur_file_id, H5F_OBJ_ATTR | H5F_OBJ_LOCAL );
        if(obj_cnt > 0)
        {
            logit<<" closing forgotten attributes: "<<obj_cnt<<std::endl;
            hid_t* objs = new hid_t[obj_cnt];
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_ATTR, -1, objs ) > -1)
            {
                for(int i=0; i<obj_cnt; i++)
                {
                    H5Aclose(objs[i]);
                }
            }
            delete [] objs;
        }
        obj_cnt = H5Fget_obj_count( _cur_file_id, H5F_OBJ_ALL | H5F_OBJ_LOCAL );
        if(obj_cnt > 1) //file is still open
        {
            logit<<"**** did not close total objects "<<obj_cnt<<std::endl;
        }


        H5Fclose(_cur_file_id);
        _cur_file_id = -1;
    }
    else
    {
        logit<<"Warning: could not close file because none is open"<<std::endl;
        return false;
    }
    _cur_filename = "";
    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_spectra_volume(const std::string path,
                                  data_struct::xrf::Spectra_Volume * spectra_volume,
                                  real_t energy_offset,
                                  real_t energy_slope,
                                  real_t energy_quad,
                                  size_t row_idx_start,
                                  int row_idx_end,
                                  size_t col_idx_start,
                                  int col_idx_end)
{
    std::lock_guard<std::mutex> lock(_mutex);

    logit<<std::endl;

    if(_cur_file_id < 0)
    {
        logit << "Error: hdf5 file was never initialized. Call start_save_seq() before this function." << std::endl;
        return false;
    }

//herr_t (*old_func)(void*);
//void *old_client_data;
//hid_t error_stack = H5Eget_current_stack();
//H5Eget_auto2(error_stack, &old_func, &old_client_data);

//H5Eset_auto2(error_stack, NULL, NULL);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    hid_t    dset_id, spec_grp_id, int_spec_grp_id, dataspace_id, memoryspace_id, memoryspace_time_id, dataspace_time_id, file_time_id, filespace_id, status, maps_grp_id, dcpl_id;
    hid_t   dset_rt_id, dset_lt_id, incnt_dset_id, outcnt_dset_id;
    herr_t   error;

    hsize_t chunk_dims[3];
    hsize_t dims_out[3];
    hsize_t offset[3];
    hsize_t count[3];
    hsize_t dims_time_out[2];
    hsize_t offset_time[2];
    hsize_t count_time[2];

    if (row_idx_end < row_idx_start || row_idx_end > spectra_volume->rows() -1)
    {
        row_idx_end = spectra_volume->rows();
    }
    if (col_idx_end < col_idx_start || col_idx_end > spectra_volume->cols() -1)
    {
        col_idx_end = spectra_volume->cols();
    }

    //get one element
    data_struct::xrf::Fit_Element_Map* element;

    //H5T_FLOAT
    dims_out[0] = spectra_volume->samples_size();
    dims_out[1] = spectra_volume->rows();
    dims_out[2] = spectra_volume->cols();
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    count[0] = dims_out[0];
    count[1] = 1;
    count[2] = 1;
    chunk_dims[0] = dims_out[0];
    chunk_dims[1] = 1;
    chunk_dims[2] = 1;


    dims_time_out[0] = spectra_volume->rows();
    dims_time_out[1] = spectra_volume->cols();
    offset_time[0] = 0;
    offset_time[1] = 0;
    count_time[0] = 1;
    count_time[1] = 1;


    memoryspace_id = H5Screate_simple(3, count, NULL);
    filespace_id = H5Screate_simple(3, dims_out, NULL);

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 3, chunk_dims);
    H5Pset_deflate (dcpl_id, 7);

    dataspace_id = H5Screate_simple (3, dims_out, NULL);

    memoryspace_time_id = H5Screate_simple(2, count_time, NULL);
    file_time_id = H5Screate_simple(2, dims_time_out, NULL);
    dataspace_time_id = H5Screate_simple (2, dims_time_out, NULL);

    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
    if(maps_grp_id < 0)
        maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(maps_grp_id < 0)
    {
        logit<<"Error creating group MAPS"<<std::endl;
        return false;
    }

    spec_grp_id = H5Gopen(maps_grp_id, "Spectra", H5P_DEFAULT);
    if(spec_grp_id < 0)
        spec_grp_id = H5Gcreate(maps_grp_id, "Spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(spec_grp_id < 0)
    {
        logit<<"Error creating group MAPS/Spectra"<<std::endl;
        return false;
    }
/*
    scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
    if(scalers_grp_id < 0)
        scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(scalers_grp_id < 0)
    {
        logit<<"Error creating group MAPS/Scalers"<<std::endl;
        return false;
    }
*/
    dset_id = H5Dcreate (spec_grp_id, path.c_str(), H5T_INTEL_R, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

    dset_rt_id = H5Dcreate (spec_grp_id, "Elapsed_Realtime", H5T_INTEL_R, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dset_lt_id = H5Dcreate (spec_grp_id, "Elapsed_Livetime", H5T_INTEL_R, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    incnt_dset_id = H5Dcreate (spec_grp_id, "Input_Counts", H5T_INTEL_R, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    outcnt_dset_id = H5Dcreate (spec_grp_id, "Output_Counts", H5T_INTEL_R, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Sselect_hyperslab (memoryspace_time_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);

    real_t real_time;
    real_t life_time;
    real_t in_cnt;
    real_t out_cnt;
    for(size_t row=row_idx_start; row < row_idx_end; row++)
    {
        offset[1] = row;
        offset_time[0] = row;
        for(size_t col=col_idx_start; col < col_idx_end; col++)
        {
            const data_struct::xrf::Spectra *spectra = &((*spectra_volume)[row][col]);
            offset[2] = col;
            offset_time[1] = col;
            H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

            status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&(*spectra)[0]);

            H5Sselect_hyperslab (file_time_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);

            real_time = spectra->elapsed_realtime();
            life_time = spectra->elapsed_lifetime();
            in_cnt = spectra->input_counts();
            out_cnt = spectra->output_counts();
            status = H5Dwrite (dset_rt_id, H5T_NATIVE_REAL, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&real_time);
            status = H5Dwrite (dset_lt_id, H5T_NATIVE_REAL, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&life_time);
            status = H5Dwrite (incnt_dset_id, H5T_NATIVE_REAL, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&in_cnt);
            status = H5Dwrite (outcnt_dset_id, H5T_NATIVE_REAL, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&out_cnt);
        }
    }


    H5Dclose(dset_id);
    H5Dclose(dset_rt_id);
    H5Dclose(dset_lt_id);
    H5Dclose(incnt_dset_id);
    H5Dclose(outcnt_dset_id);
    H5Sclose(memoryspace_time_id);
    H5Sclose(memoryspace_id);
    H5Sclose(file_time_id);
    H5Sclose(filespace_id);
    H5Sclose(dataspace_time_id);
    H5Sclose(dataspace_id);
    H5Pclose(dcpl_id);


    int_spec_grp_id = H5Gopen(spec_grp_id, "Integrated_Spectra", H5P_DEFAULT);
    if(int_spec_grp_id < 0)
        int_spec_grp_id = H5Gcreate(spec_grp_id, "Integrated_Spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //save quantification_standard integrated spectra
    data_struct::xrf::Spectra spectra = spectra_volume->integrate();
    count[0] = spectra.size();
    memoryspace_id = H5Screate_simple(1, count, NULL);
    filespace_id = H5Screate_simple(1, count, NULL);
    dataspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (int_spec_grp_id, "Spectra", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    offset[0] = 0;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&spectra[0]);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(filespace_id);
    H5Sclose(dataspace_id);


    //save real_time
    count[0] = 1;
    real_t save_val = spectra.elapsed_realtime();
    dataspace_id = H5Screate_simple (1, count, NULL);
    memoryspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (int_spec_grp_id, "Elapsed_Realtime", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save life_time
    save_val = spectra.elapsed_lifetime();
    dataspace_id = H5Screate_simple (1, count, NULL);
    memoryspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (int_spec_grp_id, "Elapsed_Livetime", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save input_counts
    save_val = spectra.input_counts();
    dataspace_id = H5Screate_simple (1, count, NULL);
    memoryspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (int_spec_grp_id, "Input_Counts", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save output_counts
    save_val = spectra.output_counts();
    dataspace_id = H5Screate_simple (1, count, NULL);
    memoryspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (int_spec_grp_id, "Output_Counts", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save energy vector
    std::vector<real_t> out_vec;
    data_struct::xrf::gen_energy_vector(spectra.size(), energy_offset, energy_slope, &out_vec);
    count[0] = out_vec.size();
    memoryspace_id = H5Screate_simple(1, count, NULL);
    filespace_id = H5Screate_simple(1, count, NULL);
    dataspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (spec_grp_id, "Energy", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    offset[0] = 0;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&out_vec[0]);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(filespace_id);
    H5Sclose(dataspace_id);

    //save energy calibration
    count[0] = 3;
    dataspace_id = H5Screate_simple (1, count, NULL);
    filespace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (spec_grp_id, "Energy_Calibration", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    count[0] = 1;
    memoryspace_id = H5Screate_simple(1, count, NULL);
    offset[0] = 0;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&energy_offset);
    offset[0] = 1;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&energy_slope);
    offset[0] = 2;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&energy_quad);
    H5Dclose(dset_id);
    H5Sclose(filespace_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save file version
    save_val = HDF5_SAVE_VERSION;
    dataspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (maps_grp_id, "version", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(dataspace_id);

    H5Gclose(int_spec_grp_id);
    H5Gclose(spec_grp_id);
    //H5Gclose(scalers_grp_id);
    H5Gclose(maps_grp_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logit << "mca_arr elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;

    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_element_fits(std::string path,
                                const data_struct::xrf::Fit_Count_Dict * const element_counts,
                                size_t row_idx_start,
                                int row_idx_end,
                                size_t col_idx_start,
                                int col_idx_end)
{
    std::lock_guard<std::mutex> lock(_mutex);

    logit<<std::endl;

    if(_cur_file_id < 0)
    {
        logit << "Error: hdf5 file was never initialized. Call start_save_seq() before this function." << std::endl;
        return false;
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    std::string xrf_grp_name = "XRF_Analyzed";
    hid_t   dset_id, dset_ch_id;
    hid_t   memoryspace, filespace, dataspace_id, dataspace_ch_id, dataspace_ch_off_id;
    hid_t   filetype, memtype, status;
    hid_t   dcpl_id;
    hid_t   xrf_grp_id, fit_grp_id, maps_grp_id;
    herr_t   error;

    //hid_t q_dataspace_id, q_dataspace_ch_id, q_memoryspace_id, q_filespace_id, q_dset_id, q_dset_ch_id, q_grp_id, q_fit_grp_id;

    dset_id = -1;
    dset_ch_id = -1;
    hsize_t dims_out[3];
    hsize_t offset[1] = {0};
    hsize_t offset2[1] = {0};
    hsize_t offset_3d[3] = {0, 0, 0};
    hsize_t count[1] = {1};
    hsize_t count_3d[3] = {1, 1, 1};
    hsize_t chunk_dims[3];

    //hsize_t q_dims_out[2];
    //fix this
    for(const auto& iter : *element_counts)
    {
		dims_out[1] = iter.second.rows();
		dims_out[2] = iter.second.cols();
        break;
    }
    //H5T_FLOAT
    dims_out[0] = element_counts->size();
    offset_3d[0] = 0;
    offset_3d[1] = 0;
    offset_3d[2] = 0;
    count_3d[0] = 1;
    count_3d[1] = dims_out[1];
    count_3d[2] = dims_out[2];
    chunk_dims[0] = 1;
    chunk_dims[1] = dims_out[1]; 
    chunk_dims[2] = dims_out[2];

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 3, chunk_dims);
    H5Pset_deflate (dcpl_id, 7);

    memoryspace = H5Screate_simple(3, count_3d, NULL);
    filespace = H5Screate_simple(3, dims_out, NULL);

    dataspace_id = H5Screate_simple (3, dims_out, NULL);

    dataspace_ch_id = H5Screate_simple (1, dims_out, NULL);
    dataspace_ch_off_id = H5Screate_simple (1, dims_out, NULL);

    H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);

    maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
    if(maps_grp_id < 0)
        maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(maps_grp_id < 0)
    {
        logit<<"Error creating group 'MAPS'"<<std::endl;
        return false;
    }

    xrf_grp_id = H5Gopen(maps_grp_id, xrf_grp_name.c_str(), H5P_DEFAULT);
    if(xrf_grp_id < 0)
        xrf_grp_id = H5Gcreate(maps_grp_id, xrf_grp_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(xrf_grp_id < 0)
    {
        logit<<"Error creating group MAPS/"<<xrf_grp_name<<std::endl;
        return false;
    }

    fit_grp_id = H5Gopen(xrf_grp_id, path.c_str(), H5P_DEFAULT);
    if(fit_grp_id < 0)
        fit_grp_id = H5Gcreate(xrf_grp_id, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(fit_grp_id < 0)
    {
        logit<<"Error creating group MAPS/"<<xrf_grp_name<<"/"<<path<<std::endl;
        return false;
    }

    dset_id = H5Dopen (fit_grp_id, "Counts_Per_Sec", H5P_DEFAULT);
    if(dset_id < 0)
        dset_id = H5Dcreate (fit_grp_id, "Counts_Per_Sec", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    if(dset_id < 0)
    {
        logit<<"Error creating dataset MAPS/"<<xrf_grp_name<<"/"<<path<<"/Counts_Per_Sec"<<std::endl;
        return false;
    }

    //filetype = H5Tcopy (H5T_FORTRAN_S1);
    filetype = H5Tcopy (H5T_C_S1);
    H5Tset_size (filetype, 256);
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, 255);

    dset_ch_id = H5Dopen (fit_grp_id, "Channel_Names", H5P_DEFAULT);
    if(dset_ch_id < 0)
        dset_ch_id = H5Dcreate (fit_grp_id, "Channel_Names", filetype, dataspace_ch_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(dset_ch_id < 0)
    {
        logit<<"Error creating dataset MAPS/"<<xrf_grp_name<<"/"<<path<<"/Channel_Names"<<std::endl;
        return false;
    }

/*
   if (row_idx_end < row_idx_start || row_idx_end > spectra_volume->rows() -1)
    {
        row_idx_end = spectra_volume->rows();
    }
    if (col_idx_end < col_idx_start || col_idx_end > spectra_volume->cols() -1)
    {
        col_idx_end = spectra_volume->cols();
    }
*/

    //create save ordered vector by element Z number with K , L, M lines
    std::vector<std::string> element_lines;
    for (std::string el_name : data_struct::xrf::Element_Symbols)
    {
        element_lines.push_back(el_name);
    }
    for (std::string el_name : data_struct::xrf::Element_Symbols)
    {
        element_lines.push_back(el_name+"_L");
    }
    for (std::string el_name : data_struct::xrf::Element_Symbols)
    {
        element_lines.push_back(el_name+"_M");
    }

    element_lines.push_back(data_struct::xrf::STR_COHERENT_SCT_AMPLITUDE);
    element_lines.push_back(data_struct::xrf::STR_COMPTON_AMPLITUDE);
    element_lines.push_back(data_struct::xrf::STR_NUM_ITR);

    int i=0;
    //save by element Z order
    //for(const auto& iter : *element_counts)
    for (std::string el_name : element_lines)
    {
        char tmp_char[255] = {0};
        if(element_counts->count(el_name) < 1 )
        {
            continue;
        }
        offset[0] = i;
        offset_3d[0] = i;

        H5Sselect_hyperslab (dataspace_ch_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Sselect_hyperslab (dataspace_ch_off_id, H5S_SELECT_SET, offset2, NULL, count, NULL);

        el_name.copy(tmp_char, 254);

        status = H5Dwrite (dset_ch_id, memtype, dataspace_ch_off_id, dataspace_ch_id, H5P_DEFAULT, (void*)tmp_char);

        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset_3d, NULL, chunk_dims, NULL);

        status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace, filespace, H5P_DEFAULT, (void*)element_counts->at(el_name).data());
        
        i++;
    }

    H5Dclose(dset_id);
    H5Dclose(dset_ch_id);
    H5Sclose(memoryspace);
    H5Sclose(filespace);
    H5Sclose(dataspace_ch_off_id);
    H5Sclose(dataspace_ch_id);
    H5Tclose(filetype);
    H5Tclose(memtype);
    H5Pclose(dcpl_id);    
    H5Sclose(dataspace_id);
    H5Gclose(fit_grp_id);
    H5Gclose(xrf_grp_id);
    H5Gclose(maps_grp_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;


    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_quantification(data_struct::xrf::Quantification_Standard * quantification_standard,
                                  size_t row_idx_start,
                                  int row_idx_end,
                                  size_t col_idx_start,
                                  int col_idx_end)
{
    std::lock_guard<std::mutex> lock(_mutex);

    logit<<std::endl;

//hid_t error_stack = H5Eget_current_stack();
//H5Eset_auto2(error_stack, NULL, NULL);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    hid_t    dset_id, memoryspace_id, filespace_id, dataspace_id, filetype, dataspace_ch_id, memtype,  dset_ch_id, q_int_spec_grp_id;
    hid_t   memtype_label, filetype_label, q_memoryspace_label_id, q_dataspace_label_id;
    herr_t   error;

    hid_t q_dataspace_id, q_memoryspace_id, q_filespace_id, q_dset_id, q_grp_id, q_fit_grp_id, maps_grp_id, status, scalers_grp_id, xrf_fits_grp_id;
    hid_t dset_labels_id;
    hsize_t offset[3];
    hsize_t count[3];

    hsize_t q_dims_out[2];

    if(_cur_file_id < 0)
    {
        logit << "Error: hdf5 file was never initialized. Call start_save_seq() before this function." << std::endl;
        return false;
    }


    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    count[0] = 1;
    count[1] = 0;
    count[2] = 0;
//    chunk_dims[0] = 1;
//    chunk_dims[1] = 1;
//    chunk_dims[2] = dims_out[2];

//    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
//    H5Pset_chunk(dcpl_id, 3, chunk_dims);
//    H5Pset_deflate (dcpl_id, 7);

//    memoryspace = H5Screate_simple(3, count, NULL);
//    filespace = H5Screate_simple(3, dims_out, NULL);

//    dataspace_id = H5Screate_simple (3, dims_out, NULL);

//    dataspace_ch_id = H5Screate_simple (1, dims_out, NULL);


//    H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
    if(maps_grp_id < 0)
        maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(maps_grp_id < 0)
    {
        logit<<"Error creating group 'MAPS'"<<std::endl;
        return false;
    }

    filetype = H5Tcopy (H5T_FORTRAN_S1);
    H5Tset_size (filetype, 256);
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, 255);


    filetype_label = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype_label, 10);
    memtype_label = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(memtype_label, 10);


    //--                        save calibration curve                  --

    if( quantification_standard != nullptr && quantification_standard->processed() == true)
    {

        q_grp_id = H5Gopen(maps_grp_id, "Quantification", H5P_DEFAULT);
        if(q_grp_id < 0)
            q_grp_id = H5Gcreate(maps_grp_id, "Quantification", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(q_grp_id < 0)
        {
            logit<<"Error creating group MAPS/Quantification"<<std::endl;
            return false;
        }

        scalers_grp_id = H5Gopen(q_grp_id, "Scalers", H5P_DEFAULT);
        if(scalers_grp_id < 0)
            scalers_grp_id = H5Gcreate(q_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(scalers_grp_id < 0)
        {
            logit<<"Error creating group MAPS/Quantification/Scalers"<<std::endl;
            return false;
        }

        xrf_fits_grp_id = H5Gopen(q_grp_id, "XRF_Analyzed", H5P_DEFAULT);
        if(xrf_fits_grp_id < 0)
            xrf_fits_grp_id = H5Gcreate(q_grp_id, "XRF_Analyzed", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(xrf_fits_grp_id < 0)
        {
            logit<<"Error creating group MAPS/Quantification/XRF_Analyzed"<<std::endl;
            return false;
        }

        //save quantification_standard element weights
        const std::unordered_map<std::string, data_struct::xrf::Element_Quant> e_weights = quantification_standard->element_weights();
        count[0] = e_weights.size();
        memoryspace_id = H5Screate_simple(1, count, NULL);
        filespace_id = H5Screate_simple(1, count, NULL);
        dataspace_id = H5Screate_simple (1, count, NULL);
        dataspace_ch_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (q_grp_id, "Element_Weights", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        dset_ch_id = H5Dcreate (q_grp_id, "Element_Weights_Names", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        offset[0] = 0;
        count[0] = 1;
        H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        int offset_idx = 0;
        for(auto itr: e_weights)
        {
            offset[0] = offset_idx;
            H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

            status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&(itr.second.weight));
            char tmp_char[255] = {0};
            itr.first.copy(tmp_char, 254);
            status = H5Dwrite (dset_ch_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)tmp_char);

            offset_idx++;
        }
        H5Dclose(dset_id);
        H5Dclose(dset_ch_id);
        H5Sclose(memoryspace_id);
        H5Sclose(filespace_id);
        H5Sclose(dataspace_id);
        H5Sclose(dataspace_ch_id);


        q_int_spec_grp_id = H5Gopen(q_grp_id, "Integrated_Spectra", H5P_DEFAULT);
        if(q_int_spec_grp_id < 0)
           q_int_spec_grp_id  = H5Gcreate(q_grp_id, "Integrated_Spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(q_int_spec_grp_id < 0)
        {
            logit<<"Error creating group MAPS/Quantification/Integrated_Spectra"<<std::endl;
            return false;
        }

        //save quantification_standard integrated spectra
        data_struct::xrf::Spectra spectra = quantification_standard->integrated_spectra();
        count[0] = spectra.size();
        memoryspace_id = H5Screate_simple(1, count, NULL);
        filespace_id = H5Screate_simple(1, count, NULL);
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Spectra", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        offset[0] = 0;
        H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&spectra[0]);
        H5Dclose(dset_id);
        H5Sclose(memoryspace_id);
        H5Sclose(filespace_id);
        H5Sclose(dataspace_id);


        //save standard name
        count[0] = 1;
        dataspace_id = H5Screate_simple (1, count, NULL);
        memoryspace_id = H5Screate_simple (1, count, NULL);
        dset_ch_id = H5Dcreate (q_grp_id, "Standard_Name", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        char tmp_char[255] = {0};
        quantification_standard->standard_filename().copy(tmp_char, 254);
        status = H5Dwrite (dset_ch_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);
        H5Dclose(dset_ch_id);
        H5Sclose(dataspace_id);

        //save sr_current
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (scalers_grp_id, "SR_Current", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quantification_standard->sr_current()));
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);

        //save us_ic
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (scalers_grp_id, "US_IC", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quantification_standard->US_IC()));
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);

        //save ds_ic
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (scalers_grp_id, "DS_IC", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quantification_standard->DS_IC()));
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);


        //save real_time
        real_t save_val = spectra.elapsed_realtime();
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Elapsed_Realtime", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);

        //save life_time
        save_val = spectra.elapsed_lifetime();
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Elapsed_Livetime", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        H5Dclose(dset_id);
        H5Sclose(memoryspace_id);
        H5Sclose(dataspace_id);

        //save input counts
        save_val = spectra.input_counts();
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Input_Counts", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        H5Dclose(dset_id);
        H5Sclose(memoryspace_id);
        H5Sclose(dataspace_id);

        //save output counts
        save_val = spectra.output_counts();
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Output_Counts", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        H5Dclose(dset_id);
        H5Sclose(memoryspace_id);
        H5Sclose(dataspace_id);

        //save calibration curves
        //data_struct::xrf::Quantifiers quantifiers = quantification_standard->calibration_curves.at(path);
        //auto shell_itr = quantification_standard->_calibration_curves.begin();

        q_dims_out[0] = 3;// shells K, L, and M
        q_dims_out[1] = 93; // elements 1 - 92

        offset[0] = 0;
        offset[1] = 0;
        offset[2] = 0;

        count[0] = 1;
        count[1] = 1;
        count[2] = 0;

        q_memoryspace_label_id = H5Screate_simple(2, count, NULL);

        count[0] = 1;
        count[1] = q_dims_out[1];
        count[2] = 0;

        q_memoryspace_id = H5Screate_simple(2, count, NULL);
        q_filespace_id = H5Screate_simple(2, q_dims_out, NULL);
        q_dataspace_label_id = H5Screate_simple (2, q_dims_out, NULL);
        q_dataspace_id = H5Screate_simple (2, q_dims_out, NULL);
        //q_dataspace_ch_id = H5Screate_simple (1, q_dims_out, NULL);


        for (auto& qitr: quantification_standard->calibration_curves)
        {

            q_fit_grp_id = H5Gopen(xrf_fits_grp_id, qitr.first.c_str(), H5P_DEFAULT);
            if(q_fit_grp_id < 0)
                q_fit_grp_id = H5Gcreate(xrf_fits_grp_id, qitr.first.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if(q_fit_grp_id < 0)
            {
                logit<<"Error creating group MAPS/Quantification/"<<qitr.first<<std::endl;
                return false;
            }

            //save quantification_standard counts
            const std::unordered_map<std::string, std::unordered_map<std::string, real_t> > all_element_counts = quantification_standard->element_counts();
            std::unordered_map<std::string, real_t> element_counts = all_element_counts.at(qitr.first);
            count[0] = element_counts.size();
            dataspace_id = H5Screate_simple (1, count, NULL);
            memoryspace_id = H5Screate_simple(1, count, NULL);
            filespace_id = H5Screate_simple(1, count, NULL);
            dset_id = H5Dcreate (q_fit_grp_id, "Counts_Per_Sec", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            dset_ch_id = H5Dcreate (q_fit_grp_id, "Channel_Names", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            //create save ordered vector by element Z number with K , L, M lines
            std::vector<std::string> element_lines;
            for (std::string el_name : data_struct::xrf::Element_Symbols)
            {
                element_lines.push_back(el_name);
            }
            for (std::string el_name : data_struct::xrf::Element_Symbols)
            {
                element_lines.push_back(el_name+"_L");
            }
            for (std::string el_name : data_struct::xrf::Element_Symbols)
            {
                element_lines.push_back(el_name+"_M");
            }

            element_lines.push_back(data_struct::xrf::STR_NUM_ITR);

            offset_idx = 0;
            offset[0] = 0;
            count[0] = 1;
            H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            //save by element Z order
            for (std::string el_name : element_lines)
            {
                if(element_counts.count(el_name) < 1 )
                {
                    continue;
                }
                offset[0] = offset_idx;
                real_t val = element_counts.at(el_name);

                H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                char tmp_char[255] = {0};
                el_name.copy(tmp_char, 254);
                status = H5Dwrite (dset_ch_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)tmp_char);
                status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&val);

                offset_idx++;
            }
            H5Dclose(dset_id);
            H5Dclose(dset_ch_id);
            H5Sclose(memoryspace_id);
            H5Sclose(filespace_id);
            H5Sclose(dataspace_id);


            dset_labels_id = H5Dopen (q_fit_grp_id, "Calibration_Curve_Labels", H5P_DEFAULT);
            if(dset_labels_id < 0)
                dset_labels_id = H5Dcreate (q_fit_grp_id, "Calibration_Curve_Labels", filetype_label, q_dataspace_label_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


            for (size_t i =0; i< qitr.second.calib_curves.size(); i++)
            {

                std::string q_dset_name = "Calibration_Curve_" + qitr.second.calib_curves[i].quantifier_name;

                q_dset_id = H5Dopen (q_fit_grp_id, q_dset_name.c_str(), H5P_DEFAULT);
                if(q_dset_id < 0)
                    q_dset_id = H5Dcreate (q_fit_grp_id, q_dset_name.c_str(), H5T_INTEL_R, q_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if(q_dset_id < 0)
                {
                    logit<<"Error creating dataset MAPS/Quantification/"<<q_dset_name<<"/"<<q_dset_name<<std::endl;
                    continue;
                }


                //for(auto& shell_itr : quant_itr.second)
                for(size_t j=0; j<3; j++)
                {
                    char label[10] = {0};
                    //int element_offset = 0;
                    //create dataset for different shell curves
                    std::vector<real_t> calibration_curve = qitr.second.calib_curves[i].shell_curves[j];
                    std::vector<std::string> calibration_curve_labels = qitr.second.calib_curves[i].shell_curves_labels[j];
                    offset[0] = j;
                    offset[1] = 0;

                    count[0] = 1;
                    count[1] = q_dims_out[1];
                    //H5Sselect_hyperslab (q_dataspace_ch_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    //H5Sselect_hyperslab (q_dataspace_ch_off_id, H5S_SELECT_SET, &offset[2], NULL, count, NULL);
                    //status = H5Dwrite (q_dset_ch_id, memtype, q_dataspace_ch_off_id, q_dataspace_ch_id, H5P_DEFAULT, (void*)(el_name.c_str()));

                    H5Sselect_hyperslab (q_filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite (q_dset_id, H5T_NATIVE_REAL, q_memoryspace_id, q_filespace_id, H5P_DEFAULT, (void*)&calibration_curve[0]);

                    for(size_t k=0; k<calibration_curve_labels.size(); k++)
                    {
                        memset(label, 0, 10);
                        calibration_curve_labels[k].copy(&label[0], 9);
                        offset[1] = k;
                        count[1] = 1;
                        status = H5Sselect_hyperslab (q_filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                        status = H5Dwrite (dset_labels_id, memtype_label, q_memoryspace_label_id, q_filespace_id, H5P_DEFAULT, (void*)&label[0]);
                    }

                }
                H5Dclose(q_dset_id);
            }

            H5Dclose(dset_labels_id);
            H5Gclose(q_fit_grp_id);
        }

        //H5Dclose(q_dset_ch_id);
        H5Sclose(q_filespace_id);
        H5Sclose(q_memoryspace_id);
        H5Sclose(q_dataspace_label_id);
        H5Sclose(q_memoryspace_label_id);
        H5Sclose(q_dataspace_id);
        //H5Sclose(q_dataspace_ch_id);
        H5Gclose(xrf_fits_grp_id);
        H5Gclose(scalers_grp_id);
        H5Gclose(q_int_spec_grp_id);
        H5Gclose(q_grp_id);



    }

    H5Tclose(filetype);
    H5Tclose(memtype);

    H5Tclose(filetype_label);
    H5Tclose(memtype_label);



//    H5Dclose(dset_id);
//    H5Dclose(dset_ch_id);
//    H5Sclose(memoryspace);
//    H5Sclose(filespace);
//    H5Sclose(dataspace_ch_off_id);
//    H5Sclose(dataspace_ch_id);
//    H5Pclose(dcpl_id);
//    H5Sclose(dataspace_id);
//    H5Gclose(xrf_grp_id);
    H5Gclose(maps_grp_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;



    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::_save_scan_meta_data(hid_t scan_grp_id, struct mda_file *mda_scalers)
{
    hid_t dataspace_id = -1, memoryspace_id = -1, filespace_id = -1;
	hid_t status;
	hid_t filetype, memtype;
    hid_t dset_id = -1;
	
	hsize_t offset[1] = { 0 };
	hsize_t count[1] = { 1 };

	bool single_row_scan = false;

    logit<<std::endl;

	try
	{

        filetype = H5Tcopy(H5T_FORTRAN_S1);
		H5Tset_size(filetype, 256);
		memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);


		memoryspace_id = H5Screate_simple(1, count, NULL);

		//save scan positions
		if (mda_scalers->scan->scan_rank > 1)
		{
			if (mda_scalers->header->data_rank == 2)
			{
				if (mda_scalers->header->dimensions[1] == 2000)
				{
					single_row_scan = true;
				}
			}

			if (single_row_scan)
			{
				count[0] = 1;
				dataspace_id = H5Screate_simple(1, count, NULL);
                dset_id = H5Dcreate(scan_grp_id, "y_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
                    logit << "Error creating dataset 'y_axis'" << std::endl;
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
					return false;
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);

                if(mda_scalers->scan->last_point == 0)
                    count[0] = 1;
                else
                    count[0] = mda_scalers->scan->last_point;
				dataspace_id = H5Screate_simple(1, count, NULL);
				filespace_id = H5Screate_simple(1, count, NULL);
                dset_id = H5Dcreate(scan_grp_id, "x_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
                    logit << "Error creating dataset 'x_axis'" << std::endl;
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
					return false;
				}
				count[0] = 1;
				for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
				{
					offset[0] = i;
					double pos = mda_scalers->scan->positioners_data[0][i];
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&pos);
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);
				H5Sclose(filespace_id);


                //save requested rows
                count[0] = 1;
                dataspace_id = H5Screate_simple(1, count, NULL);
                filespace_id = H5Screate_simple(1, count, NULL);
                dset_id = H5Dcreate(scan_grp_id, "requested_rows", H5T_INTEL_I32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (dset_id < 0)
                {
                    logit << "Error creating dataset 'requested_rows'" << std::endl;
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
                    return false;
                }
                int value = 1;
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&value);
                H5Dclose(dset_id);
                H5Sclose(dataspace_id);
                H5Sclose(filespace_id);

                //save requested cols
                count[0] = 1;
                dataspace_id = H5Screate_simple(1, count, NULL);
                filespace_id = H5Screate_simple(1, count, NULL);
                dset_id = H5Dcreate(scan_grp_id, "requested_cols", H5T_INTEL_I32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (dset_id < 0)
                {
                    logit << "Error creating dataset 'requested_cols'" << std::endl;
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
                    return false;
                }
                value = mda_scalers->header->dimensions[0];
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&value);
                H5Dclose(dset_id);
                H5Sclose(dataspace_id);
                H5Sclose(filespace_id);

			}
			else
			{
                //save y axis
                if(mda_scalers->scan->last_point == 0)
                    count[0] = 1;
                else
                    count[0] = mda_scalers->scan->last_point;
				dataspace_id = H5Screate_simple(1, count, NULL);
				filespace_id = H5Screate_simple(1, count, NULL);
                dset_id = H5Dcreate(scan_grp_id, "y_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
                    logit << "Error creating dataset 'y_axis'" << std::endl;
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
					return false;
				}
				count[0] = 1;
				for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
				{
					offset[0] = i;
					double pos = mda_scalers->scan->positioners_data[0][i];
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&pos);
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);
				H5Sclose(filespace_id);

                //save x axis
                if(mda_scalers->scan->sub_scans[0]->last_point == 0)
                    count[0] = 1;
                else
                    count[0] = mda_scalers->scan->sub_scans[0]->last_point;
				dataspace_id = H5Screate_simple(1, count, NULL);
				filespace_id = H5Screate_simple(1, count, NULL);
                dset_id = H5Dcreate(scan_grp_id, "x_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
                    logit << "Error creating dataset 'x_axis'" << std::endl;
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
					return false;
				}
				count[0] = 1;
				for (int32_t i = 0; i < mda_scalers->scan->sub_scans[0]->last_point; i++)
				{
					offset[0] = i;
					double pos = mda_scalers->scan->sub_scans[0]->positioners_data[0][i];
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&pos);
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);
				H5Sclose(filespace_id);


                //save requested rows
                count[0] = 1;
                dataspace_id = H5Screate_simple(1, count, NULL);
                filespace_id = H5Screate_simple(1, count, NULL);
                dset_id = H5Dcreate(scan_grp_id, "requested_rows", H5T_INTEL_I32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (dset_id < 0)
                {
                    logit << "Error creating dataset 'requested_rows'" << std::endl;
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
                    return false;
                }
                int value = mda_scalers->header->dimensions[0];
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&value);
                H5Dclose(dset_id);
                H5Sclose(dataspace_id);
                H5Sclose(filespace_id);

                //save requested cols
                count[0] = 1;
                dataspace_id = H5Screate_simple(1, count, NULL);
                filespace_id = H5Screate_simple(1, count, NULL);
                dset_id = H5Dcreate(scan_grp_id, "requested_cols", H5T_INTEL_I32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (dset_id < 0)
                {
                    logit << "Error creating dataset 'requested_cols'" << std::endl;
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
                    return false;
                }
                value = mda_scalers->header->dimensions[1];
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&value);
                H5Dclose(dset_id);
                H5Sclose(dataspace_id);
                H5Sclose(filespace_id);

			}

		}

		//save write date
		count[0] = 1;

		dset_id = H5Dcreate(scan_grp_id, "scan_time_stamp", filetype, memoryspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id < 0)
		{
            logit << "Error creating dataset 'scan_time_stamp'" << std::endl;
			return false;
		}
		if (mda_scalers->scan->time != nullptr)
		{
			std::string str_time = std::string(mda_scalers->scan->time);
            char tmp_char[255] = {0};
            str_time.copy(tmp_char, 254);
            status = H5Dwrite(dset_id, memtype, memoryspace_id, memoryspace_id, H5P_DEFAULT, (void*)tmp_char);
		}
        if(dset_id > -1)
        {
            H5Dclose(dset_id);
            dset_id = -1;
        }
		dset_id = H5Dcreate(scan_grp_id, "name", filetype, memoryspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id < 0)
		{
            logit << "Error creating dataset 'name'" << std::endl;
			return false;
		}
		if(mda_scalers->scan->name != nullptr)
		{
			std::string str_name = std::string(mda_scalers->scan->name);
            char tmp_char[255] = {0};
            str_name.copy(tmp_char, 254);
            status = H5Dwrite(dset_id, memtype, memoryspace_id, memoryspace_id, H5P_DEFAULT, (void*)tmp_char);
		}

        H5Tclose(filetype);
        H5Tclose(memtype);

        if(dset_id > -1)
        {
            H5Dclose(dset_id);
            dset_id = -1;
        }
        if(memoryspace_id > -1)
        {
            H5Sclose(memoryspace_id);
            memoryspace_id = -1;
        }
	}
	catch (...)
	{
        logit << "Error saving MAPS/Scan meta data" << std::endl;
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::_save_extras(hid_t scan_grp_id, struct mda_file *mda_scalers)
{

    hid_t filespace_id = -1, filespace_id2 = -1, filespace_id3 = -1, filespace_id4 = -1;
    hid_t memoryspace_id = -1, memoryspace_id2 = -1, memoryspace_id3 = -1, memoryspace_id4 = -1;
    hid_t status = -1;
	hid_t filetype, memtype;
    hid_t dset_desc_id = -1, dset_unit_id = -1, dset_id = -1, dset_val_id = -1;
    hid_t extra_grp_id = -1;
	
	hsize_t offset[1] = { 0 };
	hsize_t count[1] = { 1 };

    logit<<std::endl;

	if (mda_scalers->extra == nullptr)
	{
		return false;
	}

	try
	{

        filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, 256);
		memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);

        //logit<<"extra_pvs group"<<std::endl;
		extra_grp_id = H5Gopen(scan_grp_id, "Extra_PVs", H5P_DEFAULT);
		if (extra_grp_id < 0)
			extra_grp_id = H5Gcreate(scan_grp_id, "Extra_PVs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (extra_grp_id < 0)
		{
            logit << "Error creating group MAPS/Scan/Extra_PVs" << std::endl;
			return false;
		}

		//save extra pv's
		count[0] = (size_t)mda_scalers->extra->number_pvs;
		filespace_id = H5Screate_simple(1, count, NULL);
        //logit<<"names dataset"<<std::endl;
        dset_id = H5Dcreate(extra_grp_id, "Names", filetype, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id < 0)
		{
            logit << "Error creating dataset MAPS/Scan/Extra_PVs/Names" << std::endl;
            if(filespace_id > -1)
            {
                H5Sclose(filespace_id);
                filespace_id = -1;
            }
            if(extra_grp_id > -1)
            {
                H5Gclose(extra_grp_id);
                extra_grp_id = -1;
            }
            return false;
		}

        filespace_id2 = H5Screate_simple(1, count, NULL);
        //logit<<"values dataset"<<std::endl;
        dset_val_id = H5Dcreate(extra_grp_id, "Values", filetype, filespace_id2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_val_id < 0)
		{
            logit << "Error creating dataset MAPS/Scan/Extra_PVs/Values" << std::endl;
            if(filespace_id > -1)
            {
                H5Sclose(filespace_id);
                filespace_id = -1;
            }
            if(filespace_id2 > -1)
            {
                H5Sclose(filespace_id2);
                filespace_id2 = -1;
            }
            if(extra_grp_id > -1)
            {
                H5Gclose(extra_grp_id);
                extra_grp_id = -1;
            }
            return false;
		}

        filespace_id3 = H5Screate_simple(1, count, NULL);
        //logit<<"desc dataset"<<std::endl;
        dset_desc_id = H5Dcreate(extra_grp_id, "Description", filetype, filespace_id3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_desc_id < 0)
		{
            logit << "Error creating dataset MAPS/Scan/Extra_PVs/Description" << std::endl;
            if(filespace_id > -1)
            {
                H5Sclose(filespace_id);
                filespace_id = -1;
            }
            if(filespace_id2 > -1)
            {
                H5Sclose(filespace_id2);
                filespace_id2 = -1;
            }
            if(filespace_id3 > -1)
            {
                H5Sclose(filespace_id3);
                filespace_id3 = -1;
            }
            if(extra_grp_id > -1)
            {
                H5Gclose(extra_grp_id);
                extra_grp_id = -1;
            }
            return false;
		}

        filespace_id4 = H5Screate_simple(1, count, NULL);
        //logit<<"unit dataset"<<std::endl;
        dset_unit_id = H5Dcreate(extra_grp_id, "Unit", filetype, filespace_id4, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_unit_id < 0)
		{
            logit << "Error creating dataset MAPS/Scan/Extra_PVs/Unit" << std::endl;
            if(filespace_id > -1)
            {
                H5Sclose(filespace_id);
                filespace_id = -1;
            }
            if(filespace_id2 > -1)
            {
                H5Sclose(filespace_id2);
                filespace_id2 = -1;
            }
            if(filespace_id3 > -1)
            {
                H5Sclose(filespace_id3);
                filespace_id3 = -1;
            }
            if(filespace_id4 > -1)
            {
                H5Sclose(filespace_id4);
                filespace_id4 = -1;
            }
            if(extra_grp_id > -1)
            {
                H5Gclose(extra_grp_id);
                extra_grp_id = -1;
            }
            return false;
		}

        count[0] = 1;

		std::string str_val;
		short* s_val;
		int* i_val;
		float* f_val;
		double* d_val;

        //logit<<"memspaces"<<std::endl;
		count[0] = 1;
		memoryspace_id = H5Screate_simple(1, count, NULL);
		H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
		memoryspace_id2 = H5Screate_simple(1, count, NULL);
		H5Sselect_hyperslab(memoryspace_id2, H5S_SELECT_SET, offset, NULL, count, NULL);
		memoryspace_id3 = H5Screate_simple(1, count, NULL);
		H5Sselect_hyperslab(memoryspace_id3, H5S_SELECT_SET, offset, NULL, count, NULL);
		memoryspace_id4 = H5Screate_simple(1, count, NULL);
		H5Sselect_hyperslab(memoryspace_id4, H5S_SELECT_SET, offset, NULL, count, NULL);

		for (int16_t i = 0; i < mda_scalers->extra->number_pvs; i++)
		{
			offset[0] = i;
			struct mda_pv * pv = mda_scalers->extra->pvs[i];
            if(pv == nullptr)
            {
                continue;
            }
			switch (pv->type)
			{

			case EXTRA_PV_STRING:
				str_val = std::string(pv->values);
				break;
				//case EXTRA_PV_INT8:

				//    break;
			case EXTRA_PV_INT16:
				s_val = (short*)pv->values;
				str_val = std::to_string(*s_val);
				break;
			case EXTRA_PV_INT32:
				i_val = (int*)pv->values;
				str_val = std::to_string(*i_val);
				break;
			case EXTRA_PV_FLOAT:
				f_val = (float*)pv->values;
				str_val = std::to_string(*f_val);
				break;
			case EXTRA_PV_DOUBLE:
				d_val = (double*)pv->values;
				str_val = std::to_string(*d_val);
				break;

			}

			H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(filespace_id2, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(filespace_id3, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(filespace_id4, H5S_SELECT_SET, offset, NULL, count, NULL);

			if (pv->name != nullptr)
			{
				//need this becuase it crashes on windows if we just pass the char *
                std::string name = std::string(pv->name);
                char tmp_char[255] = {0};
                name.copy(tmp_char, 254);
                status = H5Dwrite(dset_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)tmp_char);
			}

            //logit<<"save val : "<<str_val<<" "<<str_val.length()<<std::endl;
            char tmp_char2[255] = {0};
            str_val.copy(tmp_char2, 254);
            status = H5Dwrite(dset_val_id, memtype, memoryspace_id2, filespace_id2, H5P_DEFAULT, (void*)tmp_char2);
			if (pv->description != nullptr)
			{
                std::string description = std::string(pv->description);
                char tmp_char[255] = {0};
                description.copy(tmp_char, 254);
                status = H5Dwrite(dset_desc_id, memtype, memoryspace_id3, filespace_id3, H5P_DEFAULT, (void*)tmp_char);
			}

			if (pv->unit != nullptr)
			{
                std::string unit = std::string(pv->unit);
                char tmp_char[255] = {0};
                unit.copy(tmp_char, 254);
                status = H5Dwrite(dset_unit_id, memtype, memoryspace_id4, filespace_id4, H5P_DEFAULT, (void*)tmp_char);
			}
		}

        H5Tclose(filetype);
        H5Tclose(memtype);

        if(dset_id > -1)
        {
            H5Dclose(dset_id);
            dset_id = -1;
        }
        if(dset_val_id > -1)
        {
            H5Dclose(dset_val_id);
            dset_val_id = -1;
        }
        if(dset_desc_id > -1)
        {
            H5Dclose(dset_desc_id);
            dset_desc_id = -1;
        }
        if(dset_unit_id > -1)
        {
            H5Dclose(dset_unit_id);
            dset_unit_id = -1;
        }
        if(memoryspace_id > -1)
        {
            H5Sclose(memoryspace_id);
            memoryspace_id = -1;
        }
        if(memoryspace_id2 > -1)
        {
            H5Sclose(memoryspace_id2);
            memoryspace_id2 = -1;
        }
        if(memoryspace_id3 > -1)
        {
            H5Sclose(memoryspace_id3);
            memoryspace_id3 = -1;
        }
        if(memoryspace_id4 > -1)
        {
            H5Sclose(memoryspace_id4);
            memoryspace_id4 = -1;
        }
        if(filespace_id > -1)
        {
            H5Sclose(filespace_id);
            filespace_id = -1;
        }
        if(filespace_id2 > -1)
        {
            H5Sclose(filespace_id2);
            filespace_id2 = -1;
        }
        if(filespace_id3 > -1)
        {
            H5Sclose(filespace_id3);
            filespace_id3 = -1;
        }
        if(filespace_id4 > -1)
        {
            H5Sclose(filespace_id4);
            filespace_id4 = -1;
        }
        if(extra_grp_id > -1)
        {
            H5Gclose(extra_grp_id);
            extra_grp_id = -1;
        }

	}
	catch (...)
	{
        if(dset_id > -1)
        {
            H5Dclose(dset_id);
            dset_id = -1;
        }
        if(dset_val_id > -1)
        {
            H5Dclose(dset_val_id);
            dset_val_id = -1;
        }
        if(dset_desc_id > -1)
        {
            H5Dclose(dset_desc_id);
            dset_desc_id = -1;
        }
        if(dset_unit_id > -1)
        {
            H5Dclose(dset_unit_id);
            dset_unit_id = -1;
        }
        if(memoryspace_id > -1)
        {
            H5Sclose(memoryspace_id);
            memoryspace_id = -1;
        }
        if(memoryspace_id2 > -1)
        {
            H5Sclose(memoryspace_id2);
            memoryspace_id2 = -1;
        }
        if(memoryspace_id3 > -1)
        {
            H5Sclose(memoryspace_id3);
            memoryspace_id3 = -1;
        }
        if(memoryspace_id4 > -1)
        {
            H5Sclose(memoryspace_id4);
            memoryspace_id4 = -1;
        }
        if(filespace_id > -1)
        {
            H5Sclose(filespace_id);
            filespace_id = -1;
        }
        if(filespace_id > -1)
        {
            H5Sclose(filespace_id);
            filespace_id = -1;
        }
        if(filespace_id2 > -1)
        {
            H5Sclose(filespace_id2);
            filespace_id2 = -1;
        }
        if(filespace_id3 > -1)
        {
            H5Sclose(filespace_id3);
            filespace_id3 = -1;
        }
        if(filespace_id4 > -1)
        {
            H5Sclose(filespace_id4);
            filespace_id4 = -1;
        }
        if(extra_grp_id > -1)
        {
            H5Gclose(extra_grp_id);
            extra_grp_id = -1;
        }
        logit << "Error creating group MAPS/Scan/Extra_PVs" << std::endl;
		return false;
	}
    //logit<<"done"<<std::endl;
	return true;
}

bool HDF5_IO::_save_scalers(hid_t maps_grp_id, struct mda_file *mda_scalers, size_t detector_num, data_struct::xrf::Params_Override * params_override, bool hasNetcdf)
{
    hid_t dataspace_id = -1, memoryspace_id = -1, filespace_id = -1, filespace_name_id = -1, memoryspace_str_id = -1;
    hid_t filetype, memtype;
    hid_t dset_names_id = -1;
    hid_t scalers_grp_id = -1;
    hid_t dcpl_id = -1, status;

    hsize_t offset[1] = { 0 };
    hsize_t count[1] = { 1 };

	hsize_t offset_2d[2] = { 0, 0 };
	hsize_t count_2d[2] = { 1, 1 };

    hsize_t offset_3d[3] = { 0, 0, 0 };
    hsize_t count_3d[3] = { 1, 1, 1 };

    int mda_time_scaler_idx = -1;
    real_t time_scaler_clock = 1.0;

    MDA_IO mda_io;

    bool single_row_scan = false;

    Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> scaler_mat, abs_cfg_mat, H_dpc_cfg_mat, V_dpc_cfg_mat, dia1_dpc_cfg_mat, dia2_dpc_cfg_mat;

    logit<<std::endl;

    //don't save these scalers
    std::list<std::string> ignore_scaler_strings = { "ELT1", "ERT1", "ICR1", "OCR1" };
    std::list<struct scaler_struct> scalers;

    hid_t dset_cps_id = -1;

	memoryspace_str_id = H5Screate_simple(1, count, NULL);

    try
    {

        if(params_override->time_scaler_clock.length() > 0)
        {
            time_scaler_clock = std::stod(params_override->time_scaler_clock);
        }

        filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, 256);
        memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);

        count[0] = 1;

        real_t val;
        bool save_cfg_abs = false;
        if (params_override != nullptr && mda_scalers->scan->scan_rank > 1)
        {
            int us_ic_idx = -1;
            int ds_ic_idx = -1;
            int cfg_2_idx = -1;
            int cfg_3_idx = -1;
            int cfg_4_idx = -1;
            int cfg_5_idx = -1;

            int hdf_idx = 0;
            for (auto itr : params_override->scaler_pvs)
            {
                //don't save ELT1, ERT1, ICR1, OCR1. these are saved elsewhere
                std::list<std::string>::iterator s_itr = std::find(ignore_scaler_strings.begin(), ignore_scaler_strings.end(), itr.first);
                if (s_itr != ignore_scaler_strings.end())
                    continue;

                int mda_idx = mda_io.find_scaler_index(mda_scalers, itr.second, val);
                scalers.push_back(scaler_struct(itr.first, mda_idx, hdf_idx, false));
                hdf_idx++;
                if (mda_idx > -1)
                {
                    if (itr.first == "US_IC")
                        us_ic_idx = mda_idx;
                    else if (itr.first == "DS_IC")
                        ds_ic_idx = mda_idx;
                    else if (itr.first == "CFG_2")
                        cfg_2_idx = mda_idx;
                    else if (itr.first == "CFG_3")
                        cfg_3_idx = mda_idx;
                    else if (itr.first == "CFG_4")
                        cfg_4_idx = mda_idx;
                    else if (itr.first == "CFG_5")
                        cfg_5_idx = mda_idx;
                }
            }
            for (auto itr : params_override->time_normalized_scalers)
            {
                //don't save ELT1, ERT1, ICR1, OCR1. these are saved elsewhere
                std::list<std::string>::iterator s_itr = std::find(ignore_scaler_strings.begin(), ignore_scaler_strings.end(), itr.first);
                if (s_itr != ignore_scaler_strings.end())
                    continue;

                int mda_idx = mda_io.find_scaler_index(mda_scalers, itr.second, val);
                if (mda_idx > -1)
                {
                    bool found_scaler = false;
                    for(auto& subitr : scalers)
                    {
                        if(subitr.hdf_name == itr.first)
                        {
                            subitr.mda_idx = mda_idx;
                            subitr.normalize_by_time = true;
                            found_scaler = true;
                            break;
                        }
                    }
                    if(found_scaler == false)
                    {
                        scalers.push_back(scaler_struct(itr.first, mda_idx, hdf_idx, true));
                        hdf_idx++;
                    }
                    if (itr.first == "US_IC")
                        us_ic_idx = mda_idx;
                    else if (itr.first == "DS_IC")
                        ds_ic_idx = mda_idx;
                    else if (itr.first == "CFG_2")
                        cfg_2_idx = mda_idx;
                    else if (itr.first == "CFG_3")
                        cfg_3_idx = mda_idx;
                    else if (itr.first == "CFG_4")
                        cfg_4_idx = mda_idx;
                    else if (itr.first == "CFG_5")
                        cfg_5_idx = mda_idx;
                }
            }

            //search for time scaler index
            mda_time_scaler_idx = mda_io.find_scaler_index(mda_scalers, params_override->time_scaler, val);

            scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
            if (scalers_grp_id < 0)
                scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (scalers_grp_id < 0)
            {
                logit << "Error creating group MAPS/Scalers" << std::endl;
                return false;
            }

            _save_amps(scalers_grp_id, mda_scalers, params_override);

            if (scalers.size() > 0)
            {

                if (mda_scalers->header->data_rank == 2)
                {
                    if(hasNetcdf)
                    {
                        count_3d[0] = 1;
                        if(mda_scalers->scan->last_point == 0)
                            count_3d[1] = 1;
                        else
                            count_3d[1] = mda_scalers->scan->last_point;
                        if(mda_scalers->scan->sub_scans[0]->last_point == 0)
                            count_3d[2] = 1;
                        else
                            count_3d[2] = mda_scalers->scan->sub_scans[0]->last_point;
                    }
                    else
                    {
                        if(mda_scalers->header->dimensions[1] == 2000)
                        {
                            count_3d[0] = 1;
                            count_3d[1] = 1;
                            if(mda_scalers->scan->last_point == 0)
                                count_3d[2] = 1;
                            else
                                count_3d[2] = mda_scalers->scan->last_point;
                            single_row_scan = true;
                        }
                        else
                        {
                            logit<<"Error, unknown or bad mda file"<<std::endl;
                        }
                    }
                }
                else if (mda_scalers->header->data_rank == 3)
                {
                    count_3d[0] = 1;
                    if(mda_scalers->scan->last_point == 0)
                        count_3d[1] = 1;
                    else
                        count_3d[1] = mda_scalers->scan->last_point;
                    if(mda_scalers->scan->sub_scans[0]->last_point == 0)
                        count_3d[2] = 1;
                    else
                        count_3d[2] = mda_scalers->scan->sub_scans[0]->last_point;
                }
                else
                {
                    logit << "Unsupported rank " << mda_scalers->header->data_rank << " . Skipping scalers" << std::endl;
                    return false;
                }
                dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
                H5Pset_chunk(dcpl_id, 3, count_3d);
                H5Pset_deflate(dcpl_id, 7);

                if (us_ic_idx > -1 && ds_ic_idx > -1 && cfg_2_idx > -1 && cfg_2_idx > -1 && cfg_2_idx > -1 && cfg_2_idx > -1)
                {
                    count_3d[0] = scalers.size() + 6; //abs_ic, abs_cfg, H_dpc_cfg, V_dpc_cfg, dia1_dpc_cfg, dia2_dpc_cfg
                    save_cfg_abs = true;
                }
                else
                {
                    count_3d[0] = scalers.size();
                }
                dataspace_id = H5Screate_simple(3, count_3d, NULL);
                filespace_id = H5Screate_simple(3, count_3d, NULL);

                count[0] = count_3d[0];
                filespace_name_id = H5Screate_simple(1, count, NULL);

                dset_cps_id = H5Dcreate(scalers_grp_id, "Values", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
                dset_names_id = H5Dcreate(scalers_grp_id, "Names", filetype, filespace_name_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                count_3d[0] = 1;
                count[0] = 1;

				//save scalers
				if (single_row_scan)
				{
					count_2d[0] = 1;
                    if(mda_scalers->scan->last_point == 0)
                        count_2d[1] = 1;
                    else
                        count_2d[1] = mda_scalers->scan->last_point;
				}
				else
				{
                    if(mda_scalers->scan->last_point == 0)
                        count_2d[0] = 1;
                    else
                        count_2d[0] = mda_scalers->scan->last_point;
                    if(mda_scalers->scan->sub_scans[0]->last_point == 0)
                        count_2d[1] = 1;
                    else
                        count_2d[1] = mda_scalers->scan->sub_scans[0]->last_point;
				}
				count_3d[1] = count_2d[0];
				count_3d[2] = count_2d[1];

				offset_3d[1] = 0;
				offset_3d[2] = 0;

				memoryspace_id = H5Screate_simple(2, count_2d, NULL);

				scaler_mat.resize(count_2d[0], count_2d[1]);
                scaler_mat.setZero(count_2d[0], count_2d[1]);

                for (auto &itr : scalers)
                {

                    offset[0] = itr.hdf_idx;
                    char tmp_char[255] = {0};
                    itr.hdf_name.copy(tmp_char, 254);
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

                    if(itr.mda_idx < 0)
                    {
                        continue;
                    }

                    if (single_row_scan)
                    {
                        for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
                        {
							real_t val = mda_scalers->scan->detectors_data[itr.mda_idx][i];
                            if(itr.normalize_by_time && mda_time_scaler_idx > -1)
                            {
                                real_t scaler_time_normalizer = 1.0;
                                real_t det_time = mda_scalers->scan->detectors_data[mda_time_scaler_idx][i];
                                scaler_time_normalizer = det_time / time_scaler_clock;
                                val /= scaler_time_normalizer;
                            }
							scaler_mat(0, i) = val;
                        }
                    }
                    else
                    {
                        for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
                        {
                            for (int32_t j = 0; j < mda_scalers->scan->sub_scans[0]->last_point; j++)
                            {
								real_t val = mda_scalers->scan->sub_scans[i]->detectors_data[itr.mda_idx][j];
                                if(itr.normalize_by_time && mda_time_scaler_idx > -1)
                                {
                                    real_t scaler_time_normalizer = 1.0;
                                    real_t det_time = mda_scalers->scan->sub_scans[i]->detectors_data[mda_time_scaler_idx][j];
                                    scaler_time_normalizer = det_time / time_scaler_clock;
									val /= scaler_time_normalizer;
                                }
								scaler_mat(i, j) = val;
                            }
                        }
                    }

                    offset_3d[0] = itr.hdf_idx;
                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)scaler_mat.data());

                }

                if (save_cfg_abs)
                {
					//save calculated names
                    std::string tmp_name = "abs_ic";
                    offset[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    char tmp_char[255] = {0};
                    tmp_name.copy(tmp_char, 254);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

                    tmp_name = "abs_cfg";
                    tmp_name.copy(tmp_char, 254);
                    offset[0] = hdf_idx + 1;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

                    tmp_name = "H_dpc_cfg";
                    tmp_name.copy(tmp_char, 254);
                    offset[0] = hdf_idx + 2;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

                    tmp_name = "V_dpc_cfg";
                    tmp_name.copy(tmp_char, 254);
                    offset[0] = hdf_idx + 3;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

                    tmp_name = "dia1_dpc_cfg";
                    tmp_name.copy(tmp_char, 254);
                    offset[0] = hdf_idx + 4;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

                    tmp_name = "dia2_dpc_cfg";
                    tmp_name.copy(tmp_char, 254);
                    offset[0] = hdf_idx + 5;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

					abs_cfg_mat.resize(count_2d[0], count_2d[1]);
					H_dpc_cfg_mat.resize(count_2d[0], count_2d[1]);
					V_dpc_cfg_mat.resize(count_2d[0], count_2d[1]);
					dia1_dpc_cfg_mat.resize(count_2d[0], count_2d[1]);
					dia2_dpc_cfg_mat.resize(count_2d[0], count_2d[1]);

                    abs_cfg_mat.setZero(count_2d[0], count_2d[1]);
                    H_dpc_cfg_mat.setZero(count_2d[0], count_2d[1]);
                    V_dpc_cfg_mat.setZero(count_2d[0], count_2d[1]);
                    dia1_dpc_cfg_mat.setZero(count_2d[0], count_2d[1]);
                    dia2_dpc_cfg_mat.setZero(count_2d[0], count_2d[1]);

                    if (single_row_scan)
                    {

                        for (int32_t j = 0; j < mda_scalers->scan->last_point; j++)
                        {
                            real_t us_ic = mda_scalers->scan->detectors_data[us_ic_idx][j];
                            real_t ds_ic = mda_scalers->scan->detectors_data[ds_ic_idx][j];
                            real_t t_2 = mda_scalers->scan->detectors_data[cfg_2_idx][j];
                            real_t t_3 = mda_scalers->scan->detectors_data[cfg_3_idx][j];
                            real_t t_4 = mda_scalers->scan->detectors_data[cfg_4_idx][j];
                            real_t t_5 = mda_scalers->scan->detectors_data[cfg_5_idx][j];

                            real_t t_abs = t_2 + t_3 + t_4 + t_5;
                            scaler_mat(0,j) = ds_ic / us_ic;
                            abs_cfg_mat(0,j) = t_abs / us_ic;
                            H_dpc_cfg_mat(0, j) = (t_2 - t_3 - t_4 + t_5) / t_abs;
                            V_dpc_cfg_mat(0, j) = (t_2 + t_3 - t_4 - t_5) / t_abs;
                            dia1_dpc_cfg_mat(0, j) = (t_2 - t_4) / t_abs;
                            dia2_dpc_cfg_mat(0, j) = (t_3 - t_5) / t_abs;
							/*
							if (itr.normalize_by_time)
							{
								real_t scaler_time_normalizer = 1.0;
								if (mda_time_scaler_idx > -1)
								{
                                    real_t det_time = mda_scalers->scan->sub_scans[i]->detectors_data[mda_time_scaler_idx][i];
									scaler_time_normalizer = det_time / time_scaler_clock;
								}
								val /= scaler_time_normalizer;
							}
                            */
                        }
                    }
                    else
                    {
                        for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
                        {
                            for (int32_t j = 0; j < mda_scalers->scan->sub_scans[0]->last_point; j++)
                            {
                                real_t us_ic = mda_scalers->scan->sub_scans[i]->detectors_data[us_ic_idx][j];
                                real_t ds_ic = mda_scalers->scan->sub_scans[i]->detectors_data[ds_ic_idx][j];
                                real_t t_2 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_2_idx][j];
                                real_t t_3 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_3_idx][j];
                                real_t t_4 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_4_idx][j];
                                real_t t_5 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_5_idx][j];

                                real_t t_abs = t_2 + t_3 + t_4 + t_5;
								scaler_mat(i, j) = ds_ic / us_ic;
								abs_cfg_mat(i, j) = t_abs / us_ic;
								H_dpc_cfg_mat(i, j) = (t_2 - t_3 - t_4 + t_5) / t_abs;
								V_dpc_cfg_mat(i, j) = (t_2 + t_3 - t_4 - t_5) / t_abs;
								dia1_dpc_cfg_mat(i, j) = (t_2 - t_4) / t_abs;
								dia2_dpc_cfg_mat(i, j) = (t_3 - t_5) / t_abs;
                            }
                        }
                    }
					offset_3d[0] = hdf_idx;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)scaler_mat.data());
					offset_3d[0] = hdf_idx + 1;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)abs_cfg_mat.data());
					offset_3d[0] = hdf_idx + 2;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)H_dpc_cfg_mat.data());
					offset_3d[0] = hdf_idx + 3;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)V_dpc_cfg_mat.data());
					offset_3d[0] = hdf_idx + 4;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)dia1_dpc_cfg_mat.data());
					offset_3d[0] = hdf_idx + 5;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)dia2_dpc_cfg_mat.data());
                }
            }
            else
            {
                //params_override->scaler_pvs
            }
        }

        H5Tclose(filetype);
        H5Tclose(memtype);

        if(dset_names_id > -1)
        {
            H5Dclose(dset_names_id);
            dset_names_id = -1;
        }
        if(dset_cps_id > -1)
        {
            H5Dclose(dset_cps_id);
            dset_cps_id = -1;
        }
        if(dataspace_id > -1)
        {
            H5Sclose(dataspace_id);
            dataspace_id = -1;
        }
        if(filespace_id > -1)
        {
            H5Sclose(filespace_id);
            filespace_id = -1;
        }
        if(filespace_name_id > -1)
        {
            H5Sclose(filespace_name_id);
            filespace_name_id = -1;
        }
		if (memoryspace_str_id > -1)
		{
			H5Sclose(memoryspace_str_id);
			memoryspace_str_id = -1;
		}
        if(memoryspace_id > -1)
        {
            H5Sclose(memoryspace_id);
            memoryspace_id = -1;
        }
        if(dcpl_id > -1)
        {
            H5Pclose(dcpl_id);
            dcpl_id = -1;
        }
        if(scalers_grp_id > -1)
        {
            H5Gclose(scalers_grp_id);
            scalers_grp_id = -1;
        }
    }
    catch (...)
    {
		
		if (memoryspace_str_id > -1)
			H5Sclose(memoryspace_str_id);
        if(memoryspace_id > -1)
            H5Sclose(memoryspace_id);
        if(dset_names_id > -1)
            H5Dclose(dset_names_id);
        if(dset_cps_id > -1)
            H5Dclose(dset_cps_id);
        if(dataspace_id > -1)
            H5Sclose(dataspace_id);
        if(filespace_id > -1)
            H5Sclose(filespace_id);
        if(filespace_name_id > -1)
            H5Sclose(filespace_name_id);
        if(dcpl_id > -1)
            H5Pclose(dcpl_id);
        if(scalers_grp_id > -1)
            H5Gclose(scalers_grp_id);

        logit << "Error creating group MAPS/Scalers" << std::endl;
        return false;
    }
    logit << "Done" << std::endl;
    return true;
}

//-----------------------------------------------------------------------------

void HDF5_IO::_save_amps(hid_t scalers_grp_id, struct mda_file *mda_scalers, data_struct::xrf::Params_Override * params_override)
{
    hid_t dataspace_id = -1, memoryspace_id = -1;
    hid_t dset_id = -1;
    hid_t status;
    hid_t filetype, memtype;

    logit << std::endl;

    char tmp_char[255] = {0};
    hsize_t offset[1] = { 0 };
    hsize_t count[1] = { 3 };
    io::file::MDA_IO mda_io;

    filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, 256);
    memtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(memtype, 255);

    real_t us_amp_sens_num_val = params_override->us_amp_sens_num;
    mda_io.find_scaler_index(mda_scalers, params_override->us_amp_sens_num_pv, us_amp_sens_num_val);
    real_t us_amp_sens_unit_val = params_override->us_amp_sens_unit;
    mda_io.find_scaler_index(mda_scalers, params_override->us_amp_sens_unit_pv, us_amp_sens_unit_val);

    real_t ds_amp_sens_num_val = params_override->ds_amp_sens_num;
    mda_io.find_scaler_index(mda_scalers, params_override->ds_amp_sens_num_pv, ds_amp_sens_num_val);
    real_t ds_amp_sens_unit_val = params_override->ds_amp_sens_unit;
    mda_io.find_scaler_index(mda_scalers, params_override->ds_amp_sens_unit_pv, ds_amp_sens_unit_val);

    real_t trans_us_amp_sens_num_val;
    std::string trans_us_amp_sens_unit;
    real_t trans_ds_amp_sens_num_val;
    std::string trans_ds_amp_sens_unit;

    switch((int)us_amp_sens_num_val)
    {
        case 0:
            trans_us_amp_sens_num_val = 1;
            break;
        case 1:
            trans_us_amp_sens_num_val = 2;
            break;
        case 2:
            trans_us_amp_sens_num_val = 5;
            break;
        case 3:
            trans_us_amp_sens_num_val = 10;
            break;
        case 4:
            trans_us_amp_sens_num_val = 20;
            break;
        case 5:
            trans_us_amp_sens_num_val = 50;
            break;
        case 6:
            trans_us_amp_sens_num_val = 100;
            break;
        case 7:
            trans_us_amp_sens_num_val = 200;
            break;
        case 8:
            trans_us_amp_sens_num_val = 500;
            break;
    }
    switch((int)ds_amp_sens_num_val)
    {
        case 0:
            trans_ds_amp_sens_num_val = 1;
            break;
        case 1:
            trans_ds_amp_sens_num_val = 2;
            break;
        case 2:
            trans_ds_amp_sens_num_val = 5;
            break;
        case 3:
            trans_ds_amp_sens_num_val = 10;
            break;
        case 4:
            trans_ds_amp_sens_num_val = 20;
            break;
        case 5:
            trans_ds_amp_sens_num_val = 50;
            break;
        case 6:
            trans_ds_amp_sens_num_val = 100;
            break;
        case 7:
            trans_ds_amp_sens_num_val = 200;
            break;
        case 8:
            trans_ds_amp_sens_num_val = 500;
            break;
    }

    switch((int)us_amp_sens_unit_val)
    {
        case 0:
            trans_us_amp_sens_unit= "pA/V";
            break;
        case 1:
            trans_us_amp_sens_unit= "nA/V";
            break;
        case 2:
            trans_us_amp_sens_unit= "uA/V";
            break;
        case 3:
            trans_us_amp_sens_unit= "mA/V";
            break;
    }

    switch((int)ds_amp_sens_unit_val)
    {
        case 0:
            trans_ds_amp_sens_unit= "pA/V";
            break;
        case 1:
            trans_ds_amp_sens_unit= "nA/V";
            break;
        case 2:
            trans_ds_amp_sens_unit= "uA/V";
            break;
        case 3:
            trans_ds_amp_sens_unit= "mA/V";
            break;
    }

    dataspace_id = H5Screate_simple(1, count, NULL);

    count[0] = 1;
    memoryspace_id = H5Screate_simple(1, count, NULL);

    dset_id = H5Dcreate(scalers_grp_id, "us_amp", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&us_amp_sens_num_val);

    offset[0] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&us_amp_sens_unit_val);

    offset[0] = 2;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&trans_us_amp_sens_num_val);

    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }


    offset[0] = 0;
    dset_id = H5Dcreate(scalers_grp_id, "ds_amp", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&ds_amp_sens_num_val);

    offset[0] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&ds_amp_sens_unit_val);

    offset[0] = 2;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&trans_ds_amp_sens_num_val);
    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }
    if(dataspace_id > -1)
    {
        H5Sclose(dataspace_id);
        dataspace_id = -1;
    }

    offset[0] = 0;
    count[0] = 1;
    dataspace_id = H5Screate_simple(1, count, NULL);

    dset_id = H5Dcreate(scalers_grp_id, "us_amp_num", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&trans_us_amp_sens_num_val);
    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }

    dset_id = H5Dcreate(scalers_grp_id, "ds_amp_num", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&trans_ds_amp_sens_num_val);
    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }


    trans_us_amp_sens_unit.copy(tmp_char, 254);
    dset_id = H5Dcreate(scalers_grp_id, "us_amp_unit", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);
    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }

    trans_ds_amp_sens_unit.copy(tmp_char, 254);
    dset_id = H5Dcreate(scalers_grp_id, "ds_amp_unit", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);


    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }
    if(dataspace_id > -1)
    {
        H5Sclose(dataspace_id);
        dataspace_id = -1;
    }
    if(memoryspace_id > -1)
    {
        H5Sclose(memoryspace_id);
        memoryspace_id = -1;
    }


    H5Tclose(filetype);
    H5Tclose(memtype);
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_scan_scalers(size_t detector_num,
                                struct mda_file *mda_scalers,
                                data_struct::xrf::Params_Override * params_override,
                                bool hasNetcdf,
                                size_t row_idx_start,
                                int row_idx_end,
                                size_t col_idx_start,
                                int col_idx_end)
{

    std::lock_guard<std::mutex> lock(_mutex);
	std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    hid_t scan_grp_id, maps_grp_id;

	if (mda_scalers == nullptr)
    {
        logit << "Warning: mda_scalers == nullptr. Not returning from save_scan_scalers" << std::endl;
		return false;
    }

    if(_cur_file_id < 0)
    {
        logit << "Error: hdf5 file was never initialized. Call start_save_seq() before this function." << std::endl;
        return false;
    }

    logit << "Saving scalers to hdf5"<< std::endl;

    maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
    if (maps_grp_id < 0)
        maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (maps_grp_id < 0)
    {
        logit << "Error creating group 'MAPS'" << std::endl;
        return false;
    }

    scan_grp_id = H5Gopen(maps_grp_id, "Scan", H5P_DEFAULT);
    if (scan_grp_id < 0)
        scan_grp_id = H5Gcreate(maps_grp_id, "Scan", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (scan_grp_id < 0)
    {
        logit << "Error creating group MAPS/Scan" << std::endl;
        return false;
    }

    _save_scan_meta_data(scan_grp_id, mda_scalers);
	
    _save_extras(scan_grp_id, mda_scalers);
	
    _save_scalers(maps_grp_id, mda_scalers, detector_num, params_override, hasNetcdf);

    H5Gclose(scan_grp_id);
    H5Gclose(maps_grp_id);

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	//std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_scan_scalers_confocal(std::string path,
                                size_t detector_num,
                                size_t row_idx_start,
                                int row_idx_end,
                                size_t col_idx_start,
                                int col_idx_end)
{

    std::lock_guard<std::mutex> lock(_mutex);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    hid_t scan_grp_id, maps_grp_id, scalers_grp_id, status, error, ocpypl_id;
    hid_t    file_id, src_maps_grp_id;
    hid_t    dataspace_detectors_id, dset_detectors_id, attr_detector_names_id;
    hid_t   xpos_dataspace_id, xpos_id, ypos_dataspace_id, ypos_id;
    hid_t dataspace_id, dataset_id;
    char* detector_names[256];
    int det_rank;
    hsize_t* det_dims_in;
    hsize_t n_offset[1] = {0};
    hsize_t n_count[1] = {1};
    hsize_t x_offset[3] = {0,0,0};
    hsize_t x_count[3] = {1,1,1};
    hsize_t y_offset[2] = {0,0};
    hsize_t y_count[2] = {1,1};
    void *f_data;

    if(_cur_file_id < 0)
    {
        logit << "Error: hdf5 file was never initialized. Call start_save_seq() before this function." << std::endl;
        return false;
    }

    logit << "Saving scalers to hdf5"<< std::endl;

    maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
    if (maps_grp_id < 0)
        maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (maps_grp_id < 0)
    {
        logit << "Error creating group 'MAPS'" << std::endl;
        return false;
    }

    scan_grp_id = H5Gopen(maps_grp_id, "Scan", H5P_DEFAULT);
    if (scan_grp_id < 0)
        scan_grp_id = H5Gcreate(maps_grp_id, "Scan", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (scan_grp_id < 0)
    {
        H5Gclose(maps_grp_id);
        logit << "Error creating group MAPS/Scan" << std::endl;
        return false;
    }

    scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
    if (scalers_grp_id < 0)
        scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (scalers_grp_id < 0)
    {
        H5Gclose(scan_grp_id);
        H5Gclose(scalers_grp_id);
        H5Gclose(maps_grp_id);
        logit << "Error creating group MAPS/Scalers" << std::endl;
        return false;
    }

    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
        return false;

    if ( false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "2D Scan", file_id) )
       return false;

    //_save_scan_meta_data(scan_grp_id, mda_scalers);
    if ( false == _open_h5_object(xpos_id, H5O_DATASET, close_map, "X Positions", src_maps_grp_id) )
       return false;
    xpos_dataspace_id = H5Dget_space(xpos_id);
    close_map.push({xpos_dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(ypos_id, H5O_DATASET, close_map, "Y Positions", src_maps_grp_id) )
       return false;
    ypos_dataspace_id = H5Dget_space(ypos_id);
    close_map.push({ypos_dataspace_id, H5O_DATASPACE});

    hid_t xy_type = H5Dget_type(xpos_id);
    det_rank = H5Sget_simple_extent_ndims(xpos_dataspace_id);
    det_dims_in = new hsize_t[det_rank];
    H5Sget_simple_extent_dims(xpos_dataspace_id, &det_dims_in[0], NULL);

    dataspace_id = H5Screate_simple (1, &det_dims_in[1], NULL);
    dataset_id = H5Dcreate(scan_grp_id, "x_axis", xy_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    x_count[1] = det_dims_in[1];
    f_data = malloc(sizeof(float) * det_dims_in[1]);
    H5Sselect_hyperslab (xpos_dataspace_id, H5S_SELECT_SET, x_offset, NULL, x_count, NULL);
    status = H5Dread (xpos_id, xy_type, dataspace_id, xpos_dataspace_id, H5P_DEFAULT, f_data);
    status = H5Dwrite (dataset_id, xy_type, H5S_ALL, dataspace_id, H5P_DEFAULT, f_data);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    free(f_data);


    dataspace_id = H5Screate_simple (1, &det_dims_in[0], NULL);
    dataset_id = H5Dcreate(scan_grp_id, "y_axis", xy_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    y_count[1] = det_dims_in[0];
    f_data = malloc(sizeof(float) * det_dims_in[0]);
    H5Sselect_hyperslab (ypos_dataspace_id, H5S_SELECT_SET, y_offset, NULL, y_count, NULL);
    status = H5Dread (ypos_id, xy_type, dataspace_id, ypos_dataspace_id, H5P_DEFAULT, f_data);
    status = H5Dwrite (dataset_id, xy_type, H5S_ALL, dataspace_id, H5P_DEFAULT, f_data);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    free(f_data);
    delete [] det_dims_in;

    //Save scalers
    ocpypl_id = H5Pcreate(H5P_OBJECT_COPY);
    status = H5Ocopy(src_maps_grp_id, "Detectors", scalers_grp_id, "Values", ocpypl_id, H5P_DEFAULT);

    //save the scalers names
    if ( false == _open_h5_object(dset_detectors_id, H5O_DATASET, close_map, "Detectors", src_maps_grp_id) )
       return false;
    dataspace_detectors_id = H5Dget_space(dset_detectors_id);
    close_map.push({dataspace_detectors_id, H5O_DATASPACE});
    attr_detector_names_id=H5Aopen(dset_detectors_id, "Detector Names", H5P_DEFAULT);
    close_map.push({dataspace_detectors_id, H5O_ATTRIBUTE});

    det_rank = H5Sget_simple_extent_ndims(dataspace_detectors_id);
    det_dims_in = new hsize_t[det_rank];
    H5Sget_simple_extent_dims(dataspace_detectors_id, &det_dims_in[0], NULL);

    hid_t ftype = H5Aget_type(attr_detector_names_id);
    hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
    error = H5Aread(attr_detector_names_id, type, detector_names);

    if(error == 0)
    {
        dataspace_id = H5Screate_simple (1, &det_dims_in[2], NULL);
        hid_t dataset_id = H5Dcreate(scalers_grp_id, "Names", type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, detector_names);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        for(int z = 0; z < det_dims_in[2]; z++)
        {
            free(detector_names[z]);
        }
    }
    delete [] det_dims_in;


    _close_h5_objects(close_map);

    H5Gclose(scan_grp_id);
    H5Gclose(scalers_grp_id);
    H5Gclose(maps_grp_id);


    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::generate_avg(std::string avg_filename, std::vector<std::string> files_to_avg)
{
    std::lock_guard<std::mutex> lock(_mutex);
    logit  << avg_filename << std::endl;

    hid_t ocpypl_id, status, src_maps_grp_id, src_analyzed_grp_id, dst_fit_grp_id;
    std::vector<hid_t> hdf5_file_ids;
    std::string group_name = "";

    //save to restore later
    hid_t saved_file_id = _cur_file_id;

    for(auto& filename : files_to_avg)
    {
        hid_t    rd_file_id;

        rd_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (rd_file_id > -1)
        {
            hdf5_file_ids.push_back(rd_file_id);
        }
    }

    if(hdf5_file_ids.size() > 1)
    {
        hid_t file_id = H5Fcreate(avg_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        hid_t dst_maps_grp_id = H5Gcreate(file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        ocpypl_id = H5Pcreate(H5P_OBJECT_COPY);
        status = H5Pset_copy_object(ocpypl_id, H5O_COPY_MERGE_COMMITTED_DTYPE_FLAG);

        src_maps_grp_id = H5Gopen(hdf5_file_ids[0], "MAPS", H5P_DEFAULT);
        //Copy Scan and Scalers
        if(src_maps_grp_id > -1)
        {
            status = H5Ocopy(src_maps_grp_id, "Scan", dst_maps_grp_id, "Scan", ocpypl_id, H5P_DEFAULT);
            status = H5Ocopy(src_maps_grp_id, "Scalers", dst_maps_grp_id, "Scalers", ocpypl_id, H5P_DEFAULT);
        }

        //Average XRF_Analyzed
        _generate_avg_analysis(src_maps_grp_id, dst_maps_grp_id, "MAPS", ocpypl_id, hdf5_file_ids);

        //Average Spectra
        group_name = "Spectra";
        src_analyzed_grp_id = H5Gopen(src_maps_grp_id, group_name.c_str(), H5P_DEFAULT);
        if(src_analyzed_grp_id > -1)
        {
            dst_fit_grp_id = H5Gcreate(dst_maps_grp_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            _gen_average("/MAPS/"+group_name+"/Elapsed_Livetime", "Elapsed_Livetime", src_analyzed_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
            _gen_average("/MAPS/"+group_name+"/Elapsed_Realtime", "Elapsed_Realtime", src_analyzed_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
            _gen_average("/MAPS/"+group_name+"/Input_Counts", "Input_Counts", src_analyzed_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
            _gen_average("/MAPS/"+group_name+"/Output_Counts", "Output_Counts", src_analyzed_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);

            //for now copy the first detector energy and energy_calib, in the future would be nice to normalize over 4 detectors
            status = H5Ocopy(src_analyzed_grp_id, "Energy", dst_fit_grp_id, "Energy", ocpypl_id, H5P_DEFAULT);
            status = H5Ocopy(src_analyzed_grp_id, "Energy_Calibration", dst_fit_grp_id, "Energy_Calibration", ocpypl_id, H5P_DEFAULT);
            //sum mca_arr
            _gen_average("/MAPS/"+group_name+"/mca_arr", "mca_arr", src_analyzed_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids, false);


            _generate_avg_integrated_spectra(src_analyzed_grp_id, dst_fit_grp_id, "MAPS/Spectra", ocpypl_id, hdf5_file_ids);

            H5Gclose(dst_fit_grp_id);
            H5Gclose(src_analyzed_grp_id);
        }

        //Average Quantification
        src_analyzed_grp_id = H5Gopen(src_maps_grp_id, "Quantification", H5P_DEFAULT);
        if(src_analyzed_grp_id > -1)
        {
            dst_fit_grp_id = H5Gcreate(dst_maps_grp_id, "Quantification", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            _generate_avg_analysis(src_analyzed_grp_id, dst_fit_grp_id, "MAPS/Quantification", ocpypl_id, hdf5_file_ids);

            status = H5Ocopy(src_analyzed_grp_id, "Element_Weights", dst_fit_grp_id, "Element_Weights", ocpypl_id, H5P_DEFAULT);
            status = H5Ocopy(src_analyzed_grp_id, "Element_Weights_Names", dst_fit_grp_id, "Element_Weights_Names", ocpypl_id, H5P_DEFAULT);
            status = H5Ocopy(src_analyzed_grp_id, "Scalers", dst_fit_grp_id, "Scalers", ocpypl_id, H5P_DEFAULT);
            status = H5Ocopy(src_analyzed_grp_id, "Standard_Name", dst_fit_grp_id, "Standard_Name", ocpypl_id, H5P_DEFAULT);


            _generate_avg_integrated_spectra(src_analyzed_grp_id, dst_fit_grp_id, "MAPS/Quantification", ocpypl_id, hdf5_file_ids);

            H5Gclose(dst_fit_grp_id);
            H5Gclose(src_analyzed_grp_id);
        }

        status = H5Ocopy(src_maps_grp_id, "version", dst_maps_grp_id, "version", ocpypl_id, H5P_DEFAULT);

        if(src_maps_grp_id > -1)
            H5Gclose(src_maps_grp_id);

        H5Gclose(dst_maps_grp_id);
        //H5Fflush(file_id, H5F_SCOPE_LOCAL);
        //H5Fclose(file_id);
        _cur_file_id = file_id;
        end_save_seq();
    }
    else
    {
        logit<<"Error: Could not file any files to average for dataset "<<avg_filename<<std::endl;
    }

    for(auto& f_id : hdf5_file_ids)
    {
        _cur_file_id = f_id;
        end_save_seq();
        //H5Fclose(f_id);
    }

    logit<<"closing file"<<std::endl;

    _cur_file_id = saved_file_id;

    return true;
}

//-----------------------------------------------------------------------------

void HDF5_IO::_gen_average(std::string full_hdf5_path, std::string dataset_name, hid_t src_fit_grp_id, hid_t dst_fit_grp_id, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids, bool avg)
{
    std::vector<hid_t> analysis_ids;
    hid_t error, status;

//    status = H5Ocopy(src_fit_grp_id, dataset_name.c_str(), dst_fit_grp_id, dataset_name.c_str(), ocpypl_id, H5P_DEFAULT);

//    if (status == 0)
    {
        hid_t dset_id = H5Dopen2(src_fit_grp_id, dataset_name.c_str(), H5P_DEFAULT);
        hid_t dataspace_id = H5Dget_space(dset_id);
        hid_t file_type = H5Dget_type(dset_id);
        hid_t props = H5Dget_create_plist(dset_id);
        hid_t dst_dset_id = H5Dcreate(dst_fit_grp_id, dataset_name.c_str(), file_type, dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);
        if(dst_dset_id < 1)
        {
            logit<<"Error "<<full_hdf5_path<< " " <<dataset_name<<std::endl;
            return;
        }
        int rank = H5Sget_simple_extent_ndims(dataspace_id);

        hsize_t* dims_in = new hsize_t[rank];
        unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);

        unsigned long total = 1;
        for(int i=0; i< rank; i++)
        {
            total *= dims_in[i];
        }
        //get all the other files dataset ids
        for(int k=1; k<hdf5_file_ids.size(); k++)
        {
            hid_t det_analysis_dset_id = H5Dopen2(hdf5_file_ids[k], full_hdf5_path.c_str(), H5P_DEFAULT);
            if( det_analysis_dset_id > -1 )
                analysis_ids.push_back(det_analysis_dset_id);
        }

        std::valarray<real_t> buffer1(total);
        std::valarray<real_t> buffer2(total);
        float divisor = 1.0;
        error = H5Dread(dset_id, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer1[0]);
        for(int k=0; k<analysis_ids.size(); k++)
        {
            error = H5Dread(analysis_ids[k], H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer2[0]);
            if(error > -1)
            {
                buffer1 += buffer2;
                divisor += 1.0;
            }
            else
            {
            logit<<"Error reading "<<full_hdf5_path<<" dataset "<<std::endl;
            }
        }

        if(avg)
        {
            buffer1 /= divisor;
        }

        //hid_t dst_dset_id = H5Dopen2(dst_fit_grp_id, dataset_name.c_str(), H5P_DEFAULT);
        error = H5Dwrite(dst_dset_id, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer1[0]);

        for(int k=0; k<analysis_ids.size(); k++)
        {
            H5Dclose(analysis_ids[k]);
        }

        //clean up
        delete [] dims_in;
        H5Dclose(dset_id);
        H5Dclose(dst_dset_id);
    }
}

//-----------------------------------------------------------------------------

void HDF5_IO::_generate_avg_analysis(hid_t src_maps_grp_id, hid_t dst_maps_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids)
{
    hid_t src_analyzed_grp_id = H5Gopen(src_maps_grp_id, "XRF_Analyzed", H5P_DEFAULT);
    if(src_analyzed_grp_id > -1)
    {
        hid_t dst_analyzed_grp_id = H5Gcreate(dst_maps_grp_id, "XRF_Analyzed", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for(std::string analysis_grp_name : xrf_analysis_save_names)
        {
            hid_t src_fit_grp_id = H5Gopen(src_analyzed_grp_id, analysis_grp_name.c_str(), H5P_DEFAULT);
            if(src_fit_grp_id > -1)
            {
                //std::vector<hid_t> analysis_ids;

                hid_t dst_fit_grp_id = H5Gcreate(dst_analyzed_grp_id, analysis_grp_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                //copy channel names
                hid_t status = H5Ocopy(src_fit_grp_id, "Channel_Names", dst_fit_grp_id, "Channel_Names", ocpypl_id, H5P_DEFAULT);
                status = H5Ocopy(src_fit_grp_id, "Calibration_Curve_Labels", dst_fit_grp_id, "Calibration_Curve_Labels", ocpypl_id, H5P_DEFAULT);
                _gen_average(group_name+"/XRF_Analyzed/"+analysis_grp_name+"/Counts_Per_Sec", "Counts_Per_Sec", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);

                _gen_average(group_name+"/XRF_Analyzed/"+analysis_grp_name+"/Calibration_Curve_Current", "Calibration_Curve_Current", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                _gen_average(group_name+"/XRF_Analyzed/"+analysis_grp_name+"/Calibration_Curve_DS_IC", "Calibration_Curve_DS_IC", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                _gen_average(group_name+"/XRF_Analyzed/"+analysis_grp_name+"/Calibration_Curve_US_IC", "Calibration_Curve_US_IC", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);

                H5Gclose(dst_fit_grp_id);
                H5Gclose(src_fit_grp_id);

            }
        }
        H5Gclose(dst_analyzed_grp_id);
        H5Gclose(src_analyzed_grp_id);
    }
}

//-----------------------------------------------------------------------------

void HDF5_IO::_generate_avg_integrated_spectra(hid_t src_analyzed_grp_id, hid_t dst_fit_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids)
{
    hid_t src_inner_grp_id = H5Gopen(src_analyzed_grp_id, "Integrated_Spectra", H5P_DEFAULT);
    if(src_inner_grp_id > -1)
    {
        hid_t dst_inner_grp_id = H5Gcreate(dst_fit_grp_id, "Integrated_Spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        _gen_average(group_name+"/Integrated_Spectra/Elapsed_Livetime", "Elapsed_Livetime", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);
        _gen_average(group_name+"/Integrated_Spectra/Elapsed_Realtime", "Elapsed_Realtime", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);
        _gen_average(group_name+"/Integrated_Spectra/Input_Counts", "Input_Counts", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);
        _gen_average(group_name+"/Integrated_Spectra/Output_Counts", "Output_Counts", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);

        //don't average the integrated spectra, just sum it
        _gen_average(group_name+"/Integrated_Spectra/Spectra", "Spectra", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids, false);

        H5Gclose(dst_inner_grp_id);
        H5Gclose(src_inner_grp_id);
    }
}

//-----------------------------------------------------------------------------

bool HDF5_IO::generate_stream_dataset(std::string dataset_directory,
                                      std::string dataset_name,
                                      int detector_num,
                                      size_t height,
                                      size_t width)
{

    //    std::string str_detector_num = std::to_string(detector_num);
    //    std::string full_save_path = dataset_directory+"/img.dat/"+dataset_file+".h5"+str_detector_num;
    //    io::file::HDF5_IO::inst()->set_filename(full_save_path);


    return false;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_stream_row(size_t d_hash,
                             size_t detector_num,
                             size_t row,
                             std::vector< data_struct::xrf::Spectra* >  *spectra_row)
{
    return false;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_itegrade_spectra(data_struct::xrf::Spectra * spectra)
{
    return false;
}

//-----------------------------------------------------------------------------


bool HDF5_IO::close_dataset(size_t d_hash)
{
    return false;
}

//-----------------------------------------------------------------------------

} //end namespace file
}// end namespace io
