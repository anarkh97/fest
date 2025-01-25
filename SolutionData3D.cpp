#include<SolutionData3D.h>

//------------------------------------------------------------

SolutionData3D::SolutionData3D()
{
  data_ptr      = std::make_shared<std::map<double, std::vector<Vec3D>>>();
  num_stamps    = -1;
  num_data_rows = -1;
}

//------------------------------------------------------------
/*
SolutionData3D::SolutionData3D(std::map<double, std::vector<Vec3D>> &&input_map)
{
  data_ptr      = std::make_shared<std::map<double, std::vector<Vec3D>>>(
                    std::move(input_map));
  num_stamps    = data_ptr->size();
  num_data_rows = data_ptr->begin()->second.size();
}
*/

//------------------------------------------------------------

SolutionData3D::~SolutionData3D()
{
  // smart pointer is automatically deleted.
}

//------------------------------------------------------------

//! Takes ownership of the input map. It will be destroyed in the
//! calling function.
void
SolutionData3D::move_map(std::map<double, std::vector<Vec3D>> &&input_map)
{
  data_ptr       = std::make_shared<std::map<double, std::vector<Vec3D>>>(
                     std::move(input_map));
  num_stamps     = data_ptr->size();
  num_data_rows  = data_ptr->begin()->second.size();
}

//------------------------------------------------------------

std::vector<Vec3D>&
SolutionData3D::find(double time)
{

  assert(!is_empty()); // should not be empty
  return data_ptr->at(time); // will throw an error if time is not found.

}

//------------------------------------------------------------

std::vector<Vec3D>
SolutionData3D::find(double time) const
{

  assert(!is_empty()); // should not be empty
  return data_ptr->at(time); // will throw an error if time is not found.

}

//------------------------------------------------------------

void
SolutionData3D::flatten(std::vector<double> &flat)
{

  assert(data_ptr); // should not be null

  int container_size = num_stamps + num_data_rows*3;
  if(flat.empty())
    flat.resize(container_size);
   
  int index = 0;
  for(const auto& [time, data] : *data_ptr) {

    // store time stamp
    flat[index] = time;
    index++;

    // store Vec3D data
    for(const auto& vec : data) {
      flat[index]   = vec[0];
      flat[index+1] = vec[1];
      flat[index+2] = vec[2];

      index += 3;
    }

    assert(index < container_size);

  }

}

//------------------------------------------------------------

//SolutionData3D&
void
SolutionData3D::rebuild(const std::vector<double> &other, int rows)
{

  // assign an empty map
  if(!data_ptr)
    data_ptr = std::make_shared<std::map<double, std::vector<Vec3D>>>();

  int size = other.size();
  
  for(int i=0; i<size; ++i) {

    double time;
    std::vector<Vec3D> data;

    time = other[i];
    i++;

    for(int j=0; j<rows; ++j) {
    
      data[j] = other[i];
      data[j] = other[i+1];
      data[j] = other[i+2];

      i += 3;

    }

    (*data_ptr)[time] = std::move(data);

  }

  num_stamps    = data_ptr->size();
  num_data_rows = data_ptr->begin()->second.size();

  //return *this;

}

//------------------------------------------------------------

