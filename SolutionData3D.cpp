#include<SolutionData3D.h>
#include<iterator>

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
SolutionData3D::MoveMap(std::map<double, std::vector<Vec3D>> &&input_map)
{
  data_ptr       = std::make_shared<std::map<double, std::vector<Vec3D>>>(
                     std::move(input_map));
  num_stamps     = data_ptr->size();
  num_data_rows  = data_ptr->begin()->second.size();
}

//------------------------------------------------------------

std::array<double,2>
SolutionData3D::GetTimeBounds()
{

  double tmin = data_ptr->begin()->first;
  double tmax = data_ptr->rbegin()->first;
  return std::array<double,2>{tmin, tmax};

}

//------------------------------------------------------------

std::array<double,2>
SolutionData3D::GetTimeBracket(double time)
{

  auto upp = data_ptr->lower_bound(time);
  auto low = std::prev(upp);

  double t0 = low->first;
  double t1 = upp->first;

  return std::array<double,2>{t0, t1};

}

//------------------------------------------------------------

void
SolutionData3D::Insert(double time, std::vector<Vec3D> &&data)
{

  (*data_ptr)[time] = std::move(data);
  num_stamps        = data_ptr->size();
  num_data_rows     = data_ptr->begin()->second.size();

}

//------------------------------------------------------------

std::vector<Vec3D>&
SolutionData3D::GetSolutionAtTime(double t)
{
  assert(data_ptr); // should not be null
  return data_ptr->at(t); // will throw an error if not found.
}

//------------------------------------------------------------

void
SolutionData3D::Flatten(std::vector<double> &flat)
{

  assert(data_ptr); // should not be null

  int container_size = GetSize();
  if(flat.empty())
    flat.resize(container_size);
  
  // additional check
  assert((int)flat.size() == container_size);
   
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
SolutionData3D::Rebuild(const std::vector<double> &other, int rows)
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

