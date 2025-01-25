#ifndef _SOLUTION_DATA_H_
#define _SOLUTION_DATA_H_

#include<Utils.h>
#include<Vector3D.h>
#include<map>
#include<vector>
#include<memory>
#include<utility>
#include<cassert>

/**********************************************
 * A data container which stores solutions from
 * exising FSI simulations. The data is mapped
 * based on time stamps for fast access. Methods
 * to simplify the data communication between
 * processors is also included.
 *********************************************/

class SolutionData3D {

  //! container that holds all the data. Should
  //! not copy this container as it takes a lot
  //! of memory.
  std::shared_ptr<std::map<double, std::vector<Vec3D>>> data_ptr;
  
  int num_stamps;
  int num_data_rows;

public:

  SolutionData3D();
  // move constructor
  //SolutionData3D(std::map<double, std::vector<Vec3D>> &&);
  ~SolutionData3D();

  
  bool is_empty() const { return data_ptr->empty(); }
  int get_size() const { return num_stamps + 3*num_data_rows; }
  int get_rows() const { return num_data_rows; }

  std::vector<Vec3D> &find(double time);
  std::vector<Vec3D>  find(double time) const;

  void move_map(std::map<double, std::vector<Vec3D>> &&);
  void insert(double time, std::vector<Vec3D> &&);

  // methods for communication
  void flatten(std::vector<double> &flat); 
  /*SolutionData3D&*/ void rebuild(const std::vector<double> &other, int num_rows);

};

#endif
