#ifndef __FIELD_HPP__
#define __FIELD_HPP__

#include <vector>
#include <fstream>
#include <iostream>

constexpr char* path_to_analytic_data = "..\\..\\input_for_graphs\\analytical_data.csv"; // Writing E, B to a file
constexpr char* path_to_calculated_data = "..\\..\\input_for_graphs\\my_data.csv"; // Writing E, B to a file

namespace Field {

  class ComputingField {
  public:
    ComputingField() = delete;
    ComputingField(const uint64_t _Nx, const uint64_t _Ny);
    ComputingField(const ComputingField& _field);
    ~ComputingField() = default;

    double& operator()(uint64_t i, uint64_t j);
    //double& operator()(const int64_t i, const int64_t j);
    ComputingField& operator=(const ComputingField& _field);

    void write_to_file(const double _dx, const uint64_t j = 0ull);

  private:
    uint64_t Nx, Ny;  // count of cells
    std::vector<double> field;
  };

}  // namespace Field
#endif  // !__FIELD_HPP__