#ifndef __FIELD_HPP__
#define __FIELD_HPP__

#include <vector>
#include <fstream>
#include <iostream>

constexpr char* path_to_analytic_data = "..\\..\\input_for_graphs\\analytical_data.csv"; // Writing E, B to a file
constexpr char* path_to_calculated_data = "..\\..\\input_for_graphs\\my_data.csv";

namespace Field {

  class ComputingField {
  public:
    ComputingField() = delete;
    ComputingField(const uint64_t _Nx, const uint64_t _Ny);
    ComputingField(const ComputingField& _field);
    ComputingField(ComputingField&& _field) noexcept;
    ~ComputingField() = default;

    uint64_t get_Nx() const noexcept { return Nx; }
    uint64_t get_Ny() const noexcept { return Ny; }

    double& operator()(uint64_t i, uint64_t j);
    const double& operator()(uint64_t i, uint64_t j) const;
    //double& operator()(const int64_t i, const int64_t j);
    ComputingField& operator=(const ComputingField& _field);
    ComputingField& operator=(ComputingField&& _field) noexcept;

    void write_field_to_file_OX(const char* path, const uint64_t j = 0ull);
    void write_field_to_file_OY(const char* path, const uint64_t i = 0ull);
    static void clear_file(const char* path);

  private:
    uint64_t Nx, Ny;  // count of cells
    std::vector<double> field;
  };

}  // namespace Field
#endif  // !__FIELD_HPP__