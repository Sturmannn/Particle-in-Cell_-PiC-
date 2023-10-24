#ifndef __FIELD_HPP__
#define __FIELD_HPP__

#include <vector>

namespace Field {

  class ComputingField {
  public:
    ComputingField() = delete;
    ComputingField(const uint64_t _Nx, const uint64_t _Ny);
    ComputingField(const ComputingField& _field);
    ~ComputingField() = default;

    double& operator()(const uint64_t i, const uint64_t j);
    ComputingField& operator=(const ComputingField& _field);

  private:
    uint64_t Nx, Ny;  // count of cells
    std::vector<double> field;
  };

}  // namespace Field
#endif  // !__FIELD_HPP__