#include "tests.hpp"

//void gtest::set_default_field(FDTD::FDTD& _field,
//  std::pair<double, double>& ax_bx) {
//
//  double x = ax_bx.first;
//
//  Field::ComputingField& Ex = _field.get_Ex();
//  Field::ComputingField& Ey = _field.get_Ey();
//  Field::ComputingField& Ez = _field.get_Ez();
//
//  Field::ComputingField& Bx = _field.get_Bx();
//  Field::ComputingField& By = _field.get_By();
//  Field::ComputingField& Bz = _field.get_Bz();
//
//  for (uint64_t i = 0; i < _field.get_Nx(); ++i, x += _field.get_dx())
//    for (uint64_t j = 0; j < _field.get_Ny(); ++j) {
//      Ey(i, j) = Bz(i, j) =
//        sin(2.0 * PI * (x - ax_bx.first) /
//          (ax_bx.second - ax_bx.first));
//      Ex(i, j) = Ez(i, j) = Bx(i, j) = By(i, j) = 0.0;
//    }
//}


void gtest::set_default_field_OY(FDTD::FDTD& _field, const Component E, const Component B)
{
  //if (E == )

  std::pair<double, double>ax_bx = _field.get_ax_bx();
  double x = ax_bx.first;

  Field::ComputingField& Ex = _field.get_Ex();
  Field::ComputingField& Ey = _field.get_Ey();
  Field::ComputingField& Ez = _field.get_Ez();

  Field::ComputingField& Bx = _field.get_Bx();
  Field::ComputingField& By = _field.get_By();
  Field::ComputingField& Bz = _field.get_Bz();

  for (uint64_t i = 0; i < _field.get_Nx(); ++i, x += _field.get_dx())
    for (uint64_t j = 0; j < _field.get_Ny(); ++j) {
      Ez(i, j) =
        sin(2.0 * PI * (x - ax_bx.first) /
          (ax_bx.second - ax_bx.first));
      By(i, j) =
        sin(2.0 * PI * (x - ax_bx.first) /
          (ax_bx.second - ax_bx.first));
      Ex(i, j) = Ey(i, j) = Bx(i, j) = Bz(i, j) = 0.0;
    }
}

void gtest::set_default_field(FDTD::FDTD& _field,
  std::pair<double, double>& ax_bx) {

  double x = ax_bx.first;

  Field::ComputingField& Ex = _field.get_Ex();
  Field::ComputingField& Ey = _field.get_Ey();
  Field::ComputingField& Ez = _field.get_Ez();

  Field::ComputingField& Bx = _field.get_Bx();
  Field::ComputingField& By = _field.get_By();
  Field::ComputingField& Bz = _field.get_Bz();

  for (uint64_t i = 0; i < _field.get_Nx(); ++i, x += _field.get_dx())
    for (uint64_t j = 0; j < _field.get_Ny(); ++j) {
      Ez(i, j) =
        sin(2.0 * PI * (x - ax_bx.first) /
          (ax_bx.second - ax_bx.first));
      By(i, j) = 
        sin(2.0 * PI * (x - ax_bx.first) /
        (ax_bx.second - ax_bx.first));
      Ex(i, j) = Ey(i, j) = Bx(i, j) = Bz(i, j) = 0.0;
    }
}

void gtest::analytical_default_solution(FDTD::FDTD& _field,
  std::pair<double, double>& ax_bx,
  double t) {

  double x = ax_bx.first;

  Field::ComputingField& Ex = _field.get_Ex();
  Field::ComputingField& Ey = _field.get_Ey();
  Field::ComputingField& Ez = _field.get_Ez();

  Field::ComputingField& Bx = _field.get_Bx();
  Field::ComputingField& By = _field.get_By();
  Field::ComputingField& Bz = _field.get_Bz();

  for (uint64_t i = 0ull; i < _field.get_Nx(); ++i, x += _field.get_dx())
    for (uint64_t j = 0ull; j < _field.get_Ny(); ++j) {
      Ez(i, j) =
        sin(2.0 * PI * (x - ax_bx.first + C * t) /
          (ax_bx.second - ax_bx.first));
      By(i, j) = 
        sin(2.0 * PI * (x - ax_bx.first + C * t) /
        (ax_bx.second - ax_bx.first));
      Ex(i, j) = Ey(i, j) = Bx(i, j) = Bz(i, j) = 0.0;
    }
}

double gtest::get_global_err(Field::ComputingField& field_1, Field::ComputingField& field_2)
{
  double max_err = 0.0;
  double tmp_err;
  for (uint64_t j = 0; j < field_1.get_Ny(); ++j)
    for (uint64_t i = 0; i < field_1.get_Nx(); ++i)
    {
      tmp_err = fabs(field_1(i, j) - field_2(i, j));
      if (tmp_err > max_err)
        max_err = tmp_err;
    }
  return max_err;
}



//void gtest::analytical_default_solution(FDTD::FDTD& _field,
//  std::pair<double, double>& ax_bx,
//  double t) {
//
//  double x = ax_bx.first;
//
//  Field::ComputingField& Ex = _field.get_Ex();
//  Field::ComputingField& Ey = _field.get_Ey();
//  Field::ComputingField& Ez = _field.get_Ez();
//
//  Field::ComputingField& Bx = _field.get_Bx();
//  Field::ComputingField& By = _field.get_By();
//  Field::ComputingField& Bz = _field.get_Bz();
//
//  for (uint64_t i = 0ull; i < _field.get_Nx(); ++i, x += _field.get_dx())
//    for (uint64_t j = 0ull; j < _field.get_Ny(); ++j) {
//      Ey(i, j) = Bz(i, j) =
//        sin(2.0 * PI * (x - ax_bx.first - C * t) /
//          (ax_bx.second - ax_bx.first));
//      Ex(i, j) = Ez(i, j) = Bx(i, j) = By(i, j) = 0.0;
//    }
//}

void gtest::Courant_condition_check(const double _dt, const double _dx)
{
  if (_dt <= (_dx / (C * sqrt(2)))) return;
  if (_dt <= (_dx * _dx / 2 * (C * sqrt(2)))) return;
  std::cout << "Courant's condition is not satisfied\n";
  exit(-1);
}

gtest::Test_obj::Test_obj(const Component _E, const Component _B, FDTD::FDTD&& _field) : field(std::move(_field)), analytical_field(field)
{
  Courant_condition_check();
  //std::pair<uint64_t, uint64_t> Nx_Ny = { 64ull, 64ull };
  //std::pair<double, double> ax_ay = { 0.0, 0.0 };
  //std::pair<double, double> bx_by = { 1.0, 1.0 };

  //double dt = 1e-15;
  //FDTD::FDTD copy_field(Nx_Ny, ax_ay, bx_by, dt);
  //field = copy_field;
  E = _E;
  B = _B;
}

gtest::Test_obj::Test_obj(const Test_obj& other_test_field) : field(other_test_field.field), analytical_field(field)
{
  E = other_test_field.E;
  B = other_test_field.B;
}

gtest::Test_obj& gtest::Test_obj::operator=(const Test_obj& other_test_field)
{
  if (this != &other_test_field)
  {
    field = other_test_field.field;
    E = other_test_field.E;
    B = other_test_field.B;
  }
  return *this;
}

gtest::Test_obj& gtest::Test_obj::operator=(Test_obj&& other_test_field) noexcept
{
  if (this != &other_test_field)
  {
    field = std::move(other_test_field.field);
  }
  return *this;
}

void gtest::Test_obj::analytical_default_solution_OX(const Component E, const Component B, const double t)
{
  std::pair<double, double>ax_bx = analytical_field.get_ax_bx();
  double x = ax_bx.first;

  Field::ComputingField& Ex = analytical_field.get_Ex();
  Field::ComputingField& Ey = analytical_field.get_Ey();
  Field::ComputingField& Ez = analytical_field.get_Ez();

  Field::ComputingField& Bx = analytical_field.get_Bx();
  Field::ComputingField& By = analytical_field.get_By();
  Field::ComputingField& Bz = analytical_field.get_Bz();

  double sign = set_sign(E, B);

  if (E == Component::Ey && B == Component::Bz)
  {
    for (uint64_t i = 0ull; i < field.get_Nx(); ++i, x += field.get_dx())
      for (uint64_t j = 0ull; j < field.get_Ny(); ++j) {
        Ey(i, j) =
          sin(2.0 * PI * (x - ax_bx.first - sign * C * t) /
            (ax_bx.second - ax_bx.first));
        Bz(i, j) =
          sin(2.0 * PI * (x - ax_bx.first - sign * C * t) /
            (ax_bx.second - ax_bx.first));
        Ex(i, j) = Ez(i, j) = Bx(i, j) = By(i, j) = 0.0;
      }
  }
  else if (E == Component::Ez && B == Component::By)
  {
    for (uint64_t i = 0ull; i < field.get_Nx(); ++i, x += field.get_dx())
      for (uint64_t j = 0ull; j < field.get_Ny(); ++j) {
        Ez(i, j) =
          sin(2.0 * PI * (x - ax_bx.first - sign * C * t) /
            (ax_bx.second - ax_bx.first));
        By(i, j) =
          sin(2.0 * PI * (x - ax_bx.first - sign * C * t) /
            (ax_bx.second - ax_bx.first));
        Ex(i, j) = Ey(i, j) = Bx(i, j) = Bz(i, j) = 0.0;
      }
  }
  else
  {
    std::cout << "Invalid components!\n";
    exit(-1);
  }
}

void gtest::Test_obj::analytical_default_solution_OY(const Component E, const Component B, const double t)
{
  std::pair<double, double>ay_by = analytical_field.get_ay_by();
  double y = ay_by.first;

  Field::ComputingField& Ex = analytical_field.get_Ex();
  Field::ComputingField& Ey = analytical_field.get_Ey();
  Field::ComputingField& Ez = analytical_field.get_Ez();

  Field::ComputingField& Bx = analytical_field.get_Bx();
  Field::ComputingField& By = analytical_field.get_By();
  Field::ComputingField& Bz = analytical_field.get_Bz();

  double sign = set_sign(E, B);

  if (E == Component::Ez && B == Component::Bx)
  {
    for (uint64_t j = 0ull; j < field.get_Ny(); ++j, y += field.get_dy())
      for (uint64_t i = 0ull; i < field.get_Nx(); ++i)
      {
        Ez(i, j) =
          sin(2.0 * PI * (y - ay_by.first - sign * C * t) /
            (ay_by.second - ay_by.first));
        Bx(i, j) =
          sin(2.0 * PI * (y - ay_by.first - sign * C * t) /
            (ay_by.second - ay_by.first));
        Ex(i, j) = Ey(i, j) = By(i, j) = Bz(i, j) = 0.0;
      }
  }
  else if (E == Component::Ex && B == Component::Bz)
  {
    for (uint64_t j = 0ull; j < field.get_Ny(); ++j, y += field.get_dy())
      for (uint64_t i = 0ull; i < field.get_Nx(); ++i)
      {
        Ex(i, j) =
          sin(2.0 * PI * (y - ay_by.first - sign * C * t) /
            (ay_by.second - ay_by.first));
        Bz(i, j) =
          sin(2.0 * PI * (y - ay_by.first - sign * C * t) /
            (ay_by.second - ay_by.first));
        Ey(i, j) = Ez(i, j) = Bx(i, j) = By(i, j) = 0.0;
      }
  }
  else
  {
    std::cout << "Invalid components!\n";
    exit(-1);
  }
}

void gtest::Test_obj::numerical_solution(const double t)
{
  field.field_update(t);
}

double gtest::Test_obj::get_global_err(const Component component)
{
  if (component == Component::Ex)
  {
    Field::ComputingField& field_component = field.get_Ex();
    Field::ComputingField& analytic_field_component = analytical_field.get_Ex();
    return helper_get_global_err(field_component, analytic_field_component);
  }
  else if (component == Component::Ey)
  {
    Field::ComputingField& field_component = field.get_Ey();
    Field::ComputingField& analytic_field_component = analytical_field.get_Ey();
    return helper_get_global_err(field_component, analytic_field_component);
  }
  else if (component == Component::Ez)
  {
    Field::ComputingField& field_component = field.get_Ez();
    Field::ComputingField& analytic_field_component = analytical_field.get_Ez();
    return helper_get_global_err(field_component, analytic_field_component);
  }
  else if (component == Component::Bx)
  {
    Field::ComputingField& field_component = field.get_Bx();
    Field::ComputingField& analytic_field_component = analytical_field.get_Bx();
    return helper_get_global_err(field_component, analytic_field_component);
  }
  else if (component == Component::By)
  {
    Field::ComputingField& field_component = field.get_By();
    Field::ComputingField& analytic_field_component = analytical_field.get_By();
    return helper_get_global_err(field_component, analytic_field_component);
  }
  else if (component == Component::Bz)
  {
    Field::ComputingField& field_component = field.get_Bz();
    Field::ComputingField& analytic_field_component = analytical_field.get_Bz();
    return helper_get_global_err(field_component, analytic_field_component);
  }
  else
  {
    std::cout << "Get global error: invalid components!\n";
    exit(-1);
  }
}

void gtest::Test_obj::check_ñonvergence(Test_obj& other_test)
{
  std::cout << "\n\n The 1st error is: " << this->get_global_err(E);
  std::cout << "\n The 2nd error is: " << other_test.get_global_err(E);

  std::cout << "\n Difference = " << this->get_global_err(E) / other_test.get_global_err(E) << "\n\n";
}

void gtest::Test_obj::Courant_condition_check() const noexcept
{
  if (field.get_dt() <= (field.get_dx() * field.get_dx() / 2 * (C * sqrt(2)))) return; // Courant's condition of heat conduction
  if (field.get_dt() <= (field.get_dy() * field.get_dy() / 2 * (C * sqrt(2)))) return;
  std::cout << "Courant's condition is not satisfied\n";

  //if (_dt <= (_dx / (C * sqrt(2)))) return;
  //if (_dt <= (_dx * _dx / 2 * (C * sqrt(2)))) return;
  exit(-1);
}

double gtest::Test_obj::helper_get_global_err(const Field::ComputingField& _E, const Field::ComputingField& _B)
{
  double max_err = 0.0;
  double tmp_err;
  for (uint64_t j = 0ull; j < field.get_Ny(); ++j)
    for (uint64_t i = 0ull; i < field.get_Nx(); ++i)
    {
      tmp_err = fabs(_E(i, j) - _B(i, j));
      if (tmp_err > max_err)
        max_err = tmp_err;
    }
  return max_err;
}

void gtest::Test_obj::set_default_field_OX(const Component E, const Component B)
{
  std::pair<double, double>ax_bx = field.get_ax_bx();
  double x = ax_bx.first;

  Field::ComputingField& Ex = field.get_Ex();
  Field::ComputingField& Ey = field.get_Ey();
  Field::ComputingField& Ez = field.get_Ez();

  Field::ComputingField& Bx = field.get_Bx();
  Field::ComputingField& By = field.get_By();
  Field::ComputingField& Bz = field.get_Bz();

  if (E == Component::Ey && B == Component::Bz)
  {
    for (uint64_t i = 0ull; i < field.get_Nx(); ++i, x += field.get_dx())
      for (uint64_t j = 0ull; j < field.get_Ny(); ++j) {
        Ey(i, j) =
          sin(2.0 * PI * (x - ax_bx.first) /
            (ax_bx.second - ax_bx.first));
        Bz(i, j) =
          sin(2.0 * PI * (x - ax_bx.first) /
            (ax_bx.second - ax_bx.first));
        Ex(i, j) = Ez(i, j) = Bx(i, j) = By(i, j) = 0.0;
      }
  }
  else if (E == Component::Ez && B == Component::By)
  {
    for (uint64_t i = 0ull; i < field.get_Nx(); ++i, x += field.get_dx())
      for (uint64_t j = 0ull; j < field.get_Ny(); ++j) {
        Ez(i, j) =
          sin(2.0 * PI * (x - ax_bx.first) /
            (ax_bx.second - ax_bx.first));
        By(i, j) =
          sin(2.0 * PI * (x - ax_bx.first) /
            (ax_bx.second - ax_bx.first));
        Ex(i, j) = Ey(i, j) = Bx(i, j) = Bz(i, j) = 0.0;
      }
  }
  else
  {
    std::cout << "Invalid components!\n";
    exit(-1);
  }
}

void gtest::Test_obj::set_default_field_OY(const Component E, const Component B)
{
  std::pair<double, double>ay_by = field.get_ay_by();

  double y = ay_by.first;

  Field::ComputingField& Ex = field.get_Ex();
  Field::ComputingField& Ey = field.get_Ey();
  Field::ComputingField& Ez = field.get_Ez();

  Field::ComputingField& Bx = field.get_Bx();
  Field::ComputingField& By = field.get_By();
  Field::ComputingField& Bz = field.get_Bz();

  if (E == Component::Ez && B == Component::Bx)
  {
    for (uint64_t j = 0ull; j < field.get_Ny(); ++j, y += field.get_dy())
      for (uint64_t i = 0ull; i < field.get_Nx(); ++i)
      {
        Ez(i, j) =
          sin(2.0 * PI * (y - ay_by.first) /
            (ay_by.second - ay_by.first));
        Bx(i, j) =
          sin(2.0 * PI * (y - ay_by.first) /
            (ay_by.second - ay_by.first));
        Ex(i, j) = Ey(i, j) = By(i, j) = Bz(i, j) = 0.0;
      }
  }
  else if (E == Component::Ex && B == Component::Bz)
  {
    for (uint64_t j = 0ull; j < field.get_Ny(); ++j, y += field.get_dy())
      for (uint64_t i = 0ull; i < field.get_Nx(); ++i)
      {
        Ex(i, j) =
          sin(2.0 * PI * (y - ay_by.first) /
            (ay_by.second - ay_by.first));
        Bz(i, j) =
          sin(2.0 * PI * (y - ay_by.first) /
            (ay_by.second - ay_by.first));
        Ey(i, j) = Ez(i, j) = Bx(i, j) = By(i, j) = 0.0;
      }
  }
  else
  {
    std::cout << "Invalid components!\n";
    exit(-1);
  }
}

double gtest::Test_obj::set_sign(const Component E, const Component B)
{
  if (E == Component::Ey && B == Component::Bz || E == Component::Ez && B == Component::Bx) return 1.0;
  if (E == Component::Ez && B == Component::By || E == Component::Ex && B == Component::Bz) return -1.0;
  std::cout << "Invalid components E and B\n";
  exit(-1);
}
