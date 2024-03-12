#include "tests.hpp"

gtest::Test_obj::Test_obj(const Component _E, const Component _B, FDTD::FDTD& _field) : field(_field), analytical_field(field)
{
  E = _E;
  B = _B;
}

gtest::Test_obj::Test_obj(const Component _E, const Component _B, FDTD::FDTD&& _field) : field(std::move(_field)), analytical_field(field)
{
  //Courant_condition_check();

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
    this->gtest::Test_obj::Test_obj(other_test_field);
  }
  return *this;
}

gtest::Test_obj& gtest::Test_obj::operator=(const FDTD::FDTD& other_test_field)
{
  if (&(this->field) != &other_test_field)
  {
    field = other_test_field;
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

//void gtest::Test_obj::analytical_default_solution_OX(const Component E, const Component B, const double t, const Shift _shift)
//{
//  std::pair<double, double>ax_bx = analytical_field.get_ax_bx();
//  double x = ax_bx.first;
//  double coeff = (_shift == Shift::shifted) ? 0.5 : 0.0;
//
//  Field::ComputingField& Ex = analytical_field.get_Ex();
//  Field::ComputingField& Ey = analytical_field.get_Ey();
//  Field::ComputingField& Ez = analytical_field.get_Ez();
//
//  Field::ComputingField& Bx = analytical_field.get_Bx();
//  Field::ComputingField& By = analytical_field.get_By();
//  Field::ComputingField& Bz = analytical_field.get_Bz();
//
//  double sign = set_sign(E, B);
//
//  if (E == Component::Ey && B == Component::Bz)
//  {
//    for (uint64_t i = 0ull; i < field.get_Nx(); ++i, x += field.get_dx())
//      for (uint64_t j = 0ull; j < field.get_Ny(); ++j) {
//        Ey(i, j) =
//          sin(2.0 * PI * (x - ax_bx.first - sign * C * t) /
//            (ax_bx.second - ax_bx.first));
//        Bz(i, j) =
//          sin(2.0 * PI * (x + field.get_dx() * coeff - ax_bx.first - sign * C * t) /
//            (ax_bx.second - ax_bx.first));
//        Ex(i, j) = Ez(i, j) = Bx(i, j) = By(i, j) = 0.0;
//      }
//  }
//  else if (E == Component::Ez && B == Component::By)
//  {
//    for (uint64_t i = 0ull; i < field.get_Nx(); ++i, x += field.get_dx())
//      for (uint64_t j = 0ull; j < field.get_Ny(); ++j) {
//        Ez(i, j) =
//          sin(2.0 * PI * (x - ax_bx.first - sign * C * t) /
//            (ax_bx.second - ax_bx.first));
//        By(i, j) =
//          sin(2.0 * PI * (x + field.get_dx() * coeff - ax_bx.first - sign * C * t) /
//            (ax_bx.second - ax_bx.first));
//        Ex(i, j) = Ey(i, j) = Bx(i, j) = Bz(i, j) = 0.0;
//      }
//  }
//  else
//  {
//    std::cout << "Invalid components!\n";
//    exit(-1);
//  }
//}

void gtest::Test_obj::analytical_default_solution(const Component E, const Component B, const double t, const Shift _shift)
{
  //std::pair<double, double>ay_by = analytical_field.get_ay_by();
  //double y = ay_by.first;
  double coeff = (_shift == Shift::shifted) ? 0.5 : 0.0;

  Field::ComputingField& Ex = analytical_field.get_Ex();
  Field::ComputingField& Ey = analytical_field.get_Ey();
  Field::ComputingField& Ez = analytical_field.get_Ez();

  Field::ComputingField& Bx = analytical_field.get_Bx();
  Field::ComputingField& By = analytical_field.get_By();
  Field::ComputingField& Bz = analytical_field.get_Bz();


  //double sign = set_sign(E, B);

  //FDTD::FDTD& analytical_field = this->analytical_field;
  auto get_E = [this, &E]() -> Field::ComputingField& {
    if (E == Component::Ex) return analytical_field.get_Ex();
    else if (E == Component::Ey) return analytical_field.get_Ey();
    else if (E == Component::Ez) return analytical_field.get_Ez();
    else {
      std::cout << "Error: Analytical default solution. Wrong E - field!\n";
      exit(-1);
    }
  };
  auto get_B = [this, &B]() -> Field::ComputingField& {
    if (B == Component::Bx) return analytical_field.get_Bx();
    else if (B == Component::By) return analytical_field.get_By();
    else if (B == Component::Bz) return analytical_field.get_Bz();
    else {
      std::cout << "Error: Analytical default solution. Wrong B - field!\n";
      exit(-1);
    }
  };

  std::pair<double, double> ai_bi;
  double coordinate = 0.0;
  double delta_coordinate = 0.0;
  double sign = 0.0;

  
  auto helper_set_data = [&] \
    (std::pair<double, double>&_ai_bi, double delta, double _sign) {
    ai_bi = _ai_bi;
    coordinate = ai_bi.first;
    delta_coordinate = delta;
    sign = _sign;
  };

  Field::ComputingField& E_field = get_E();
  Field::ComputingField& B_field = get_B();

  auto loop_function = [&](std::tuple<Axis, uint64_t, uint64_t> axis_1, std::tuple<Axis, uint64_t, uint64_t> axis_2, std::tuple<Axis, uint64_t, uint64_t> axis_3) {
    uint64_t* i, * j, * k;

    //std::vector<std::tuple<Axis, uint64_t, uint64_t>> params = {axis_1, axis_2, axis_3};
    //for (auto& axis : params)
    //{
    //  switch (std::get<0>(axis))
    //  {
    //  case Axis::Ox: i = &std::get<1>(axis); break;
    //  case Axis::Oy: j = &std::get<1>(axis); break;
    //  case Axis::Oz: k = &std::get<1>(axis); break;
    //  default: break;
    //  }
    //}

    uint64_t axis_1_counter = std::get<1>(axis_1);
    uint64_t axis_2_counter = std::get<1>(axis_2);
    uint64_t axis_3_counter = std::get<1>(axis_3);

    switch (std::get<0>(axis_1))
    {
    case Axis::Ox: i = &axis_1_counter; break;
    case Axis::Oy: j = &axis_1_counter; break;
    case Axis::Oz: k = &axis_1_counter; break;
    default: break;
    }
    switch (std::get<0>(axis_2))
    {
    case Axis::Ox: i = &axis_2_counter; break;
    case Axis::Oy: j = &axis_2_counter; break;
    case Axis::Oz: k = &axis_2_counter; break;
    default: break;
    }
    switch (std::get<0>(axis_3))
    {
    case Axis::Ox: i = &axis_3_counter; break;
    case Axis::Oy: j = &axis_3_counter; break;
    case Axis::Oz: k = &axis_3_counter; break;
    default: break;
    }

    for (axis_1_counter = 0ull; axis_1_counter < std::get<2>(axis_1); ++(axis_1_counter), coordinate += delta_coordinate)
      for (axis_2_counter = 0ull; axis_2_counter < std::get<2>(axis_2); ++(axis_2_counter))
        for (axis_3_counter = 0ull;  axis_3_counter < std::get<2>(axis_3); ++(axis_3_counter))
        {
          Ex(*i, *j, *k) = Ey(*i, *j, *k) = Ez(*i, *j, *k) = Bx(*i, *j, *k) = By(*i, *j, *k) = Bz(*i, *j, *k) = 0.0;

          E_field(*i, *j, *k) =
            sin(2.0 * PI * (coordinate - ai_bi.first - sign * C * t) /
              (ai_bi.second - ai_bi.first));
          B_field(*i, *j, *k) =
            sin(2.0 * PI * (coordinate + delta_coordinate * coeff - ai_bi.first - sign * C * t) /
              (ai_bi.second - ai_bi.first));
        }
  };

  auto set_computational_data = [this, &helper_set_data, &loop_function, &E, &B]() {
    // OX
    if (E == Component::Ey && B == Component::Bz)
    {
      helper_set_data(field.get_ax_bx(), field.get_dx(), 1.0);
      loop_function(std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Oz, 0, field.get_Nz()));
    }
    else if (E == Component::Ez && B == Component::By)
    {
      helper_set_data(field.get_ax_bx(), field.get_dx(), -1.0);
      loop_function(std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Oz, 0, field.get_Nz()));
    }

    // OY
    else if (E == Component::Ez && B == Component::Bx)
    {
      helper_set_data(field.get_ay_by(), field.get_dy(), 1.0);
      loop_function(std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oz, 0, field.get_Nz()));
    }
    else if (E == Component::Ex && B == Component::Bz)
    {
      helper_set_data(field.get_ay_by(), field.get_dy(), -1.0);
      loop_function(std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oz, 0, field.get_Nz()));
    }

    // OZ
    else if (E == Component::Ex && B == Component::By)
    {
      helper_set_data(field.get_az_bz(), field.get_dz(), 1.0);
      loop_function(std::make_tuple(Axis::Oz, 0, field.get_Nz()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()));
    }
    else if (E == Component::Ey && B == Component::Bx)
    {
      helper_set_data(field.get_az_bz(), field.get_dz(), -1.0);
      loop_function(std::make_tuple(Axis::Oz, 0, field.get_Nz()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()));
    }
    else {
      std::cout << "Error: Analytical solution. Invalid E/B components!\n";
      exit(-1);
    }
  };
  
  set_computational_data();


  //if (E == Component::Ez && B == Component::Bx)
  //{
    //for (uint64_t j = 0ull; j < field.get_Ny(); ++j, y += field.get_dy())
    //{
    //  for (uint64_t i = 0ull; i < field.get_Nx(); ++i)
    //    for (uint64_t k = 0ull; k < field.get_Nz(); ++k)
    //  {
    //    Ex(i, j, k) = Ey(i, j, k) = Ez(i, j, k) = Bx(i, j, k) = By(i, j, k) = Bz(i, j, k) = 0.0;
    //   
    //    E_field(i, j, k) =
    //      sin(2.0 * PI * (y - ay_by.first - sign * C * t) /
    //        (ay_by.second - ay_by.first));
    //    B_field(i, j, k) =
    //      sin(2.0 * PI * (y + field.get_dy() * coeff - ay_by.first - sign * C * t) /
    //        (ay_by.second - ay_by.first));

    //    //Ez(i, j) =
    //    //  sin(2.0 * PI * (y - ay_by.first - sign * C * t) /
    //    //    (ay_by.second - ay_by.first));
    //    //Bx(i, j) =
    //    //  sin(2.0 * PI * (y + field.get_dy() * coeff - ay_by.first - sign * C * t) /
    //    //    (ay_by.second - ay_by.first));
    //    //Ex(i, j) = Ey(i, j) = By(i, j) = Bz(i, j) = 0.0;
    //  }
    //}
  //}

  

  /*else if (E == Component::Ex && B == Component::Bz)
  {
    for (uint64_t j = 0ull; j < field.get_Ny(); ++j, y += field.get_dy())
    {

      for (uint64_t i = 0ull; i < field.get_Nx(); ++i)
      {
        Ex(i, j) =
          sin(2.0 * PI * (y - ay_by.first - sign * C * t) /
            (ay_by.second - ay_by.first));
        Bz(i, j) =
          sin(2.0 * PI * (y + field.get_dy() * coeff - ay_by.first - sign * C * t) /
            (ay_by.second - ay_by.first));
        Ey(i, j) = Ez(i, j) = Bx(i, j) = By(i, j) = 0.0;
      }
    }
  }*/
  //else
  //{
  //  std::cout << "Invalid components!\n";
  //  exit(-1);
  //}
}

void gtest::Test_obj::numerical_solution(const double t, const Shift _shift)
{
  Courant_condition_check(_shift);
  (_shift == Shift::shifted) ? field.shifted_field_update(t) : field.field_update(t);
}

void gtest::Test_obj::numerical_solution(const uint64_t t, const Shift _shift)
{
  Courant_condition_check(_shift);
  (_shift == Shift::shifted) ? field.shifted_field_update(t) : field.field_update(t);
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

void gtest::Test_obj::print_convergence(Test_obj& other_test)
{
  std::cout << "\n\n================================\n The 1st error is: " << this->get_global_err(E);
  std::cout << "\n The 2nd error is: " << other_test.get_global_err(E);

  std::cout << "\n Difference(E) = " << this->get_global_err(E) / other_test.get_global_err(E);
  std::cout << "\n Difference(B) = " << this->get_global_err(B) / other_test.get_global_err(B) << "\n================================" << "\n\n";
}

void gtest::Test_obj::write_convergence_to_file(const char* path, std::vector<double>& data)
{
  std::ofstream outfile;
  outfile.open(path);
  if (!outfile.is_open())
  {
    std::cout << "The file can't be opened!" << std::endl;
    exit(-1);
  }
  uint64_t i = 0ull;

  for (i = 0ull; i < data.size() - 1ull; ++i)
  {
    outfile << data[i] << ';';
  }
  outfile << data[i] << std::endl;

  //for (const auto& it : data)
  //{
  //  outfile << it << ' ';
  //}
  //outfile << std::endl;
  
  outfile.close();
}

void gtest::Test_obj::Courant_condition_check(const Shift _shift) const noexcept
{
  if (_shift == Shift::unshifted)
  {
    //if (field.get_dt() <= (field.get_dx() * field.get_dx() / 2 * (C * sqrt(2)))) return; // Courant's condition of heat conduction
    //if (field.get_dt() <= (field.get_dy() * field.get_dy() / 2 * (C * sqrt(2)))) return;
    if (field.get_dt() <= (field.get_dx() * field.get_dx() / 2 * (C * sqrt(2))) && field.get_dt() <= (field.get_dy() * field.get_dy() / 2 * (C * sqrt(2)))) return;
  }
  else
  {
    double tmp = (field.get_dx() / (C * sqrt(2)));
    //if (field.get_dt() <= (field.get_dx() / (C * sqrt(2)))) return;
    //if (field.get_dt() <= (field.get_dy() / (C * sqrt(2)))) return;
    if (field.get_dt() <= (field.get_dx() / (C * sqrt(2))) && field.get_dt() <= (field.get_dy() / (C * sqrt(2)))) return;
  }
  std::cout << "Courant's condition is not satisfied\n";

  //if (_dt <= (_dx / (C * sqrt(2)))) return;
  //if (_dt <= (_dx * _dx / 2 * (C * sqrt(2)))) return;
  exit(-1);
}

double gtest::Test_obj::helper_get_global_err(const Field::ComputingField& field_1, const Field::ComputingField& field_2)
{
  double max_err = 0.0;
  //double tmp_err;
  for (uint64_t j = 0ull; j < field.get_Ny(); ++j)
    for (uint64_t i = 0ull; i < field.get_Nx(); ++i)
    {
      //tmp_err = fabs(field_1(i, j) - field_2(i, j));
      //if (tmp_err > max_err)
      //  max_err = tmp_err;
      max_err = std::max(max_err, fabs(field_1(i, j) - field_2(i, j)));
    }
  return max_err;
}

void gtest::Test_obj::set_default_field_OX(const Component E, const Component B, const Shift _shift)
{
  std::pair<double, double>ax_bx = field.get_ax_bx();
  double x = ax_bx.first;
  double coeff = (_shift == Shift::shifted) ? 0.5 : 0.0;

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
          sin(2.0 * PI * (x + field.get_dx() * coeff - ax_bx.first) /
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
          sin(2.0 * PI * (x + field.get_dx() * coeff - ax_bx.first) /
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

void gtest::Test_obj::set_default_field(const Component E, const Component B, const Shift _shift)
{
  std::pair<double, double>ay_by = field.get_ay_by();

  double y = ay_by.first;
  double coeff = (_shift == Shift::shifted) ? 0.5 : 0.0;

  Field::ComputingField& Ex = field.get_Ex();
  Field::ComputingField& Ey = field.get_Ey();
  Field::ComputingField& Ez = field.get_Ez();

  Field::ComputingField& Bx = field.get_Bx();
  Field::ComputingField& By = field.get_By();
  Field::ComputingField& Bz = field.get_Bz();


  auto get_E = [this, &E]() -> Field::ComputingField& {
    if (E == Component::Ex) return field.get_Ex();
    else if (E == Component::Ey) return field.get_Ey();
    else if (E == Component::Ez) return field.get_Ez();
    else {
      std::cout << "Error: Set default data. Wrong E - field!\n";
      exit(-1);
    }
  };
  auto get_B = [this, &B]() -> Field::ComputingField& {
    if (B == Component::Bx) return field.get_Bx();
    else if (B == Component::By) return field.get_By();
    else if (B == Component::Bz) return field.get_Bz();
    else {
      std::cout << "Error: Set default data. Wrong B - field!\n";
      exit(-1);
    }
  };

  std::pair<double, double> ai_bi;
  double coordinate = 0.0;
  double delta_coordinate = 0.0;
  double sign = 0.0;


  auto helper_set_data = [&] \
    (std::pair<double, double>&_ai_bi, double delta, double _sign) {
    ai_bi = _ai_bi;
    coordinate = ai_bi.first;
    delta_coordinate = delta;
    sign = _sign;
  };

  Field::ComputingField& E_field = get_E();
  Field::ComputingField& B_field = get_B();

  auto loop_function = [&](std::tuple<Axis, uint64_t, uint64_t> axis_1, std::tuple<Axis, uint64_t, uint64_t> axis_2, std::tuple<Axis, uint64_t, uint64_t> axis_3) {
    uint64_t* i, * j, * k;

    uint64_t axis_1_counter = std::get<1>(axis_1);
    uint64_t axis_2_counter = std::get<1>(axis_2);
    uint64_t axis_3_counter = std::get<1>(axis_3);

    switch (std::get<0>(axis_1))
    {
    case Axis::Ox: i = &axis_1_counter; break;
    case Axis::Oy: j = &axis_1_counter; break;
    case Axis::Oz: k = &axis_1_counter; break;
    default: break;
    }
    switch (std::get<0>(axis_2))
    {
    case Axis::Ox: i = &axis_2_counter; break;
    case Axis::Oy: j = &axis_2_counter; break;
    case Axis::Oz: k = &axis_2_counter; break;
    default: break;
    }
    switch (std::get<0>(axis_3))
    {
    case Axis::Ox: i = &axis_3_counter; break;
    case Axis::Oy: j = &axis_3_counter; break;
    case Axis::Oz: k = &axis_3_counter; break;
    default: break;
    }

    for (axis_1_counter = 0ull; axis_1_counter < std::get<2>(axis_1); ++(axis_1_counter), coordinate += delta_coordinate)
      for (axis_2_counter = 0ull; axis_2_counter < std::get<2>(axis_2); ++(axis_2_counter))
        for (axis_3_counter = 0ull; axis_3_counter < std::get<2>(axis_3); ++(axis_3_counter))
        {
          Ex(*i, *j, *k) = Ey(*i, *j, *k) = Ez(*i, *j, *k) = Bx(*i, *j, *k) = By(*i, *j, *k) = Bz(*i, *j, *k) = 0.0;

          E_field(*i, *j, *k) =
            sin(2.0 * PI * (coordinate - ai_bi.first) /
              (ai_bi.second - ai_bi.first));
          B_field(*i, *j, *k) =
            sin(2.0 * PI * (coordinate + delta_coordinate * coeff - ai_bi.first) /
              (ai_bi.second - ai_bi.first));
        }
  };

  auto set_computational_data = [this, &helper_set_data, &loop_function, &E, &B]() {
    // OX
    if (E == Component::Ey && B == Component::Bz)
    {
      helper_set_data(field.get_ax_bx(), field.get_dx(), 1.0);
      loop_function(std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Oz, 0, field.get_Nz()));
    }
    else if (E == Component::Ez && B == Component::By)
    {
      helper_set_data(field.get_ax_bx(), field.get_dx(), -1.0);
      loop_function(std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Oz, 0, field.get_Nz()));
    }

    // OY
    else if (E == Component::Ez && B == Component::Bx)
    {
      helper_set_data(field.get_ay_by(), field.get_dy(), 1.0);
      loop_function(std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oz, 0, field.get_Nz()));
    }
    else if (E == Component::Ex && B == Component::Bz)
    {
      helper_set_data(field.get_ay_by(), field.get_dy(), -1.0);
      loop_function(std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oz, 0, field.get_Nz()));
    }

    // OZ
    else if (E == Component::Ex && B == Component::By)
    {
      helper_set_data(field.get_az_bz(), field.get_dz(), 1.0);
      loop_function(std::make_tuple(Axis::Oz, 0, field.get_Nz()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()));
    }
    else if (E == Component::Ey && B == Component::Bx)
    {
      helper_set_data(field.get_az_bz(), field.get_dz(), -1.0);
      loop_function(std::make_tuple(Axis::Oz, 0, field.get_Nz()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()));
    }
    else {
      std::cout << "Error: Set default data. Invalid E/B components!\n";
      exit(-1);
    }
  };

  set_computational_data();



  //if (E == Component::Ez && B == Component::Bx)
  //{
  //  for (uint64_t j = 0ull; j < field.get_Ny(); ++j, y += field.get_dy())
  //    for (uint64_t i = 0ull; i < field.get_Nx(); ++i)
  //      for (uint64_t k = 0ull; k < field.get_Nz(); ++k)
  //    {
  //      Ez(i, j, k) =
  //        sin(2.0 * PI * (y - ay_by.first) /
  //          (ay_by.second - ay_by.first));
  //      Bx(i, j, k) =
  //        sin(2.0 * PI * (y + field.get_dy() * coeff - ay_by.first) /
  //          (ay_by.second - ay_by.first));
  //      Ex(i, j, k) = Ey(i, j, k) = By(i, j, k) = Bz(i, j, k) = 0.0;
  //    }
  //}
  //else if (E == Component::Ex && B == Component::Bz)
  //{
  //  for (uint64_t j = 0ull; j < field.get_Ny(); ++j, y += field.get_dy())
  //    for (uint64_t i = 0ull; i < field.get_Nx(); ++i)
  //    {
  //      Ex(i, j) =
  //        sin(2.0 * PI * (y - ay_by.first) /
  //          (ay_by.second - ay_by.first));
  //      Bz(i, j) =
  //        sin(2.0 * PI * (y + field.get_dy() * coeff - ay_by.first) /
  //          (ay_by.second - ay_by.first));
  //      Ey(i, j) = Ez(i, j) = Bx(i, j) = By(i, j) = 0.0;
  //    }
  //}
  //else
  //{
  //  std::cout << "Invalid components!\n";
  //  exit(-1);
  //}
}

//double gtest::Test_obj::set_sign(const Component E, const Component B)
//{
//  if (E == Component::Ey && B == Component::Bz || E == Component::Ez && B == Component::Bx) return 1.0;
//  if (E == Component::Ez && B == Component::By || E == Component::Ex && B == Component::Bz) return -1.0;
//  std::cout << "Invalid components E and B\n";
//  exit(-1);
//}
