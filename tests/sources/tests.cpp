#include "tests.hpp"

gtest::Test_obj::Test_obj(const Component _E, const Component _B, FDTD::FDTD& _field) : field(_field), analytical_field(field)
{
  E = _E;
  B = _B;
}

gtest::Test_obj::Test_obj(const Component _E, const Component _B, FDTD::FDTD&& _field) : field(std::move(_field)), analytical_field(field)
{
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

void gtest::Test_obj::analytical_default_solution(const Component E, const Component B, const double t, const Shift _shift)
{
  double coeff = (_shift == Shift::shifted) ? 0.5 : 0.0;

  Field::ComputingField& Ex = analytical_field.get_Ex();
  Field::ComputingField& Ey = analytical_field.get_Ey();
  Field::ComputingField& Ez = analytical_field.get_Ez();

  Field::ComputingField& Bx = analytical_field.get_Bx();
  Field::ComputingField& By = analytical_field.get_By();
  Field::ComputingField& Bz = analytical_field.get_Bz();


  auto get_E = [this, &E]() -> Field::ComputingField& {
    if (E == Component::Ex) return analytical_field.get_Ex();
    else if (E == Component::Ey) return analytical_field.get_Ey();
    else if (E == Component::Ez) return analytical_field.get_Ez();
    else {
      std::cout << "\nError: Analytical default solution. Wrong E - field!\n";
      exit(-1);
    }
  };
  auto get_B = [this, &B]() -> Field::ComputingField& {
    if (B == Component::Bx) return analytical_field.get_Bx();
    else if (B == Component::By) return analytical_field.get_By();
    else if (B == Component::Bz) return analytical_field.get_Bz();
    else {
      std::cout << "\nError: Analytical default solution. Wrong B - field!\n";
      exit(-1);
    }
  };

  std::pair<double, double> ai_bi{};
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
    int64_t* i = nullptr;
    int64_t* j = nullptr;
    int64_t* k = nullptr;


    int64_t axis_1_counter = std::get<1>(axis_1);
    int64_t axis_2_counter = std::get<1>(axis_2);
    int64_t axis_3_counter = std::get<1>(axis_3);


    switch (std::get<0>(axis_1))
    {
    case Axis::Ox: i = &axis_1_counter; j = &axis_1_counter; k = &axis_1_counter; break;
    case Axis::Oy: j = &axis_1_counter; i = &axis_1_counter; k = &axis_1_counter; break;
    case Axis::Oz: k = &axis_1_counter; i = &axis_1_counter; j = &axis_1_counter; break;
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

    std::cout << std::get<2>(axis_1) << '\n';
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
          //std::cout << "i = " << *i << " j = " << *j << " k = " << *k << '\n';
        }
    // Обновляю границу - низ 2-х мерной системы
    //std::memcpy(E_field.data() + 1, E_field.data() + E_field.get_Nx() * (E_field.get_Ny() - 1), E_field.get_Nx() * sizeof(double));
    //std::memcpy(B_field.data() + 1, B_field.data() + B_field.get_Nx() * (B_field.get_Ny() - 1), B_field.get_Nx() * sizeof(double));

    //int sc = 0;
    //for (int i = 0; i < E_field.get_Nx() - 1; ++i)
    //{
    //  for (int j = 0; j < E_field.get_Ny() - 1; ++j)
    //  {
    //    E_field(i, j) = sc++;
    //  }
    //}

    //int c = 0;
    //for (int x = 0; x < E_field.field.size(); x++)
    //{
    //  if (c == 10) std::cout << '\n';
    //  std::cout << E_field.field[x] << '\t';
    //}

    //for (int j = 0; j < E_field.get_Ny() + 1; ++j)
    //{
    //  int c = 0;
    //  for (int i = 0; i < E_field.get_Nx() + 1; ++i)
    //  {
    //    std::cout << E_field(i, j) << '\t';
    //    if (c == 10) std::cout << '\n';
    //    c++;
    //  }
    //}
    std::cout << "\n\n";

    // Обновляю границу - верх 2-х мерной системы
    for (uint64_t x = 0; x < E_field.get_Nx(); ++x)
    {
      *(E_field.data() + x + 1) = E_field(x, E_field.get_Ny() - 1); // снизу
      E_field(x, E_field.get_Ny()) = E_field(x, 0); // сверху

      *(B_field.data() + x + 1) = B_field(x, B_field.get_Ny() - 1); // снизу
      B_field(x, B_field.get_Ny()) = B_field(x, 0); // сверху
    }
    for (uint64_t y = 0; y < E_field.get_Ny(); ++y)
    {
      E_field(-1, y) = E_field(E_field.get_Nx() - 1, y); // слева
      E_field(E_field.get_Nx(), y) = E_field(0, y); // справа

      B_field(-1, y) = B_field(B_field.get_Nx() - 1, y); // слева
      B_field(B_field.get_Nx(), y) = B_field(0, y); // справа
    }
    E_field(-1, -1) = E_field(E_field.get_Nx() - 1, E_field.get_Ny() - 1); // левый нижний узел
    E_field(-1, E_field.get_Ny()) = E_field(E_field.get_Nx() - 1, 0); // левый верхний узел
    E_field(E_field.get_Nx(), E_field.get_Ny()) = E_field(0, 0); // правый верхний узел
    *(E_field.data() + E_field.get_Nx() + 1) = E_field(0, E_field.get_Ny() - 1); // правый нижний узел

    B_field(-1, -1) = B_field(B_field.get_Nx() - 1, B_field.get_Ny() - 1); // левый нижний узел
    B_field(-1, B_field.get_Ny()) = B_field(B_field.get_Nx() - 1, 0); // левый верхний узел
    B_field(B_field.get_Nx(), B_field.get_Ny()) = B_field(0, 0); // правый верхний узел
    *(B_field.data() + B_field.get_Nx() + 1) = B_field(0, B_field.get_Ny() - 1); // правый нижний узел
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
  auto get_err = [this](Field::ComputingField& numerical_field, Field::ComputingField& analytical_field) {
    double max_err = 0.0;
    for (uint64_t i = 0ull; i < field.get_Nx(); ++i)
      for (uint64_t j = 0ull; j < field.get_Ny(); ++j)
        for (uint64_t k = 0ull; k < field.get_Nz(); ++k)
          max_err = std::max(max_err, fabs(numerical_field(i, j, k) - analytical_field(i, j, k)));
    
    return max_err;
  };
  switch (component)
  {
  case Component::Ex:
    get_err(field.get_Ex(), analytical_field.get_Ex());
    break;
  case Component::Ey:
    get_err(field.get_Ey(), analytical_field.get_Ey());
    break;
  case Component::Ez:
    get_err(field.get_Ez(), analytical_field.get_Ez());
    break;
  case Component::Bx:
    get_err(field.get_Bx(), analytical_field.get_Bx());
    break;
  case Component::By:
    get_err(field.get_By(), analytical_field.get_By());
    break;
  case Component::Bz:
    get_err(field.get_Bz(), analytical_field.get_Bz());
    break;
  default:
    std::cout << "\nGet global error: Error! Wrong component!\n";
    exit(-1);
  }
}

void gtest::Test_obj::print_convergence(Test_obj& other_test)
{
  double this_E_error = 5;
  this_E_error = this->get_global_err(this->E);
  double this_B_error = this->get_global_err(this->B);

  double other_E_error = other_test.get_global_err(other_test.E);
  double other_B_error = other_test.get_global_err(other_test.B);

  std::cout << "this_E_error = " << this_E_error << '\n';
  std::cout << "\n\n================================\n The 1st (E) error is: " << this_E_error;
  std::cout << "\n The 2nd (E) error is: " << other_E_error;

  std::cout << "\n Difference(E) = " << this_E_error / other_E_error;
  std::cout << "\n Difference(B) = " << this_B_error / other_B_error << "\n================================" << "\n\n";
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
      std::cout << "\nError: Set default data. Wrong E - field!\n";
      exit(-1);
    }
  };
  auto get_B = [this, &B]() -> Field::ComputingField& {
    if (B == Component::Bx) return field.get_Bx();
    else if (B == Component::By) return field.get_By();
    else if (B == Component::Bz) return field.get_Bz();
    else {
      std::cout << "\nError: Set default data. Wrong B - field!\n";
      exit(-1);
    }
  };

  std::pair<double, double> ai_bi{};
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
    uint64_t* i = nullptr, * j = nullptr, * k = nullptr;

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

    // Обновляю границу - верх 2-х мерной системы
    for (uint64_t x = 0; x < E_field.get_Nx(); ++x)
    {
      *(E_field.data() + x + 1) = E_field(x, E_field.get_Ny() - 1); // снизу
      E_field(x, E_field.get_Ny()) = E_field(x, 0); // сверху

      *(B_field.data() + x + 1) = B_field(x, B_field.get_Ny() - 1); // снизу
      B_field(x, B_field.get_Ny()) = B_field(x, 0); // сверху
    }
    for (uint64_t y = 0; y < E_field.get_Ny(); ++y)
    {
      E_field(-1, y) = E_field(E_field.get_Nx() - 1, y); // слева
      E_field(E_field.get_Nx(), y) = E_field(0, y); // справа

      B_field(-1, y) = B_field(B_field.get_Nx() - 1, y); // слева
      B_field(B_field.get_Nx(), y) = B_field(0, y); // справа
    }
    E_field(-1, -1) = E_field(E_field.get_Nx() - 1, E_field.get_Ny() - 1); // левый нижний узел
    E_field(-1, E_field.get_Ny()) = E_field(E_field.get_Nx() - 1, 0); // левый верхний узел
    E_field(E_field.get_Nx(), E_field.get_Ny()) = E_field(0, 0); // правый верхний узел
    *(E_field.data() + E_field.get_Nx() + 1) = E_field(0, E_field.get_Ny() - 1); // правый нижний узел

    B_field(-1, -1) = B_field(B_field.get_Nx() - 1, B_field.get_Ny() - 1); // левый нижний узел
    B_field(-1, B_field.get_Ny()) = B_field(B_field.get_Nx() - 1, 0); // левый верхний узел
    B_field(B_field.get_Nx(), B_field.get_Ny()) = B_field(0, 0); // правый верхний узел
    *(B_field.data() + B_field.get_Nx() + 1) = B_field(0, B_field.get_Ny() - 1); // правый нижний узел
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
      std::cout << "\nError: Set default data. Invalid E/B components!\n";
      exit(-1);
    }
  };

  set_computational_data();
}
