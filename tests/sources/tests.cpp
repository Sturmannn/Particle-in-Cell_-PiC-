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
  _field.clear_fields();
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
    // this->gtest::Test_obj::Test_obj(other_test_field);
    field = other_test_field.field;
    analytical_field = other_test_field.analytical_field;
    E = other_test_field.E;
    B = other_test_field.B;
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
    analytical_field = std::move(other_test_field.analytical_field);
    
    // Clear fields
    other_test_field.field.clear_fields();
    other_test_field.analytical_field.clear_fields();
  }
  return *this;
}

void gtest::Test_obj::analytical_default_solution(const Component E, const Component B, const double t, const Shift _shift)
{
  // Условие куранта проверяется на всех процессах, и там есть exit(-1) в случае ошибки
  Courant_condition_check(_shift);

  double coeff = (_shift == Shift::shifted) ? 0.5 : 0.0;

  Field::ComputingField& Ex = analytical_field.get_Ex();
  Field::ComputingField& Ey = analytical_field.get_Ey();
  Field::ComputingField& Ez = analytical_field.get_Ez();

  Field::ComputingField& Bx = analytical_field.get_Bx();
  Field::ComputingField& By = analytical_field.get_By();
  Field::ComputingField& Bz = analytical_field.get_Bz();


  auto get_E = [this, &E]() -> Field::ComputingField& {
    if (E == Component::Ex) return analytical_field.get_Ex();
    if (E == Component::Ey) return analytical_field.get_Ey();
    if (E == Component::Ez) return analytical_field.get_Ez();
    std::cout << "\nError: Analytical default solution. Wrong E - field!\n";
    exit(-1);
  };
  auto get_B = [this, &B]() -> Field::ComputingField& {
    if (B == Component::Bx) return analytical_field.get_Bx();
    if (B == Component::By) return analytical_field.get_By();
    if (B == Component::Bz) return analytical_field.get_Bz();
    std::cout << "\nError: Analytical default solution. Wrong B - field!\n";
    exit(-1);
  };

  std::pair<double, double> ai_bi{};
  double coordinate = 0.0;
  double delta_coordinate = 0.0;
  double sign = 0.0;

  int rank = 0;
  int world_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  auto helper_set_data = [&] \
    (std::pair<double, double>&&_ai_bi, double delta, double _sign) {
    // (std::pair<double, double>&_ai_bi, double delta, double _sign) {
    ai_bi = _ai_bi;
    
    std::cout << "rank = " << rank << " Nx = " << field.get_Nx() << " Ny = " << field.get_Ny() << " Nz = " << field.get_Nz() << " delta = " << delta << std::endl;

    // coordinate = ai_bi.first;
    delta_coordinate = delta;

    // Вектор значений соответствующих размеров (Nx или Ny или Nz) для каждого процесса
    std::vector<int64_t> all_N(world_size);
    int64_t shift = 0;

    // ЗДЕСЬ РАЗБИЕНИЕ ТОЛЬКО ПО СТРОКАМ (OY) В ДВУМЕРНОМ СЛУЧАЕ И ТОЛЬКО ПРИ ДВИЖЕНИИ ПО OY. ДАЛЬШЕ НЕОБХОДИМО БУДЕТ МЕНЯТЬ, ЕСЛИ ПО РАЗНЫМ ОСЯМ РАЗБИВАТЬ
    Axis axis = FDTD::FDTD::get_axis(E, B);
    // Опасное место, так как шаг может быть недостаточно сдвинут.
    // Насколько я понял, shift и coordinate должны вычисляться по-разному, в зависимости от выбранной оси ДВИЖЕНИЯ
    // Это необходимо, так как на данный момент разбиение по процессам происходит только по OY
    if (axis == Axis::Ox) {
      int64_t local_Nx = field.get_Nx();
      MPI_Allgather(&local_Nx, 1, MPI_INT64_T, all_N.data(), 1, MPI_INT64_T, MPI_COMM_WORLD);
      for (int i = 0; i < rank; ++i) shift += all_N[i];
      // coordinate = ai_bi.first + delta_coordinate * shift/* + delta_coordinate*/; // Для нескольких процессов
    }
    else if (axis == Axis::Oy) {
      int64_t local_Ny = field.get_Ny();
      MPI_Allgather(&local_Ny, 1, MPI_INT64_T, all_N.data(), 1, MPI_INT64_T, MPI_COMM_WORLD);
      for (int i = 0; i < rank; ++i) shift += all_N[i];
      coordinate = ai_bi.first + delta_coordinate * shift/* + delta_coordinate*/; // Для нескольких процессов
    }
    // С OZ пока вообще непонятно, так как пока что не предполагается движение по OZ
    else if (axis == Axis::Oz) {
      int64_t local_Nz = field.get_Nz();
      MPI_Allgather(&local_Nz, 1, MPI_INT64_T, all_N.data(), 1, MPI_INT64_T, MPI_COMM_WORLD);
      for (int i = 0; i < rank; ++i) shift += all_N[i];
      coordinate = ai_bi.first + delta_coordinate * shift/* + delta_coordinate*/; // Для нескольких процессов
    }

    // std::cout << "rank = " << rank << " coordinate = " << coordinate << std::endl;
    sign = _sign;
  };

  Field::ComputingField& E_field = get_E();
  Field::ComputingField& B_field = get_B();
  
  // Функция синхронизации всех границ для 2D и 3D
  std::function<void(void)> boundary_synchronization{};
  if (this->analytical_field.get_Nz() > 1)
    boundary_synchronization = std::bind(&FDTD::FDTD::boundary_synchronization_3D, &analytical_field);
  else
    boundary_synchronization = std::bind(&FDTD::FDTD::boundary_synchronization, &analytical_field);

  auto loop_function = [&](std::tuple<Axis, int64_t, int64_t> axis_1, std::tuple<Axis, int64_t, int64_t> axis_2, std::tuple<Axis, int64_t, int64_t> axis_3) {
    int64_t* i = nullptr;
    int64_t* j = nullptr;
    int64_t* k = nullptr;


    int64_t axis_1_counter = std::get<1>(axis_1);
    int64_t axis_2_counter = std::get<1>(axis_2);
    int64_t axis_3_counter = std::get<1>(axis_3);

    switch (std::get<0>(axis_1))
    {
    case Axis::Ox: 
      i = &axis_1_counter; 
      // axis_1_counter += rank * field.get_Nx();
      break;
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

    for (axis_1_counter = std::get<1>(axis_1); axis_1_counter < std::get<2>(axis_1); ++(axis_1_counter), coordinate += delta_coordinate)
      for (axis_2_counter = std::get<1>(axis_2); axis_2_counter < std::get<2>(axis_2); ++(axis_2_counter))
        for (axis_3_counter = std::get<1>(axis_3); axis_3_counter < std::get<2>(axis_3); ++(axis_3_counter))
        {
          Ex(*i, *j, *k) = Ey(*i, *j, *k) = Ez(*i, *j, *k) = Bx(*i, *j, *k) = By(*i, *j, *k) = Bz(*i, *j, *k) = 0.0;

          E_field(*i, *j, *k) =
            sin(2.0 * PI * (coordinate - ai_bi.first - sign * C * t) /
              (ai_bi.second - ai_bi.first));
          B_field(*i, *j, *k) =
            sin(2.0 * PI * (coordinate + delta_coordinate * coeff - ai_bi.first - sign * C * t) /
              (ai_bi.second - ai_bi.first));
        }
    boundary_synchronization();
  };

  auto set_computational_data = [&]() {

    // Все эти _Nz пока неизменные, для MPI возможно потребуется изменить, если нужно будет разбивать оси по Oz
    // Может быть эти _Nz здесь не нужны, так как проверка в самом конструкторе. После готовой реализации MPI следует проверить.
    int64_t _Nz = (field.get_Nz() > 1) ? field.get_Nz() : 0; // Учитываем случай 2D
    int64_t _Nz_counter = (_Nz > 1) ? 0 : -1; // Это эквивалентно -1 для 2D

    // OX
    if (E == Component::Ey && B == Component::Bz)
    {
      helper_set_data(field.get_ax_bx(), field.get_dx(), 1.0);
      loop_function(std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }
    else if (E == Component::Ez && B == Component::By)
    {
      helper_set_data(field.get_ax_bx(), field.get_dx(), -1.0);
      loop_function(std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }

    // OY
    else if (E == Component::Ez && B == Component::Bx)
    {
      helper_set_data(field.get_ay_by(), field.get_dy(), 1.0);
      loop_function(std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }
    else if (E == Component::Ex && B == Component::Bz)
    {
      helper_set_data(field.get_ay_by(), field.get_dy(), -1.0);
      loop_function(std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }

    // OZ
    // Возможно нужно учесть get_dz() в случае 2D
    else if (E == Component::Ex && B == Component::By)
    {
      helper_set_data(field.get_az_bz(), field.get_dz(), 1.0);
      loop_function(std::make_tuple(Axis::Oz, _Nz_counter, _Nz),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()));
    }
    else if (E == Component::Ey && B == Component::Bx)
    {
      helper_set_data(field.get_az_bz(), field.get_dz(), -1.0);
      loop_function(std::make_tuple(Axis::Oz, _Nz_counter, _Nz),
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

void gtest::Test_obj::numerical_solution(const int64_t t, const Shift _shift)
{
  Courant_condition_check(_shift);
  (_shift == Shift::shifted) ? field.shifted_field_update(t) : field.field_update(t);
}

double gtest::Test_obj::get_delta_space(void) const
{
  Axis axis = field.get_axis(E, B);
  switch (axis)
  {
  case Axis::Ox: return field.get_dx();
  case Axis::Oy: return field.get_dy();
  case Axis::Oz: return field.get_dz();
  default:
    std::cout << "\nGet delta space: Error! Wrong axis!\n";
    exit(-1);
  }
  return -1.0; // Error code 
}

double gtest::Test_obj::get_global_err(const Component component)
{
  auto get_err = [this](Field::ComputingField& numerical_field, Field::ComputingField& analytical_field) {
    double max_err = 0.0;
    int64_t _Nz = (field.get_Nz() > 1) ? field.get_Nz() : 0; // Учитываем случай 2D
    int64_t _Nz_counter = (_Nz > 1) ? 0 : _Nz - 1; // Это эквивалентно -1 для 2D
    for (int64_t i = 0; i < field.get_Nx(); ++i)
      for (int64_t j = 0; j < field.get_Ny(); ++j)
        for (int64_t k = _Nz_counter; k < _Nz; ++k)
          max_err = std::max(max_err, fabs(numerical_field(i, j, k) - analytical_field(i, j, k)));
    
    return max_err;
  };
  switch (component)
  {
  case Component::Ex:
    return get_err(field.get_Ex(), analytical_field.get_Ex());
  case Component::Ey:
    return get_err(field.get_Ey(), analytical_field.get_Ey());
  case Component::Ez:
    return get_err(field.get_Ez(), analytical_field.get_Ez());
  case Component::Bx:
    return get_err(field.get_Bx(), analytical_field.get_Bx());
  case Component::By:
    return get_err(field.get_By(), analytical_field.get_By());
  case Component::Bz:
    return get_err(field.get_Bz(), analytical_field.get_Bz());
  default:
    std::cout << "\nGet global error: Error! Wrong component!\n";
    exit(-1);
  }
  double error = 0.0; // Error code
  return error;
}

void gtest::Test_obj::print_convergence(Test_obj& other_test)
{
  double this_E_error = this->get_global_err(this->E);
  double this_B_error = this->get_global_err(this->B);

  double other_E_error = other_test.get_global_err(other_test.E);
  double other_B_error = other_test.get_global_err(other_test.B);

  std::cout << "\n=================================\n The 1st (E) error is: " << this_E_error;
  std::cout << "\n The 2nd (E) error is: " << other_E_error;

  std::cout << "\n Difference(E) = " << this_E_error / other_E_error;
  std::cout << "\n Difference(B) = " << this_B_error / other_B_error << "\n=================================" << std::endl;
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
  int64_t i = 0ull;

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

void gtest::Test_obj::set_default_field_or_analytical_default_solution(const Component E, const Component B, const Shift _shift, std::function<void(std::tuple<Axis, int64_t, int64_t>, std::tuple<Axis, int64_t, int64_t>, std::tuple<Axis, int64_t, int64_t>)> loop_function)
{
  Courant_condition_check(_shift);
  double coeff = (_shift == Shift::shifted) ? 0.5 : 0.0;

  Field::ComputingField& Ex = analytical_field.get_Ex();
  Field::ComputingField& Ey = analytical_field.get_Ey();
  Field::ComputingField& Ez = analytical_field.get_Ez();

  Field::ComputingField& Bx = analytical_field.get_Bx();
  Field::ComputingField& By = analytical_field.get_By();
  Field::ComputingField& Bz = analytical_field.get_Bz();


  auto get_E = [this, &E]() -> Field::ComputingField& {
    if (E == Component::Ex) return analytical_field.get_Ex();
    if (E == Component::Ey) return analytical_field.get_Ey();
    if (E == Component::Ez) return analytical_field.get_Ez();
    std::cout << "\nError: Analytical default solution. Wrong E - field!\n";
    exit(-1);
  };
  auto get_B = [this, &B]() -> Field::ComputingField& {
    if (B == Component::Bx) return analytical_field.get_Bx();
    if (B == Component::By) return analytical_field.get_By();
    if (B == Component::Bz) return analytical_field.get_Bz();
    std::cout << "\nError: Analytical default solution. Wrong B - field!\n";
    exit(-1);
  };

  std::pair<double, double> ai_bi{};
  double coordinate = 0.0;
  double delta_coordinate = 0.0;
  double sign = 0.0;

  
  auto helper_set_data = [&] \
    (std::pair<double, double>&&_ai_bi, double delta, double _sign) {
    // (std::pair<double, double>&_ai_bi, double delta, double _sign) {
    ai_bi = _ai_bi;
    coordinate = ai_bi.first;
    delta_coordinate = delta;
    sign = _sign;
  };

  Field::ComputingField& E_field = get_E();
  Field::ComputingField& B_field = get_B();
  
  // Функция синхронизации всех границ для 2D и 3D
  std::function<void(void)> boundary_synchronization{};
  if (this->analytical_field.get_Nz() > 1)
    boundary_synchronization = std::bind(&FDTD::FDTD::boundary_synchronization_3D, &analytical_field);
  else
    boundary_synchronization = std::bind(&FDTD::FDTD::boundary_synchronization, &analytical_field);

  // auto set_computational_data = [this, &helper_set_data, &loop_function, &E, &B]() {
  auto set_computational_data = [&]() {

    int64_t _Nz = (field.get_Nz() > 1) ? field.get_Nz() : 0; // Учитываем случай 2D
    // int64_t _Nz_counter = (_Nz > 1) ? 0 : _Nz - 1; // Это эквивалентно -1 для 2D
    int64_t _Nz_counter = (_Nz > 1) ? 0 : -1; // Это эквивалентно -1 для 2D

    // OX
    if (E == Component::Ey && B == Component::Bz)
    {
      helper_set_data(field.get_ax_bx(), field.get_dx(), 1.0);
      loop_function(std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }
    else if (E == Component::Ez && B == Component::By)
    {
      helper_set_data(field.get_ax_bx(), field.get_dx(), -1.0);
      loop_function(std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }

    // OY
    else if (E == Component::Ez && B == Component::Bx)
    {
      helper_set_data(field.get_ay_by(), field.get_dy(), 1.0);
      loop_function(std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }
    else if (E == Component::Ex && B == Component::Bz)
    {
      helper_set_data(field.get_ay_by(), field.get_dy(), -1.0);
      loop_function(std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }

    // OZ
    // Возможно нужно учесть get_dz() в случае 2D
    else if (E == Component::Ex && B == Component::By)
    {
      helper_set_data(field.get_az_bz(), field.get_dz(), 1.0);
      loop_function(std::make_tuple(Axis::Oz, _Nz_counter, _Nz),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()));
    }
    else if (E == Component::Ey && B == Component::Bx)
    {
      helper_set_data(field.get_az_bz(), field.get_dz(), -1.0);
      loop_function(std::make_tuple(Axis::Oz, _Nz_counter, _Nz),
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

void gtest::Test_obj::set_default_field(const Component E, const Component B, const Shift _shift)
{
  double coeff = (_shift == Shift::shifted) ? 0.5 : 0.0;

  Field::ComputingField& Ex = field.get_Ex();
  Field::ComputingField& Ey = field.get_Ey();
  Field::ComputingField& Ez = field.get_Ez();

  Field::ComputingField& Bx = field.get_Bx();
  Field::ComputingField& By = field.get_By();
  Field::ComputingField& Bz = field.get_Bz();


  auto get_E = [this, &E]() -> Field::ComputingField& {
    if (E == Component::Ex) return field.get_Ex();
    if (E == Component::Ey) return field.get_Ey();
    if (E == Component::Ez) return field.get_Ez();
    std::cout << "\nError: Set default data. Wrong E - field!\n";
    exit(-1);
  };
  auto get_B = [this, &B]() -> Field::ComputingField& {
    if (B == Component::Bx) return field.get_Bx();
    if (B == Component::By) return field.get_By();
    if (B == Component::Bz) return field.get_Bz();
    std::cout << "\nError: Set default data. Wrong B - field!\n";
    exit(-1);
  };

  std::pair<double, double> ai_bi{};
  double coordinate = 0.0;
  double delta_coordinate = 0.0;
  double sign = 0.0;

  int rank = 0;
  int world_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);


  auto helper_set_data = [&] \
    (std::pair<double, double>&&_ai_bi, double delta, double _sign) {
    // (std::pair<double, double>&_ai_bi, double delta, double _sign) {
    ai_bi = _ai_bi;
    
    // coordinate = ai_bi.first;
    delta_coordinate = delta;

    // Вектор значений соответствующих размеров (Nx или Ny или Nz) для каждого процесса
    std::vector<int64_t> all_N(world_size);
    int64_t shift = 0;

    Axis axis = FDTD::FDTD::get_axis(E, B);
    // Опасное место, так как шаг может быть недостаточно сдвинут.
    if (axis == Axis::Ox) {
      int64_t local_Nx = field.get_Nx();
      MPI_Allgather(&local_Nx, 1, MPI_INT64_T, all_N.data(), 1, MPI_INT64_T, MPI_COMM_WORLD);
      for (int i = 0; i < rank; ++i) shift += all_N[i];
      // coordinate = ai_bi.first + delta_coordinate * shift/* + delta_coordinate*/; // Для нескольких процессов
    }
    else if (axis == Axis::Oy) {
      int64_t local_Ny = field.get_Ny();
      MPI_Allgather(&local_Ny, 1, MPI_INT64_T, all_N.data(), 1, MPI_INT64_T, MPI_COMM_WORLD);
      for (int i = 0; i < rank; ++i) shift += all_N[i];
      coordinate = ai_bi.first + delta_coordinate * shift/* + delta_coordinate*/; // Для нескольких процессов
    }
    else if (axis == Axis::Oz) {
      int64_t local_Nz = field.get_Nz();
      MPI_Allgather(&local_Nz, 1, MPI_INT64_T, all_N.data(), 1, MPI_INT64_T, MPI_COMM_WORLD);
      for (int i = 0; i < rank; ++i) shift += all_N[i];
      coordinate = ai_bi.first + delta_coordinate * shift/* + delta_coordinate*/; // Для нескольких процессов
    }

    sign = _sign;
  };

  Field::ComputingField& E_field = get_E();
  Field::ComputingField& B_field = get_B();

  // Функция синхронизации всех границ для 2D и 3D
  std::function<void(void)> boundary_synchronization{};
  if (this->field.get_Nz() > 1)
    boundary_synchronization = std::bind(&FDTD::FDTD::boundary_synchronization_3D, &field);
  else
    boundary_synchronization = std::bind(&FDTD::FDTD::boundary_synchronization, &field);

  auto loop_function = [&](std::tuple<Axis, int64_t, int64_t> axis_1, std::tuple<Axis, int64_t, int64_t> axis_2, std::tuple<Axis, int64_t, int64_t> axis_3) {
    int64_t* i = nullptr, * j = nullptr, * k = nullptr;

    int64_t axis_1_counter = std::get<1>(axis_1);
    int64_t axis_2_counter = std::get<1>(axis_2);
    int64_t axis_3_counter = std::get<1>(axis_3);

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

    for (axis_1_counter = std::get<1>(axis_1); axis_1_counter < std::get<2>(axis_1); ++(axis_1_counter), coordinate += delta_coordinate)
      for (axis_2_counter = std::get<1>(axis_2); axis_2_counter < std::get<2>(axis_2); ++(axis_2_counter))
        for (axis_3_counter = std::get<1>(axis_3); axis_3_counter < std::get<2>(axis_3); ++(axis_3_counter))
        {
          Ex(*i, *j, *k) = Ey(*i, *j, *k) = Ez(*i, *j, *k) = Bx(*i, *j, *k) = By(*i, *j, *k) = Bz(*i, *j, *k) = 0.0;

          E_field(*i, *j, *k) =
            sin(2.0 * PI * (coordinate - ai_bi.first) /
              (ai_bi.second - ai_bi.first));
          B_field(*i, *j, *k) =
            sin(2.0 * PI * (coordinate + delta_coordinate * coeff - ai_bi.first) /
              (ai_bi.second - ai_bi.first));

        }
    boundary_synchronization();
    // if (rank == 0) {
    //     for (int i = 0; i < Bz.get_field().size(); ++i)
    //         Bz.get_field()[i] = i;
    // }
    // else {
    //     for (int i = 0; i < Bz.get_field().size(); ++i) {
    //         Bz.get_field()[i] = 7;
    //     }
    // }
    // for (int proc = 0; proc < 2; ++proc) {
    //     if (rank == proc) {
    //         // for (int i = 0; i < Ex.get_field().size(); ++i) {
    //         //     E_field.get_field()[i] = rank * E_field.get_Ny() * (E_field.get_Nx() + 2) + i;
    //         // }
    //         // Вывод значений для проверки
    //         std::cout << "Process " << rank << " values: ";
    //         for (int i = 0; i < Bz.get_field().size(); ++i) {
    //             std::cout << Bz.get_field()[i] << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     // Синхронизация процессов
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    // boundary_synchronization();
    //     for (int proc = 0; proc < 2; ++proc) {
    //     if (rank == proc) {
    //         // Вывод значений для проверки
    //         std::cout << "Process " << rank << " values: ";
    //         for (int i = 0; i < Bz.get_field().size(); ++i) {
    //             std::cout << Bz.get_field()[i] << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     // Синхронизация процессов
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    // exit(-100);
  };

  auto set_computational_data = [&]() {

    int64_t _Nz = (field.get_Nz() > 1) ? field.get_Nz() : 0; // Учитываем случай 2D
    int64_t _Nz_counter = (_Nz > 1) ? 0 : _Nz - 1;

    // OX
    if (E == Component::Ey && B == Component::Bz)
    {
      helper_set_data(field.get_ax_bx(), field.get_dx(), 1.0);
      loop_function(std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }
    else if (E == Component::Ez && B == Component::By)
    {
      helper_set_data(field.get_ax_bx(), field.get_dx(), -1.0);
      loop_function(std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }

    // OY
    else if (E == Component::Ez && B == Component::Bx)
    {
      helper_set_data(field.get_ay_by(), field.get_dy(), 1.0);
      loop_function(std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }
    else if (E == Component::Ex && B == Component::Bz)
    {
      helper_set_data(field.get_ay_by(), field.get_dy(), -1.0);
      loop_function(std::make_tuple(Axis::Oy, 0, field.get_Ny()),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oz, _Nz_counter, _Nz));
    }

    // OZ
    // Возможно нужно учесть get_dz() в случае 2D
    else if (E == Component::Ex && B == Component::By)
    {
      helper_set_data(field.get_az_bz(), field.get_dz(), 1.0);
      loop_function(std::make_tuple(Axis::Oz, _Nz_counter, _Nz),
        std::make_tuple(Axis::Ox, 0, field.get_Nx()),
        std::make_tuple(Axis::Oy, 0, field.get_Ny()));
    }
    else if (E == Component::Ey && B == Component::Bx)
    {
      helper_set_data(field.get_az_bz(), field.get_dz(), -1.0);
      loop_function(std::make_tuple(Axis::Oz, _Nz_counter, _Nz),
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
