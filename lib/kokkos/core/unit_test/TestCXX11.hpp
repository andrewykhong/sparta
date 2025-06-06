//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include <Kokkos_Core.hpp>

namespace TestCXX11 {

template <class DeviceType>
struct FunctorAddTest {
  using view_type       = Kokkos::View<double**, DeviceType>;
  using execution_space = DeviceType;
  using team_member = typename Kokkos::TeamPolicy<execution_space>::member_type;

  view_type a_, b_;

  FunctorAddTest(view_type& a, view_type& b) : a_(a), b_(b) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const {
    b_(i, 0) = a_(i, 1) + a_(i, 2);
    b_(i, 1) = a_(i, 0) - a_(i, 3);
    b_(i, 2) = a_(i, 4) + a_(i, 0);
    b_(i, 3) = a_(i, 2) - a_(i, 1);
    b_(i, 4) = a_(i, 3) + a_(i, 4);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member& dev) const {
    const int begin = dev.league_rank() * 4;
    const int end   = begin + 4;
    for (int i = begin + dev.team_rank(); i < end; i += dev.team_size()) {
      b_(i, 0) = a_(i, 1) + a_(i, 2);
      b_(i, 1) = a_(i, 0) - a_(i, 3);
      b_(i, 2) = a_(i, 4) + a_(i, 0);
      b_(i, 3) = a_(i, 2) - a_(i, 1);
      b_(i, 4) = a_(i, 3) + a_(i, 4);
    }
  }
};

template <class DeviceType, bool PWRTest>
double AddTestFunctor() {
  using policy_type = Kokkos::TeamPolicy<DeviceType>;

  Kokkos::View<double**, DeviceType> a("A", 100, 5);
  Kokkos::View<double**, DeviceType> b("B", 100, 5);
  typename Kokkos::View<double**, DeviceType>::HostMirror h_a =
      Kokkos::create_mirror_view(a);
  typename Kokkos::View<double**, DeviceType>::HostMirror h_b =
      Kokkos::create_mirror_view(b);

  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 5; j++) {
      h_a(i, j) = 0.1 * i / (1.1 * j + 1.0) + 0.5 * j;
    }
  }
  Kokkos::deep_copy(a, h_a);

  if (PWRTest == false) {
    Kokkos::parallel_for(100, FunctorAddTest<DeviceType>(a, b));
  } else {
    Kokkos::parallel_for(policy_type(25, Kokkos::AUTO),
                         FunctorAddTest<DeviceType>(a, b));
  }
  Kokkos::deep_copy(h_b, b);

  double result = 0;
  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 5; j++) {
      result += h_b(i, j);
    }
  }

  return result;
}

template <class DeviceType, bool PWRTest>
double AddTestLambda() {
  Kokkos::View<double**, DeviceType> a("A", 100, 5);
  Kokkos::View<double**, DeviceType> b("B", 100, 5);
  typename Kokkos::View<double**, DeviceType>::HostMirror h_a =
      Kokkos::create_mirror_view(a);
  typename Kokkos::View<double**, DeviceType>::HostMirror h_b =
      Kokkos::create_mirror_view(b);

  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 5; j++) {
      h_a(i, j) = 0.1 * i / (1.1 * j + 1.0) + 0.5 * j;
    }
  }
  Kokkos::deep_copy(a, h_a);

  if (PWRTest == false) {
    Kokkos::parallel_for(
        100, KOKKOS_LAMBDA(const int& i) {
          b(i, 0) = a(i, 1) + a(i, 2);
          b(i, 1) = a(i, 0) - a(i, 3);
          b(i, 2) = a(i, 4) + a(i, 0);
          b(i, 3) = a(i, 2) - a(i, 1);
          b(i, 4) = a(i, 3) + a(i, 4);
        });
  } else {
    using policy_type = Kokkos::TeamPolicy<DeviceType>;
    using team_member = typename policy_type::member_type;

    policy_type policy(25, Kokkos::AUTO);

    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(const team_member& dev) {
          const unsigned int begin = dev.league_rank() * 4;
          const unsigned int end   = begin + 4;
          for (unsigned int i = begin + dev.team_rank(); i < end;
               i += dev.team_size()) {
            b(i, 0) = a(i, 1) + a(i, 2);
            b(i, 1) = a(i, 0) - a(i, 3);
            b(i, 2) = a(i, 4) + a(i, 0);
            b(i, 3) = a(i, 2) - a(i, 1);
            b(i, 4) = a(i, 3) + a(i, 4);
          }
        });
  }
  Kokkos::deep_copy(h_b, b);

  double result = 0;
  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 5; j++) {
      result += h_b(i, j);
    }
  }

  return result;
}

template <class DeviceType>
struct FunctorReduceTest {
  using view_type       = Kokkos::View<double**, DeviceType>;
  using execution_space = DeviceType;
  using value_type      = double;
  using team_member = typename Kokkos::TeamPolicy<execution_space>::member_type;

  view_type a_;

  FunctorReduceTest(view_type& a) : a_(a) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i, value_type& sum) const {
    sum += a_(i, 1) + a_(i, 2);
    sum += a_(i, 0) - a_(i, 3);
    sum += a_(i, 4) + a_(i, 0);
    sum += a_(i, 2) - a_(i, 1);
    sum += a_(i, 3) + a_(i, 4);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member& dev, value_type& sum) const {
    const int begin = dev.league_rank() * 4;
    const int end   = begin + 4;
    for (int i = begin + dev.team_rank(); i < end; i += dev.team_size()) {
      sum += a_(i, 1) + a_(i, 2);
      sum += a_(i, 0) - a_(i, 3);
      sum += a_(i, 4) + a_(i, 0);
      sum += a_(i, 2) - a_(i, 1);
      sum += a_(i, 3) + a_(i, 4);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& update) const { update = 0.0; }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& update, value_type const& input) const {
    update += input;
  }
};

template <class DeviceType, bool PWRTest>
double ReduceTestFunctor() {
  using policy_type = Kokkos::TeamPolicy<DeviceType>;
  using view_type   = Kokkos::View<double**, DeviceType>;
  using unmanaged_result =
      Kokkos::View<double, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

  view_type a("A", 100, 5);
  typename view_type::HostMirror h_a = Kokkos::create_mirror_view(a);

  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 5; j++) {
      h_a(i, j) = 0.1 * i / (1.1 * j + 1.0) + 0.5 * j;
    }
  }
  Kokkos::deep_copy(a, h_a);

  double result = 0.0;
  if (PWRTest == false) {
    Kokkos::parallel_reduce(100, FunctorReduceTest<DeviceType>(a),
                            unmanaged_result(&result));
  } else {
    Kokkos::parallel_reduce(policy_type(25, Kokkos::AUTO),
                            FunctorReduceTest<DeviceType>(a),
                            unmanaged_result(&result));
  }
  Kokkos::fence();

  return result;
}

template <class DeviceType, bool PWRTest>
double ReduceTestLambda() {
  using policy_type = Kokkos::TeamPolicy<DeviceType>;
  using view_type   = Kokkos::View<double**, DeviceType>;
  using unmanaged_result =
      Kokkos::View<double, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

  view_type a("A", 100, 5);
  typename view_type::HostMirror h_a = Kokkos::create_mirror_view(a);

  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 5; j++) {
      h_a(i, j) = 0.1 * i / (1.1 * j + 1.0) + 0.5 * j;
    }
  }
  Kokkos::deep_copy(a, h_a);

  double result = 0.0;

  if (PWRTest == false) {
    Kokkos::parallel_reduce(
        100,
        KOKKOS_LAMBDA(const int& i, double& sum) {
          sum += a(i, 1) + a(i, 2);
          sum += a(i, 0) - a(i, 3);
          sum += a(i, 4) + a(i, 0);
          sum += a(i, 2) - a(i, 1);
          sum += a(i, 3) + a(i, 4);
        },
        unmanaged_result(&result));
  } else {
    using team_member = typename policy_type::member_type;
    Kokkos::parallel_reduce(
        policy_type(25, Kokkos::AUTO),
        KOKKOS_LAMBDA(const team_member& dev, double& sum) {
          const unsigned int begin = dev.league_rank() * 4;
          const unsigned int end   = begin + 4;
          for (unsigned int i = begin + dev.team_rank(); i < end;
               i += dev.team_size()) {
            sum += a(i, 1) + a(i, 2);
            sum += a(i, 0) - a(i, 3);
            sum += a(i, 4) + a(i, 0);
            sum += a(i, 2) - a(i, 1);
            sum += a(i, 3) + a(i, 4);
          }
        },
        unmanaged_result(&result));
  }
  Kokkos::fence();

  return result;
}

template <class DeviceType>
double TestVariantLambda(int test) {
  switch (test) {
    case 1: return AddTestLambda<DeviceType, false>();
    case 2: return AddTestLambda<DeviceType, true>();
    case 3: return ReduceTestLambda<DeviceType, false>();
    case 4: return ReduceTestLambda<DeviceType, true>();
    default: Kokkos::abort("unreachable");
  }

  return 0;
}

template <class DeviceType>
double TestVariantFunctor(int test) {
  switch (test) {
    case 1: return AddTestFunctor<DeviceType, false>();
    case 2: return AddTestFunctor<DeviceType, true>();
    case 3: return ReduceTestFunctor<DeviceType, false>();
    case 4: return ReduceTestFunctor<DeviceType, true>();
    default: Kokkos::abort("unreachable");
  }

  return 0;
}

template <class DeviceType>
bool Test(int test) {
  double res_functor = TestVariantFunctor<DeviceType>(test);
  double res_lambda  = TestVariantLambda<DeviceType>(test);

  char testnames[5][256] = {" ", "AddTest", "AddTest TeamPolicy", "ReduceTest",
                            "ReduceTest TeamPolicy"};
  bool passed            = true;

  auto a = res_functor;
  auto b = res_lambda;
  // use a tolerant comparison because functors and lambdas vectorize
  // differently https://github.com/trilinos/Trilinos/issues/3233
  auto rel_err = (std::abs(b - a) / std::max(std::abs(a), std::abs(b)));
  auto tol     = 1e-14;
  if (rel_err > tol) {
    passed = false;

    std::cout << "CXX11 ( test = '" << testnames[test]
              << "' FAILED : relative error " << rel_err << " > tolerance "
              << tol << std::endl;
  }

  return passed;
}

}  // namespace TestCXX11

namespace Test {
TEST(TEST_CATEGORY, cxx11) {
  if (std::is_same_v<Kokkos::DefaultExecutionSpace, TEST_EXECSPACE>) {
    ASSERT_TRUE((TestCXX11::Test<TEST_EXECSPACE>(1)));
    ASSERT_TRUE((TestCXX11::Test<TEST_EXECSPACE>(2)));
    ASSERT_TRUE((TestCXX11::Test<TEST_EXECSPACE>(3)));
    ASSERT_TRUE((TestCXX11::Test<TEST_EXECSPACE>(4)));
  }
}

}  // namespace Test
