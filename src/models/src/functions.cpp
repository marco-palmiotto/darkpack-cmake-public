#pragma once

#include <marty.h>

csl::Expr axisRotation(char axis, const csl::Expr& angle)
{
  csl::Expr rotMatrix{};
  switch (axis)
  {
  case 'x':
    rotMatrix = csl::matrix_s({{CSL_1, CSL_0, CSL_0},
                               {CSL_0, csl::cos_s(angle), -csl::sin_s(angle)},
                               {CSL_0, csl::sin_s(angle), csl::cos_s(angle)}});
    break;
  case 'y':
    rotMatrix = csl::matrix_s({{csl::cos_s(angle), CSL_0, -csl::sin_s(angle)},
                               {CSL_0, CSL_1, CSL_0},
                               {csl::sin_s(angle), CSL_0, csl::cos_s(angle)}});
    break;
  case 'z':
    rotMatrix = csl::matrix_s({{csl::cos_s(angle), -csl::sin_s(angle), CSL_0},
                               {csl::sin_s(angle), csl::cos_s(angle), CSL_0},
                               {CSL_0, CSL_0, CSL_1}});
    break;
  }
  return rotMatrix;
}

csl::Expr mtxProd(const csl::Expr& A, const csl::Expr& B)
{
  const int N = static_cast<int>(csl::Size(A));
  csl::Expr AB = csl::matrix_s(N, N);
  for (size_t i = 0; i != N; ++i)
  {
    for (size_t j = 0; j != N; ++j)
    {
      csl::Expr AB_ij = CSL_0;
      for (size_t k = 0; k != N; ++k)
      {
        AB_ij += A[i][k] * B[k][j];
      }
      AB[i][j] = AB_ij;
    }
  }
  return AB;
}

csl::Expr hermitian(const csl::Expr& A)
{
  const int N = static_cast<int>(csl::Size(A));
  csl::Expr A_T = csl::matrix_s(N, N);
  for (size_t i = 0; i != N; ++i)
  {
    for (size_t j = 0; j != N; ++j)
    {
      A_T[i][j] = csl::GetComplexConjugate(A[j][i]);
    }
  }
  return A_T;
}

csl::Expr mtxProd(const std::vector<csl::Expr>& factors)
{
  if (std::size(factors) == 1)
  {
    return factors[0];
  }
  else
  {
    return mtxProd(factors[0], mtxProd(std::vector<csl::Expr>(factors.begin() + 1, factors.end())));
  }
}

csl::Expr rotation(const csl::Expr& t_x, const csl::Expr& t_y, const csl::Expr& t_z)
{
  csl::Expr V_x = axisRotation('x', t_x);
  csl::Expr V_y = axisRotation('y', t_y);
  csl::Expr V_z = axisRotation('z', t_z);
  return mtxProd({V_z, V_x, V_y});
}

csl::Tensor Q(const std::string& name, const csl::Expr& T, const csl::Space* space, const csl::Expr& t_12,
              const csl::Expr& t_23, const csl::Expr& t_13)
{
  csl::Expr V_mtx = rotation(t_23, t_13, t_12);
  csl::Tensor Q = csl::tensor_s(name, {space, space}, mtxProd({hermitian(V_mtx), T, V_mtx}));
  return Q;
}