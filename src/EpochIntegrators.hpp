/*
 * Authors: Gustavo V. Barroso
 * Created: 12/09/2025
 * Last modified: 18/09/2025
 *
 * Standalone Crank–Nicolson integrators for double and mpfr::mpreal
 */

#pragma once


#include "Matrix.hpp"
#include "Vector.hpp"

#include <mpreal.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/LU>
#include <cmath>
#include <algorithm>


//— Fixed‐step Crank–Nicolson for double ———————————————
inline Vector<double> integrateDoubleCN(
    const Matrix<double>& A,
    const Vector<double>& y0,
    double dt,
    double totalTime)
{
  int steps = std::max(1, int(std::ceil(totalTime / dt)));
  const auto& As = A.eigen();
  int n = As.rows();
  Eigen::SparseMatrix<double> I(n, n);
  I.setIdentity();

  Eigen::SparseMatrix<double> M1 = I - (dt * 0.5) * As;
  Eigen::SparseMatrix<double> M2 = I + (dt * 0.5) * As;
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(M1);

  Eigen::VectorXd y = y0.eigen();
  for (int k = 0; k < steps; ++k)
    y = solver.solve(M2 * y);

  return Vector<double>(std::move(y));
}


//— Fixed‐step Crank–Nicolson for mpfr::mpreal ———————————
inline Vector<mpfr::mpreal> integrateMpfrCN(
    const Matrix<mpfr::mpreal>& A,
    const Vector<mpfr::mpreal>& y0,
    double dt,
    double totalTime)
{
  int steps = std::max(1, int(std::ceil(totalTime / dt)));
  auto Ad = A.eigen();
  int n = Ad.rows();

  // Build identity in full (dense) form for mpfr
  Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> Iden =
    Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic>::Identity(n, n);

  Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> M1 =
    Iden - (dt * mpfr::mpreal(0.5)) * Ad;
  Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> M2 =
    Iden + (dt * mpfr::mpreal(0.5)) * Ad;
  Eigen::PartialPivLU<decltype(M1)> solver(M1);

  auto y = y0.eigen();
  for (int k = 0; k < steps; ++k)
    y = solver.solve(M2 * y);

  return Vector<mpfr::mpreal>(std::move(y));
}


//— Adaptive two‐half‐step Crank–Nicolson for double ———————
inline Vector<double> integrateAdaptiveDoubleCN(
    const Matrix<double>& A,
    const Vector<double>& y0,
    double dt,
    double totalTime,
    double tolerance,
    double dtMin,
    double dtMax)
{
  Eigen::VectorXd y = y0.eigen();
  double t = 0.0, h = dt;

  const auto& As = A.eigen();
  int n = As.rows();
  Eigen::SparseMatrix<double> I(n, n);
  I.setIdentity();

  while (t < totalTime) {
    h = std::clamp(h, dtMin, dtMax);
    if (t + h > totalTime) h = totalTime - t;
    double h2 = h * 0.5;

    // Full‐step
    Eigen::SparseMatrix<double> M1  = I - (h * 0.5) * As;
    Eigen::SparseMatrix<double> M2  = I + (h * 0.5) * As;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solverFull(M1);
    Eigen::VectorXd yFull = solverFull.solve(M2 * y);

    // Two half‐steps
    Eigen::SparseMatrix<double> M1h = I - (h2 * 0.5) * As;
    Eigen::SparseMatrix<double> M2h = I + (h2 * 0.5) * As;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solverHalf(M1h);
    Eigen::VectorXd yHalf = solverHalf.solve(M2h * y);
    yHalf = solverHalf.solve(M2h * yHalf);

    double err = (yFull - yHalf).lpNorm<Eigen::Infinity>();

    if (err <= tolerance) {
      y = yHalf;
      t += h;
      h = std::clamp(h * std::sqrt(tolerance / (err + 1e-16)), dtMin, dtMax);
    } else {
      h = std::max(h * 0.5, dtMin);
    }
  }

  return Vector<double>(std::move(y));
}


//— Adaptive two‐half‐step Crank–Nicolson for mpfr::mpreal ————
inline Vector<mpfr::mpreal> integrateAdaptiveMpfrCN(
    const Matrix<mpfr::mpreal>& A,
    const Vector<mpfr::mpreal>& y0,
    double dt,
    double totalTime,
    double tolerance,
    double dtMin,
    double dtMax)
{
  auto y = y0.eigen();
  double t = 0.0, h = dt;

  auto Ad = A.eigen();
  int n = Ad.rows();
  Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> Iden =
    Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic>::Identity(n, n);

  while (t < totalTime) {
    h = std::clamp(h, dtMin, dtMax);
    if (t + h > totalTime) h = totalTime - t;
    double h2 = h * 0.5;

    // Full‐step
    auto M1 = Iden - (h * mpfr::mpreal(0.5)) * Ad;
    auto M2 = Iden + (h * mpfr::mpreal(0.5)) * Ad;
    Eigen::PartialPivLU<decltype(M1)> solverFull(M1);
    Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> yFull =
      solverFull.solve(M2 * y);

    // Two half‐steps
    auto M1h = Iden - (h2 * mpfr::mpreal(0.5)) * Ad;
    auto M2h = Iden + (h2 * mpfr::mpreal(0.5)) * Ad;
    Eigen::PartialPivLU<decltype(M1h)> solverHalf(M1h);
    Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> yHalf =
      solverHalf.solve(M2h * y);
    yHalf = solverHalf.solve(M2h * yHalf);

    mpfr::mpreal err = (yFull - yHalf).lpNorm<Eigen::Infinity>();

    if (err <= tolerance) {
      y = yHalf;
      t += h;
      h = std::clamp(
        h * std::pow(tolerance / (static_cast<double>(err) + 1e-16), 0.5),
        dtMin, dtMax
      );
    } else {
      h = std::max(h * 0.5, dtMin);
    }
  }

  return Vector<mpfr::mpreal>(std::move(y));
}
