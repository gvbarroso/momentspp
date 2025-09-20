#ifndef MATRIXENGINE_HPP
#define MATRIXENGINE_HPP

// MatrixEngine.hpp
// Precision‐aware wrapper around Eigen sparse matrices and dense vectors,
// supporting both double and mpfr::mpreal via std::variant.

#include <variant>
#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <mpreal.h>

#include "Matrix.hpp"
#include "Vector.hpp"
#include <Bpp/Exceptions.h>

class MatrixEngine {
public:
  // Scalar‐typed wrappers
  using MatDouble = Matrix<double>;
  using MatMP     = Matrix<mpfr::mpreal>;
  using VecDouble = Vector<double>;
  using VecMP     = Vector<mpfr::mpreal>;

  // Variants over double vs mpfr::mpreal
  using MatrixVariant      = std::variant<MatDouble, MatMP>;
  using VectorVariant      = std::variant<VecDouble, VecMP>;
  using MatrixEigenVariant = std::variant<
    Eigen::SparseMatrix<double>,
    Eigen::SparseMatrix<mpfr::mpreal>
  >;
  using VectorEigenVariant = std::variant<
    Eigen::VectorXd,
    Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>
  >;

private:
  MatrixVariant matrix_;
  VectorVariant vector_;
  bool use_mpflag_{false};

public:
  // Constructors
  explicit MatrixEngine(bool useMP = false)
    : use_mpflag_(useMP) {}

  explicit MatrixEngine(MatrixVariant m)
    : matrix_(std::move(m)),
      use_mpflag_(std::holds_alternative<MatMP>(matrix_))
  {}

  // Scalar operations
  MatrixEngine& operator*=(double s) {
    std::visit([s](auto& m){ m *= s; }, matrix_);
    std::visit([s](auto& v){ v *= s; }, vector_);
    return *this;
  }

  MatrixEngine& operator/=(double s) {
    if (s == 0.0) throw bpp::Exception("Division by zero");
    std::visit([s](auto& m){ m /= s; }, matrix_);
    std::visit([s](auto& v){ v /= s; }, vector_);
    return *this;
  }

  // Accessors
  MatrixVariant&       getMatrixVariant()       { return matrix_; }
  const MatrixVariant& getMatrixVariant() const { return matrix_; }
  VectorVariant&       getVectorVariant()       { return vector_; }
  const VectorVariant& getVectorVariant() const { return vector_; }
  bool                 useMPReal()       const { return use_mpflag_; }

  // Eigen conversions
  MatrixEigenVariant toEigenMatrixVariant() const {
    MatrixEigenVariant out;
    std::visit([&](auto const& m){ out = m.eigen(); }, matrix_);
    return out;
  }

  VectorEigenVariant toEigenVectorVariant() const {
    VectorEigenVariant out;
    std::visit([&](auto const& v){ out = v.eigen(); }, vector_);
    return out;
  }

  // Setters
  void setMatrix(const MatrixVariant& m) {
    matrix_ = m;
    use_mpflag_ = std::holds_alternative<MatMP>(matrix_);
  }

  void setMatrix(MatrixVariant&& m) {
    matrix_ = std::move(m);
    use_mpflag_ = std::holds_alternative<MatMP>(matrix_);
  }

  void setVector(const VectorVariant& v) {
    vector_ = v;
  }

  void setVector(VectorVariant&& v) {
    vector_ = std::move(v);
  }

  void resetVector() {
    std::visit([](auto& v){ v.setZero(); }, vector_);
  }

  void setMatrixFromEigen(const MatrixEigenVariant& ev) {
    std::visit([&](auto const& me){
      using Scalar = typename std::decay_t<decltype(me)>::Scalar;
      matrix_ = MatrixVariant{ Matrix<Scalar>(me) };
      use_mpflag_ = std::is_same_v<Scalar, mpfr::mpreal>;
    }, ev);
  }

  void setMatrixFromEigen(MatrixEigenVariant&& ev) {
    setMatrixFromEigen(ev);
  }

  void setMatrixFromDense(const Eigen::MatrixXd& dense) {
    std::visit([&](auto& m){
      using Scalar = typename std::decay_t<decltype(m)>::Scalar;
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> tmp(
        dense.rows(), dense.cols()
      );
      for (int i = 0; i < dense.rows(); ++i)
        for (int j = 0; j < dense.cols(); ++j)
          tmp(i,j) = static_cast<Scalar>(dense(i,j));
      m.eigen() = tmp.sparseView();
    }, matrix_);
  }

  void setMatrixFromTriplets(
    const std::vector<Eigen::Triplet<double>>& triplets
  ) {
    std::visit([&](auto& m){
      using Scalar = typename std::decay_t<decltype(m)>::Scalar;
      std::vector<Eigen::Triplet<Scalar>> conv;
      conv.reserve(triplets.size());
      for (auto const& t: triplets)
        conv.emplace_back(
          t.row(), t.col(), static_cast<Scalar>(t.value())
        );
      m.setFromTriplets(conv);
    }, matrix_);
  }

  void setVectorFromEigen(const VectorEigenVariant& ev) {
    std::visit([&](auto const& ve){
      using Scalar = typename std::decay_t<decltype(ve)>::Scalar;
      Vector<Scalar> w(static_cast<size_t>(ve.size()));
      w.eigen() = ve;
      vector_ = VectorVariant{ std::in_place_type<Vector<Scalar>>, std::move(w) };
    }, ev);
  }

  void setVectorFromEigen(VectorEigenVariant&& ev) {
    setVectorFromEigen(ev);
  }

  void setVectorFromMPReal(
    const Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>& mpv
  ) {
    std::visit([&](auto& v){
      using Scalar = typename std::decay_t<decltype(v)>::Scalar;
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> tmp(mpv.size());
      for (int i = 0; i < mpv.size(); ++i)
        tmp(i) = static_cast<Scalar>(mpv(i));
      v.eigen() = tmp;
    }, vector_);
  }

  // Matrix operations
  void addIdentityInPlace() {
    std::visit([](auto& m){
      auto& E = m.eigen();
      for (Eigen::Index i = 0; i < E.rows(); ++i)
        E.coeffRef(i,i) += typename std::decay_t<decltype(m)>::Scalar(1);
    }, matrix_);
  }

  void pruneInPlace() {
    std::visit([](auto& m){
      m.eigen().prune(typename std::decay_t<decltype(m)>::Scalar(0));
    }, matrix_);
  }

  void pruneInPlace(double tol) {
    std::visit([&](auto& m){
      m.eigen().prune(
        static_cast<typename std::decay_t<decltype(m)>::Scalar>(tol)
      );
    }, matrix_);
  }

  void compressInPlace() {
    std::visit([](auto& m){ m.eigen().makeCompressed(); }, matrix_);
  }

  void makeMatrixCompressed() {
    compressInPlace();
  }

  void scaleMatrix(double s) {
    std::visit([&](auto& m){
      m.scale(
        static_cast<typename std::decay_t<decltype(m)>::Scalar>(s)
      );
    }, matrix_);
  }

  void zeroMatrixNegatives() {
    std::visit([](auto& m){ m.zeroNegatives(); }, matrix_);
  }

  void insertMatrixValue(size_t row, size_t col, double val) {
    std::visit([&](auto& m){
      m.insert(
        row, col,
        static_cast<typename std::decay_t<decltype(m)>::Scalar>(val)
      );
    }, matrix_);
  }

  // Matrix–Matrix operations
  void addMatrixInPlace(const MatrixVariant& other) {
    if (matrix_.index() != other.index())
      throw bpp::Exception("MatrixEngine::scalar-type mismatch");
    if (matrix_.index() == 0) {
      std::get<MatDouble>(matrix_).addInPlace(
        std::get<MatDouble>(other)
      );
    } else {
      std::get<MatMP>(matrix_).addInPlace(
        std::get<MatMP>(other)
      );
    }
  }

  MatrixVariant addMatrix(const MatrixVariant& other) const {
    if (matrix_.index() != other.index())
      throw bpp::Exception("MatrixEngine::scalar-type mismatch");
    MatrixVariant out;
    if (matrix_.index() == 0) {
      auto ptr = std::get<MatDouble>(matrix_).add(
        std::get<MatDouble>(other)
      );
      out = *ptr;
    } else {
      auto ptr = std::get<MatMP>(matrix_).add(
        std::get<MatMP>(other)
      );
      out = *ptr;
    }
    return out;
  }

  MatrixVariant multiplyMatrix(const MatrixVariant& other) const {
    if (matrix_.index() != other.index())
      throw bpp::Exception("MatrixEngine::scalar-type mismatch");
    MatrixVariant out;
    if (matrix_.index() == 0) {
      auto ptr = std::get<MatDouble>(matrix_).multiply(
        std::get<MatDouble>(other)
      );
      out = *ptr;
    } else {
      auto ptr = std::get<MatMP>(matrix_).multiply(
        std::get<MatMP>(other)
      );
      out = *ptr;
    }
    return out;
  }

  MatrixVariant fetchScaledMatrix(double s) const {
    MatrixVariant out;
    std::visit([&](auto const& m){
      out = *m.fetchScaled(
        static_cast<typename std::decay_t<decltype(m)>::Scalar>(s)
      );
    }, matrix_);
    return out;
  }

  MatrixVariant identityMatrix() const {
    MatrixVariant out;
    std::visit([&](auto const& m){
      out = *m.identity();
    }, matrix_);
    return out;
  }

  // Matrix–Vector operations
  VectorVariant multiplyMatrixVector() const {
    if (matrix_.index() != vector_.index())
      throw bpp::Exception("MatrixEngine::scalar-type mismatch");
    if (matrix_.index() == 0) {
      auto ptr = std::get<MatDouble>(matrix_).multiply(
        std::get<VecDouble>(vector_)
      );
      return *ptr;
    } else {
      auto ptr = std::get<MatMP>(matrix_).multiply(
        std::get<VecMP>(vector_)
      );
      return *ptr;
    }
  }

  VectorVariant solveSystem() const {
    if (matrix_.index() != vector_.index())
      throw bpp::Exception("MatrixEngine::scalar-type mismatch");
    if (matrix_.index() == 0) {
      auto ptr = std::get<MatDouble>(matrix_).solve(
        std::get<VecDouble>(vector_)
      );
      return *ptr;
    } else {
      auto ptr = std::get<MatMP>(matrix_).solve(
        std::get<VecMP>(vector_)
      );
      return *ptr;
    }
  }

  // Diagnostics
  void logMatrixStats(const std::string& label) const {
    std::visit([&](auto const& m){
      std::cout
        << label << ": "
        << m.rows() << "x" << m.cols()
        << ", nonzeros=" << m.nonZeros()
        << "\n";
    }, matrix_);
  }

  void printMatrix(const std::string& filename) const {
    std::visit([&](auto const& m){ m.print(filename); }, matrix_);
  }

  size_t matrixRows() const {
    return std::visit([](auto const& m){ return m.rows(); }, matrix_);
  }

  size_t matrixCols() const {
    return std::visit([](auto const& m){ return m.cols(); }, matrix_);
  }

  size_t vectorSize() const {
    return std::visit([](auto const& v){ return v.size(); }, vector_);
  }

  void printVector() const {
    std::visit([](auto const& v){ v.print(); }, vector_);
  }

  // Vector cloning
  VectorVariant cloneVector() const {
    VectorVariant out;
    std::visit([&](auto const& v){
      out = *v.clone();
    }, vector_);
    return out;
  }

  VectorVariant cloneVectorWithSize(size_t size) const {
    VectorVariant out;
    std::visit([&](auto const& v){
      out = *v.cloneWithSize(size);
    }, vector_);
    return out;
  }

  // Dense conversion
  template<typename Scalar>
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> getMatrixAsDense() const {
    return std::visit([](auto const& m){
      using MatScalar = typename std::decay_t<decltype(m)>::Scalar;
      static_assert(
        std::is_same_v<MatScalar, Scalar>,
        "Scalar mismatch in getMatrixAsDense"
      );
      return Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(m.eigen());
    }, matrix_);
  }

  // Initialization
  void initialize(size_t rows, size_t cols, size_t vecSize) {
    if (use_mpflag_) {
      matrix_ = MatMP(rows, cols);
      vector_ = VecMP(vecSize);
    } else {
      matrix_ = MatDouble(rows, cols);
      vector_ = VecDouble(vecSize);
    }
  }

  // Deep clone
  std::unique_ptr<MatrixEngine> clone() const {
    auto copy = std::make_unique<MatrixEngine>(use_mpflag_);
    copy->matrix_ = std::visit([](auto const& m){
      return MatrixVariant{ *m.clone() };
    }, matrix_);
    copy->vector_ = std::visit([](auto const& v){
      return VectorVariant{ *v.clone() };
    }, vector_);
    return copy;
  }
};

#endif // MATRIXENGINE_HPP
