#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <stdexcept>
#include <vector>

class Matrix {
private:
  std::vector<std::vector<double> > data;
  size_t rows;
  size_t columns;

public:
  // Constructor
  Matrix(size_t rows = 0, size_t columns = 0)
      : data(rows, std::vector<double>(columns, 0.0)), rows(rows),
        columns(columns) {}

  // Conversion constructor from Matrix to Matrix
  Matrix(const Matrix &matrix)
      : data(matrix.data), rows(matrix.rows), columns(matrix.columns) {}

  // Function to access matrix elements
  double &operator()(size_t i, size_t j) { return data[i][j]; }

  // Const version of operator() for read-only access to matrix elements
  const double &operator()(size_t i, size_t j) const { return data[i][j]; }

  // Function to display the matrix
  void display() const {
    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < columns; ++j) {
        std::cout << data[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }

  // Matrix addition
  Matrix operator+(const Matrix &other) const {
    return add_matrices(*this, other);
  }

  // Matrix subtraction
  Matrix operator-(const Matrix &other) const {
    return subtract_matrices(*this, other);
  }

  // Matrix multiplication
  Matrix operator*(const Matrix &other) const {
    return multiply_matrices(*this, other);
  }

  // function for calculating the determinant of the matrix
  double determinant() const {
    if (rows != columns) { // checks if rows are equal to columns, if not, an
                           // error is printed
      std::cerr << "Error: Determinant is undefined for non-square matrices."
                << std::endl;
      return 0.0;
    }

    if (rows == 1) { // outputs the number itself because the determinant of a
                     // 1x1 matrix is itself
      return data[0][0];
    }

    double det = 0.0;
    int sign = 1;

    for (size_t j = 0; j < columns; ++j) {
      Matrix minor = getMinor(0, j);
      det += sign * data[0][j] * minor.determinant();
      sign = -sign;
    }

    return det;
  }

  // function to get the minor matrix, which we do by excluding a specified row
  // and column
  Matrix getMinor(size_t row, size_t col) const {
    if (rows != columns || row >= rows || col >= columns) {
      throw std::out_of_range("Invalid row or column index for getMinor.");
    }

    Matrix minor(rows - 1,
                 columns - 1); // creates a new matrix, not visible to the user,
                               // which is used to get the matrix

    for (size_t i = 0, r = 0; i < rows; ++i) {
      if (i == row) {
        continue;
      }

      for (size_t j = 0, c = 0; j < columns; ++j) {
        if (j == col) {
          continue;
        }

        minor(r, c) = data[i][j];
        ++c;
      }

      ++r;
    }

    return minor;
  }

  // Other matrix operations can be added as needed

  // This section transposes the matrix
  Matrix transpose() const {
    Matrix result(columns, rows); // Creates a new matrix 'result'

    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < columns; ++j) {
        result(j, i) = data[i][j];
      }
    }

    return result;
  }

  // This section transforms the original matrix into an adjoint matrix
  Matrix adjoint() const {
    Matrix result(
        rows,
        columns); // Creates a new matrix which is a copy of the original matrix

    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < columns; ++j) {
        Matrix minor = getMinor(i, j); // Gets the minor

        // Calculate the cofactor using the formula: (-1)^(i+j) * determinant of
        // the minor and uses pow from cmath library
        double cofactor = std::pow(-1, i + j) * minor.determinant();
        result(i, j) = cofactor;
      }
    }

    return result
        .transpose(); // Finally transposes the cofactor to obtain the adjoint
  }

  // This section transforms the original matrix into an inverse matrix
  Matrix inverse() const {
    double determinant_value = determinant();

    if (determinant_value == 0) {
      std::cerr << "Error: Matrix is singular, cannot find inverse."
                << std::endl;

      return Matrix(rows, columns); //calcultes the determinant and checks wether matrix if singular or not and if it is singular, it outputs an error message.
    }

    Matrix adj = adjoint();
    Matrix result(rows, columns);

    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < columns; ++j) {
        result(i, j) = adj(i, j) / determinant_value; //gets the adjoint matrix and divides it by the determinant to get the inverse matrix
      }
    }

    return result; //displays the inverse matrix
  }

//This section scales the original matrix by a specified scalar value
  Matrix scale(double factor) const {
    Matrix result(rows, columns);

    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < columns; ++j) {
        result(i, j) = factor * data[i][j]; //multiplies the original matrix by the factor inputted by the user
      }
    }

    return result; //displays the scaled matrix
  }

//This section translates the matrix by a specified vector
  Matrix translate(double tx, double ty) const {
    Matrix result(rows, columns);

    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < columns; ++j) {
        if (i == 0) {
          result(i, j) = data[i][j] + tx;
        } else {
          result(i, j) = data[i][j] + ty;
        }
      }
    } //The user is prompted to choose whether they want to translate the matrix horizontally or vertically. The user is then prompted to input the x and y coordinates of the translation. The x and y coordinates are then added to the original matrix to get the translated matrix.

    return result; //This displays the translated matrix
  }

//This section rotates the matrix by a specified angle
  Matrix rotate(double angle) const {
    if (rows != 2 || columns != 2) {
      std::cerr << "Error: Rotation is only defined for 2x2 matrices."
                << std::endl;

      return *this; //If the matrix is not 2x2, an error message is displayed and the original matrix is returned
    }

    Matrix result(rows, columns);

    double cosTheta = std::cos(angle);
    double sinTheta = std::sin(angle);

    result(0, 0) = cosTheta * data[0][0] - sinTheta * data[0][1];
    result(0, 1) = sinTheta * data[0][0] + cosTheta * data[0][1];
    result(1, 0) = cosTheta * data[1][0] - sinTheta * data[1][1];
    result(1, 1) = sinTheta * data[1][0] + cosTheta * data[1][1];

    return result; //If the matrix is a 2x2, each position in the matrix undergoes an operation by multiplying the cosine and sine of the angle. The result is then returned.
  }
// Function to perform Gaussian Elimination on a matrix
void performGaussianElimination() {
    // Loop over each row (excluding the last)
    for (size_t i = 0; i < rows - 1; ++i) {
        // Get the pivot element from the current row
        double pivot = data[i][i];
        // Check if pivot is zero (to avoid division by zero)
        if (pivot == 0.0) {
            // If pivot is zero, need to handle it by swapping rows,
            // skipping the iteration, or using a different strategy
            continue;
        }
        // Eliminate the elements below the pivot
        for (size_t j = i + 1; j < rows; ++j) {
            // Calculate the ratio for row reduction
            double ratio = data[j][i] / pivot;
            // Subtract the scaled row from the current row
            for (size_t k = 0; k < columns; ++k) {
                data[j][k] -= ratio * data[i][k];
            }
        }
    }
}

// Function to perform Gauss-Jordan Elimination on a matrix
void performGaussJordanElimination() {
    // Small number to check for near-zero values (floating point precision)
    const double epsilon = 1e-10;

    // Loop over each row
    for (size_t i = 0; i < rows; ++i) {
        size_t pivotRow = i;
        // Find the row with the maximum element in the current column
        for (size_t j = i + 1; j < rows; ++j) {
            if (std::abs(data[j][i]) > std::abs(data[pivotRow][i])) {
                pivotRow = j;
            }
        }

        // Swap the current row with the row having the maximum element
        std::swap(data[i], data[pivotRow]);

        // Get the pivot element
        double pivot = data[i][i];
        // Check if pivot is near zero to avoid division by a very small number
        if (std::abs(pivot) < epsilon) {
            continue;
        }

        // Normalize the pivot row by dividing by the pivot element
        for (size_t j = 0; j < columns; ++j) {
            data[i][j] /= pivot;
        }

        // Eliminate all other elements in the current column
        for (size_t k = 0; k < rows; ++k) {
            if (k != i) {
                // Calculate the ratio for row reduction
                double ratio = data[k][i];
                // Subtract the scaled row from the current row
                for (size_t j = 0; j < columns; ++j) {
                    data[k][j] -= ratio * data[i][j];
                }
            }
        }
    }
}

  // Second functions for matrix operations
  friend Matrix add_matrices(const Matrix &mat1, const Matrix &mat2);
  friend Matrix subtract_matrices(const Matrix &mat1, const Matrix &mat2);
  friend Matrix multiply_matrices(const Matrix &mat1, const Matrix &mat2);
};

// function for matrix addition
Matrix add_matrices(const Matrix &mat1, const Matrix &mat2) {
  // a check to see if the dimensions of the two matrices are the same
  if (mat1.rows != mat2.rows || mat1.columns != mat2.columns) {
    throw std::invalid_argument(
        "Matrices dimensions do not match for addition");
  }

  // create a new matrix to store the result of the addition
  Matrix result(mat1.rows, mat1.columns);

  // iterate over each element and add corresponding elements from mat1 and mat2
  for (size_t i = 0; i < mat1.rows; ++i) {
    for (size_t j = 0; j < mat1.columns; ++j) {
      result(i, j) = mat1(i, j) + mat2(i, j);
    }
  }

  // return the resulting matrix which can then be outputted
  return result;
}

// function for matrix subtraction (follows the same structure as matrix addition)
Matrix subtract_matrices(const Matrix &mat1, const Matrix &mat2) {
  if (mat1.rows != mat2.rows || mat1.columns != mat2.columns) {
    throw std::invalid_argument(
        "Matrices dimensions do not match for subtraction");
  }

  Matrix result(mat1.rows, mat1.columns);

  for (size_t i = 0; i < mat1.rows; ++i) {
    for (size_t j = 0; j < mat1.columns; ++j) {
      result(i, j) = mat1(i, j) - mat2(i, j);
    }
  }

  return result;
}

// function for matrix multiplication (follows the same structure as matrix addition)
Matrix multiply_matrices(const Matrix &mat1, const Matrix &mat2) {
  if (mat1.columns != mat2.rows) {
    throw std::invalid_argument(
        "Matrices dimensions do not match for multiplication");
  }

  Matrix result(mat1.rows, mat2.columns);

  for (size_t i = 0; i < result.rows; ++i) {
    for (size_t j = 0; j < result.columns; ++j) {
      for (size_t k = 0; k < mat1.columns; ++k) {
        result(i, j) += mat1(i, k) * mat2(k, j);
      }
    }
  }

  return result;
}
/// Function to input or generate a matrix
static Matrix Input_Matrix(bool useRandom = false) {
    int rows, columns;

    // Getting matrix dimensions from user
    std::cout << "Enter the number of rows: ";
    std::cin >> rows;
    std::cout << "Enter the number of columns: ";
    std::cin >> columns;

    Matrix mat(rows, columns);  // Initializing matrix with given dimensions

    if (useRandom) {
        // Filling matrix with random numbers if useRandom is true
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                mat(i, j) = rand() % 100;
            }
        }
    } else {
        // Allowing user to input matrix elements if useRandom is false
        std::cout << "Enter the matrix elements:" << std::endl;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                std::cout << "Enter element at position (" << i + 1 << ", " << j + 1 << "): ";
                std::cin >> mat(i, j);
            }
        }
    }

    return mat;  // Returning the created/filled matrix
}

#endif // MATRIX_H