#include <iostream>
#include <vector>
#include <cstdlib>
#include <eigen-3.4.0/Eigen/Dense>

class Matrix1 {
private:
    std::vector<std::vector<double> > data;  // 2D vector to store matrix data
    size_t size;  // Size of the square matrix

public:
    // Constructor for Matrix1 class
    Matrix1(size_t size, int choice) : size(size), data(size, std::vector<double>(size)) {
        if (choice == 1) {
            // Initialize the matrix with random values if choice is 1
            for (size_t i = 0; i < size; ++i) {
                for (size_t j = 0; j < size; ++j) {
                    data[i][j] = static_cast<double>(rand() % 100);
                }
            }
        } else {
            // Manually input matrix values if choice is not 1
            std::cout << "Enter the elements of the matrix row by row:" << std::endl;
            for (size_t i = 0; i < size; ++i) {
                for (size_t j = 0; j < size; ++j) {
                    // Validate input
                    while (!(std::cin >> data[i][j])) {
                        std::cout << "Invalid input. Please enter a valid number: ";
                        std::cin.clear();
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                }
            }
        }
    }

    // Method to display the matrix
    void display() const {
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                std::cout << data[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // Method to convert the matrix to Eigen::MatrixXd format
    Eigen::MatrixXd toEigenMatrix() const {
        Eigen::MatrixXd m(size, size);
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(size); ++i) {
            for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(size); ++j) {
                m(i, j) = data[i][j];
            }
        }
        return m;
    }

    // Method to check if the matrix is symmetric
    bool isSymmetric() const {
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = i + 1; j < size; ++j) {
                // Check symmetry condition
                if (data[i][j] != data[j][i]) {
                    return false;
                }
            }
        }
        return true;
    }
};
