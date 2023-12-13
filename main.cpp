#include <iostream>
#include "matrix.h"
#include "EigenMatrix.h"


int main() {
    int choice;
    Matrix result1, result2, result_final;

    do {
        std::cout << "Choose an option:" << std::endl;
        std::cout << "1. Enter matrix elements manually" << std::endl;
        std::cout << "2. Generate matrix elements randomly" << std::endl;
        std::cout << "3. Exit" << std::endl;
        std::cout << "Enter your choice: ";
        std::cin >> choice;

        switch (choice) {
            case 1:
                do {
                    result1 = Input_Matrix(false); // Manually enter matrix elements
                    std::cout << "First matrix:" << std::endl;
                    result1.display();

                    std::cout << "\nChoose an option:" << std::endl;
                    std::cout << "1. Enter a different matrix" << std::endl;
                    std::cout << "2. Proceed to operations" << std::endl;
                    std::cout << "Enter your choice: ";
                    std::cin >> choice;
                } while (choice == 1);
                break;

            case 2:
                do {
                    result1 = Input_Matrix(true); // Generate matrix elements randomly
                    std::cout << " matrix:" << std::endl;
                    result1.display();

                    std::cout << "\nChoose an option:" << std::endl;
                    std::cout << "1. Enter a different matrix" << std::endl;
                    std::cout << "2. Proceed to operations" << std::endl;
                    std::cout << "Enter your choice: ";
                    std::cin >> choice;
                } while (choice == 1);
                break;

            case 3:
                std::cout << "Exiting Matrix Calculator. Goodbye!" << std::endl;
                return 0;

            default:
                std::cerr << "Invalid choice. Please enter a valid option." << std::endl;
                break;
        }
    } while (choice != 2);


    do {
        std::cout << "\n------------------ Matrix Calculator ------------------\n";
        std::cout << "| 1. Basic operations                                  |\n";
        std::cout << "| 2. Matrix Inversion and Transformation               |\n";
        std::cout << "| 3. Gaussian Elimination                              |\n";
        std::cout << "| 4. Specialised Matrix Operations                     |\n";
        std::cout << "| 5. Exit                                              |\n";
        std::cout << "--------------------------------------------------------\n";

        std::cout << "Enter your choice: ";
        std::cin >> choice;
        switch (choice) {
            int basicChoice;
    case 1: {
      do{
      std::cout << "\n------ Basic Operations -----\n"; //had to play around to get the correct spacing to make it a rectangular shape
      std::cout << "|\t1. Addition\t\t\t\t|\n";
      std::cout << "|\t2. Subtraction\t\t\t|\n";
      std::cout << "|\t3. Multiplication\t\t|\n";
      std::cout << "|\t4. Determinant\t\t\t|\n";
      std::cout << "|\t5. Back to Main Menu\t|\n";
      std::cout << "-----------------------------\n";
      std::cout << "Enter your choice: ";
      std::cin >> basicChoice;

      switch (basicChoice) {
      case 1:
      case 2:
      case 3: {
        if (choice == 1 || choice == 2) {
          std::cout << "\nChoose an option for the second matrix:" << std::endl;
          std::cout << "1. Enter matrix elements manually" << std::endl;
          std::cout << "2. Generate matrix elements randomly" << std::endl;
          int secondMatrixChoice;
          std::cin >> secondMatrixChoice;

          if (secondMatrixChoice == 1) {
            result2 = Input_Matrix(
                false); // Manually enter matrix elements for the second matrix
          } else if (secondMatrixChoice == 2) {
            result2 = Input_Matrix(true); // Generate matrix elements randomly
                                          // for the second matrix
          } else {
            std::cerr << "Invalid choice for the second matrix. Exiting."
                      << std::endl;
            return 1;
          }

          // Display the second matrix
          std::cout << "Second matrix:" << std::endl;
          result2.display();
        }

        // Perform addition, subtraction, or multiplication
        if (choice == 1) { //simplified the method for basic operations in the main.cpp to make it compact and not make the main.cpp cluttered
          result_final = result1 + result2;
        } else if (choice == 2) {
          result_final = result1 - result2;
        } else {
          result_final = result1 * result2;
        }

        // Display the result of basic operations
        std::cout << "Resulting matrix:" << std::endl;
        result_final.display();
        break;
      }
      case 4: {
        // Determinant of the first matrix
        double det = result1.determinant();
        std::cout << "Determinant of the first matrix: " << det << std::endl;
        break;
      }
      case 5: {
      }
        while (basicChoice != 5)
          ;

        break;
      }
    }
      while (basicChoice != 5)
        ;
      break;
    }
    case 2: {
      // Matrix Inversion and Transformation operations
      int transformChoice; //defining transform choice 
      do {
        std::cout << "\n------- Matrix Inversion and Transformation -------\n";
        std::cout << "|\t1. Transpose\t\t\t\t|\n";
        std::cout << "|\t2. Adjoint Matrix\t\t\t|\n";
        std::cout << "|\t3. Inverse Matrix\t\t\t|\n";
        std::cout << "|\t4. Transformations\t\t\t|\n";
        std::cout << "|\t5. Back to Main Menu\t\t\t|\n";
        std::cout << "-----------------------------------------------\n";
        std::cout << "Enter your choice: ";
        std::cin >> transformChoice; //input of trasnform choice

        switch (transformChoice) {
        case 1: {
          // Transpose operation
          result_final = result1.transpose(); // Transpose operation
          std::cout << "Transposed matrix:" << std::endl;
          result_final.display(); // Displays the result
          break;
        }
        case 2: {
          // Adjoint Matrix operation
          result_final = result1.adjoint(); // Adjoint Matrix operation
          std::cout << "Adjoint matrix:" << std::endl;
          result_final.display(); // Displays the result
          break;
        }
        case 3: {
          // Inverse Matrix operation
          result_final = result1.inverse(); // Inverse Matrix operation
          std::cout << "Inverse matrix:" << std::endl;
          result_final.display(); // Displays the result
          break;
        }
        case 4: {
          // Transformation menu
          std::cout
              << "\n----------- Transformation menu -----------\n";
          std::cout << "|\t1. Scale\t\t|\n";
          std::cout << "|\t2. Translate\t\t|\n";
          std::cout << "|\t3. Rotate\t\t|\n";
          std::cout << "|\t4. Back Inversion Menu\t|\n";
          std::cout << "---------------------------------------\n";
          std::cout << "Enter your choice: ";
          std::cin >> choice; //input of transform choice

          switch (choice) {
          case 1: {
            // Scale Matrix operation
            double scale_factor; //defining scale factor
            std::cout << "Enter the scaling factor: ";
            std::cin >> scale_factor; //input of scale factor
            result_final = result1.scale(scale_factor); // Scale Matrix operation
            std::cout << "Scaled matrix:" << std::endl;
            result_final.display(); // Displays the result
            break;
          }
          case 2: {
            // Translate Matrix operation
            double tx, ty; //defining translation coordinates
            std::cout << "Enter the translation in x-axis: ";
            std::cin >> tx; //input of horizontal translation coordinate
            std::cout << "Enter the translation in y-axis: ";
            std::cin >> ty; //input of vertical translation coordinate
            result_final = result1.translate(tx, ty); //
            std::cout << "Translated matrix:" << std::endl;
            result_final.display(); // Displays the result
            break;
          }
          case 3: {
            // Rotate Matrix operation
            double angle; //defining rotation angle
            std::cout << "Enter the rotation angle (in radians): ";
            std::cin >> angle; //input of rotation angle
            result_final = result1.rotate(angle); // Rotate Matrix operation
            std::cout << "Rotated matrix:" << std::endl;
            result_final.display(); // Displays the result
            break;
          }
          case 4: {
            // Back to Transformations submenu
            break; // Exits the Transformations submenu
          }
          default:
            std::cerr << "Invalid choice. Exiting." << std::endl;
            return 1; // Exits the program with error code 1
          }
          break; // Exits the case
        }
        case 5: {
        }
          while (transformChoice != 5) // Exits the Matrix Inversion and Transformation submenu
            ;

          break; // Exits the case
        }
      } while (transformChoice != 5);
      break; // Exits the case
    }
// Case 3: Gaussian Elimination Menu
case 3: {
    do {
        // Displaying the menu options for Gaussian Elimination
        std::cout << "\n------- Gaussian Elimination -------\n";
        std::cout << "|\t1. Perform Gaussian Elimination\t|\n";
        std::cout << "|\t2. Perform Gauss-Jordan Elimination\t|\n";
        std::cout << "|\t3. Back to main menu\t\t|\n";
        std::cout << "-------------------------------------\n";
        std::cout << "Enter your choice: ";
        std::cin >> choice;

        // Switch case to handle user's choice
        switch (choice) {
            case 1:
                // Gaussian Elimination
                // Call to perform Gaussian Elimination on the matrix
                result1.performGaussianElimination();
                std::cout << "Matrix after Gaussian Elimination:" << std::endl;
                // Display the matrix after Gaussian Elimination
                result1.display();
                break;
            case 2:
                // Gauss-Jordan Elimination
                // Call to perform Gauss-Jordan Elimination on the matrix
                result1.performGaussJordanElimination();
                std::cout << "Matrix after Gauss-Jordan Elimination:" << std::endl;
                // Display the matrix after Gauss-Jordan Elimination
                result1.display();
                break;
            case 3:
                // Option to go back to the main menu
                break;
            default:
                // Error handling for invalid choice
                std::cerr << "Invalid choice. Exiting." << std::endl;
                return 1;
        }
    } while (choice != 3); // Loop to stay in the menu until 'Back to main menu' is selected

    break;
}

// Case 4: Eigenvalue and Eigenvector Calculations
case 4: {
    srand(time(0)); // Seed the random number generator with current time for randomness

    char continueChoice;

    // User input for square matrix size
    size_t size;
    std::cout << "Enter the size of the square matrix: ";
    while (!(std::cin >> size) || size <= 0) {
        // Validating input for matrix size
        std::cout << "Invalid input. Please enter a positive integer: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    int choice;
    // Asking user for matrix initialization method
    std::cout << "Do you want to enter the matrix elements manually (0) or use random numbers (1)? ";
    while (!(std::cin >> choice) || (choice != 0 && choice != 1)) {
        // Validating input for matrix initialization choice
        std::cout << "Invalid input. Please enter 0 or 1: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // Creating and displaying the matrix
    Matrix1 matrix1(size, choice);
    std::cout << "Matrix:" << std::endl;
    // Display the created matrix
    matrix1.display();

    // Preparing for eigenvalue and eigenvector calculations
    Eigen::MatrixXd m = matrix1.toEigenMatrix();
    Eigen::EigenSolver<Eigen::MatrixXd> solver(m);

    do {
        // Menu for specialized matrix operations
    std::cout << "\n----------------- Specialised Matrix Operations -----------------\n";
    std::cout << "| 1. Eigenvalues                                                 |\n";
    std::cout << "| 2. Eigenvectors                                                |\n";
    std::cout << "| 3. Diagonalized matrix (if applicable)                         |\n";
    std::cout << "| 4. Exit Program                                                |\n";
    std::cout << "-----------------------------------------------------------------\n";


        std::cout << "Enter your choice: ";
        std::cin >> choice;

        switch (choice) {
            // Case 1: Compute and display eigenvalues
            case 1: {
                // Symmetric matrix yields real eigenvalues
                if (matrix1.isSymmetric()) {
                    Eigen::VectorXd eigenvalues = solver.eigenvalues().real();
                    std::cout << "Eigenvalues:" << std::endl << eigenvalues << std::endl;
                } else {
                    // Non-symmetric matrix yields complex eigenvalues
                    Eigen::VectorXcd eigenvalues = solver.eigenvalues();
                    std::cout << "Eigenvalues (Complex):" << std::endl << eigenvalues << std::endl;
                    // Display real and imaginary parts separately
                    std::cout << "Real Part:" << std::endl << eigenvalues.real() << std::endl;
                    std::cout << "Imaginary Part:" << std::endl << eigenvalues.imag() << std::endl;
                }
                break;
            }

            // Case 2: Compute and display eigenvectors
            case 2: {
                // Symmetric matrix yields real eigenvectors
                if (matrix1.isSymmetric()) {
                    Eigen::MatrixXd eigenvectors = solver.eigenvectors().real();
                    std::cout << "Eigenvectors:" << std::endl << eigenvectors << std::endl;
                } else {
                    // Non-symmetric matrix yields complex eigenvectors
                    Eigen::MatrixXcd eigenvectors = solver.eigenvectors();
                    std::cout << "Eigenvectors (Complex):" << std::endl << eigenvectors << std::endl;
                    // Display real and imaginary parts separately
                    std::cout << "Real Part:" << std::endl << eigenvectors.real() << std::endl;
                    std::cout << "Imaginary Part:" << std::endl << eigenvectors.imag() << std::endl;
                }
                break;
            }

            // Case 3: Diagonalize the matrix if applicable
            case 3: {
                // Symmetric matrix can be diagonalized directly
                if (matrix1.isSymmetric()) {
                    Eigen::MatrixXd D = solver.eigenvalues().real().asDiagonal();
                    Eigen::MatrixXd V_inverse = solver.eigenvectors().real().inverse();
                    Eigen::MatrixXd diagonalizedMatrix = solver.eigenvectors().real() * D * V_inverse;
                    std::cout << "Diagonalized Matrix:" << std::endl << diagonalizedMatrix << std::endl;
                } else {
                    // Diagonalization process for non-symmetric matrix
                    Eigen::MatrixXcd D = solver.eigenvalues().asDiagonal();
                    Eigen::MatrixXcd V_inverse = solver.eigenvectors().inverse();
                    Eigen::MatrixXcd diagonalizedMatrix = solver.eigenvectors() * D * V_inverse;
                    std::cout << "Diagonalized Matrix (Complex):" << std::endl << diagonalizedMatrix << std::endl;
                    // Display real and imaginary parts separately
                    std::cout << "Real Part:" << std::endl << diagonalizedMatrix.real() << std::endl;
                    std::cout << "Imaginary Part:" << std::endl << diagonalizedMatrix.imag() << std::endl;
                }
                break;
            }

            case 4:{
             // Exit the program
                std::cout << "Program exited." << std::endl;
                break;
            }
            default:{
             // Handling invalid choice
                std::cerr << "Invalid choice. Please choose a valid option." << std::endl;
                break;
            }
        }

        // Asking user if they want to perform another calculation
        std::cout << "Do you want to perform another calculation? (Y/N): ";
        std::cin >> continueChoice;
        while (std::cin.fail() || (continueChoice != 'Y' && continueChoice != 'y' && continueChoice != 'N' && continueChoice != 'n')) {
            // Validating user's choice for continuation
            std::cout << "Invalid input. Please enter Y or N." << std::endl;
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cin >> continueChoice;
        }

    } while (continueChoice == 'Y' or continueChoice == 'y'); // Loop for repeated operations

    break;
}


    
            case 5: {
                // Exit the program
                std::cout << "Exiting Matrix Calculator. Goodbye!" << std::endl;
                return 0;
            }
            default:
                std::cerr << "Invalid choice. Exiting." << std::endl;
                return 1;
        }

        } while (true); // Keep in the menu until the user chooses to exit
    } 