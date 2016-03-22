#pragma once

#include <Eigen/Dense>
#include <string>
#include <json/json.h>

// Matrix handling utilities
class MatrixUtil
{
private:
	MatrixUtil();
	~MatrixUtil();

public:
	// Removes the row at index rowToRemove from the given matrix
	static void matrixRemoveRow(Eigen::MatrixXd & matrix, unsigned int rowToRemove);
	
	// Removes the column at index colToRemove from the given matrix
	static void matrixRemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);

	// Loads a matrix from a file
	// The first line of the file must contain the number of dimensions of the matrix
	static Eigen::MatrixXd loadMatrix(std::string filename, char delimiter = ' ');

	// Saves mat to a file
	// The first line of the file contains the number of dimensions of mat
	static void saveMatrix(const Eigen::MatrixXd & mat, std::string filename, char delimiter = ' ');

	// Converts mat to a JSON representation
	static Json::Value matrixToJSON(const Eigen::MatrixXd & mat);
	
};