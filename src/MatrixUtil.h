#pragma once

#include <eigen3/Eigen/Dense>
#include <string>
#include <json/json.h>

class MatrixUtil
{
private:
	MatrixUtil();
	~MatrixUtil();

public:
	static void matrixRemoveRow(Eigen::MatrixXd & matrix, unsigned int rowToRemove);
	
	static void matrixRemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);

	static Eigen::MatrixXd loadMatrix(std::string filename, char delimiter = ' ');

	static void saveMatrix(const Eigen::MatrixXd & mat, std::string filename, char delimiter = ' ');

	static Json::Value matrixToJSON(const Eigen::MatrixXd & mat);
	
};