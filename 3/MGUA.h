#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <cmath>

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

namespace matrix {
	class Matrix;
}

namespace mgua
{

	using matrix::Matrix;

	/// <summary>
	/// @author Inna
	/// </summary>
	class MGUA
	{

		/// <param name="args"> the command line arguments </param>
		
	public:/// <exception cref="matrix.NoSquareException"> </exception>
		static void main();

	private:
		static void writeData(const std::wstring &name, std::vector<std::vector<double> > &data);

		static std::vector<std::vector<double> > readData(const std::wstring &name);

		static std::vector<std::vector<double> > getData();

		static std::vector<double> getCriterion(std::vector<std::vector<int> > &models, Matrix *mx, Matrix *my);

		static std::vector<double> getCriterionReg(std::vector<std::vector<int> > &models, Matrix *mx, Matrix *my, Matrix *mxB, Matrix *myB);

		static std::vector<double> getCriterionUnbiesedness(std::vector<std::vector<int> > &models, Matrix *mx, Matrix *my, Matrix *mxB, Matrix *myB);

		static double getSquaresOfY(Matrix *mx, Matrix *my);

		static double getCriterion(std::vector<int> &model, Matrix *mx, Matrix *my);

		static double getCriterionReg(std::vector<int> &model, Matrix *mx, Matrix *my, Matrix *mxB, Matrix *myB);

		static double getCriterionUnbiesedness(std::vector<int> &model, Matrix *mx, Matrix *my, Matrix *mxB, Matrix *myB);

		//Calculation of unbiesedness criterion with given extrapolation. In matrix mxC  -
		// errors may be...
		static double getCriterionUnbiesedness(std::vector<int> &model, Matrix *mx, Matrix *my, Matrix *mxB, Matrix *myB, Matrix *mxC);

		static Matrix *regressParam(Matrix *x, Matrix *y);

		static void print(std::vector<std::vector<double> > &x);

		static void print(std::vector<double> &x);

		static std::vector<std::vector<double> > dataX(int n, int m);

		static std::vector<double> dataY(std::vector<std::vector<double> > &dataX);

		static std::vector<double> dataY(std::vector<std::vector<double> > &dataX, std::vector<double> &b);

		static Matrix *subMatrix(std::vector<int> &model, Matrix *X);

		static std::vector<std::vector<int> > setOfModels(int q);

		static std::vector<int> convertIntToBinary(int q, int r);

		static int pow2(int n);

	};

}
