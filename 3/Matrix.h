#pragma once

#include <string>
#include <vector>
#include <iostream>

namespace matrix
{

	class Matrix
	{

	private:
		int nrows;
		int ncols;
		std::vector<std::vector<double>> data;

	public:
		Matrix(std::vector<std::vector<double>> &dat);

		Matrix(int nrow, int ncol);


			virtual void print();
		virtual int getNrows();

		virtual void setNrows(int nrows);

		virtual int getNcols();

		virtual void setNcols(int ncols);

		virtual std::vector<std::vector<double>> getValues();

		virtual void setValues(std::vector<std::vector<double>> &values);

		virtual void setValueAt(int row, int col, double value);

		virtual double getValueAt(int row, int col);

		virtual bool isSquare();

		virtual int size();

		virtual Matrix *multiplyByConstant(double constant);
		virtual Matrix *insertColumnWithValue1();
	};

}
