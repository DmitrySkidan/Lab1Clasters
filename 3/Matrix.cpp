#include "Matrix.h"
#include "RectangularVectors.h"

namespace matrix
{

	Matrix::Matrix(std::vector<std::vector<double>> &dat)
	{
		this->data = dat;
		this->nrows = dat.size();
		this->ncols = dat[0].size();
	}

	Matrix::Matrix(int nrow, int ncol)
	{
		this->nrows = nrow;
		this->ncols = ncol;

		data = RectangularVectors::ReturnRectangularDoubleVector(nrow, ncol);
	}

	void Matrix::print()
	{
		for (int i = 0; i < this->nrows; i++)
		{
			for (int j = 0; j < this->ncols; j++)
			{
				std::wcout << data[i][j] << std::wstring(L"\t");
			}
			std::wcout << std::endl;
		}
	}

	int Matrix::getNrows()
	{
		return nrows;
	}

	void Matrix::setNrows(int nrows)
	{
		this->nrows = nrows;
	}

	int Matrix::getNcols()
	{
		return ncols;
	}

	void Matrix::setNcols(int ncols)
	{
		this->ncols = ncols;
	}

	std::vector<std::vector<double>> Matrix::getValues()
	{
		return data;
	}

	void Matrix::setValues(std::vector<std::vector<double>> &values)
	{
		this->data = values;
	}

	void Matrix::setValueAt(int row, int col, double value)
	{
		data[row][col] = value;
	}

	double Matrix::getValueAt(int row, int col)
	{
		return data[row][col];
	}

	bool Matrix::isSquare()
	{
		return nrows == ncols;
	}

	int Matrix::size()
	{
		if (isSquare())
		{
			return nrows;
		}
		return -1;
	}

	Matrix *Matrix::multiplyByConstant(double constant)
	{
		Matrix *mat = new Matrix(nrows, ncols);
		for (int i = 0; i < nrows; i++)
		{
			for (int j = 0; j < ncols; j++)
			{
				mat->setValueAt(i, j, data[i][j] * constant);
			}
		}
		return mat;
	}

	Matrix *Matrix::insertColumnWithValue1()
	{
		Matrix *X_ = new Matrix(this->getNrows(), this->getNcols() + 1);
		for (int i = 0;i < X_->getNrows();i++)
		{
			for (int j = 0;j < X_->getNcols();j++)
			{
				if (j == 0)
				{
					X_->setValueAt(i, j, 1.0);
				}
				else
				{
					X_->setValueAt(i, j, this->getValueAt(i, j - 1));
				}

			}
		}
		return X_;
	}
}
