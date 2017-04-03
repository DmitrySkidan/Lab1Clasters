#define _USE_MATH_DEFINES
#include "MGUA.h"
#include "RectangularVectors.h"
#include "MatrixMathematics.h"
#include "Matrix.h"
#include "mpi.h"
#include <random>
#include <fstream>
#include <cmath>

#define MASTER 0
#define TAG 0

namespace mgua
{
	using matrix::Matrix;
	using matrix::MatrixMathematics;
	using namespace std;
	
	void MGUA::main()
	{
		int numtasks;
		int taskid;
		MPI_Init(0, 0);
		MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
		MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
		int numworkers = numtasks - 1;

		auto t1 = MPI_Wtime();
		if (taskid == MASTER)
		{
			int n;
			int m;

			std::vector<std::vector<double>> data = readData(L"test_data_2003.txt");
				//getData();

			n = data.size();
			m = data[0].size() - 1;

			std::wcout << std::wstring(L"data was red") << std::endl;
			print(data);

			int nA = n - n / 2;
			int nB = n / 2;

			std::vector<std::vector<double>> xA = RectangularVectors::ReturnRectangularDoubleVector(nA, m);
			std::vector<double> yA(nA);

			std::vector<std::vector<double>> xB = RectangularVectors::ReturnRectangularDoubleVector(nB, m);
			std::vector<double> yB(nB);
			for (int j = 0; j < m; j++)
			{
				int ii = 0;
				while (2 * ii < n)
				{
					xA[ii][j] = data[2 * ii][j];
					ii++;
				}
				ii = 0;
				while (2 * ii + 1 < n)
				{
					xB[ii][j] = data[2 * ii + 1][j];
					ii++;
				}
			}
			int ii = 0;
			while (2 * ii < n)
			{
				yA[ii] = data[2 * ii][m];
				ii++;
			}
			ii = 0;
			while (2 * ii + 1 < n)
			{
				yB[ii] = data[2 * ii + 1][m];
				ii++;
			}

			std::vector<std::vector<double>> X = RectangularVectors::ReturnRectangularDoubleVector(nA, m + 1);
			for (int j = 0; j < nA; j++)
			{
				X[j][0] = 1;
				for (int i = 0; i < m; i++)
				{
					X[j][i + 1] = xA[j][i];
				}

			}
			std::vector<std::vector<double>> Y = RectangularVectors::ReturnRectangularDoubleVector(nA, 1);
			for (int j = 0; j < nA; j++)
			{
				Y[j][0] = yA[j];
			}
			std::vector<std::vector<double>> XB = RectangularVectors::ReturnRectangularDoubleVector(nB, m + 1);
			for (int j = 0; j < nB; j++)
			{
				XB[j][0] = 1;
				for (int i = 0; i < m; i++)
				{
					XB[j][i + 1] = xB[j][i];
				}

			}

			std::vector<std::vector<double>> YB = RectangularVectors::ReturnRectangularDoubleVector(nB, 1);
			for (int j = 0; j < nB; j++)
			{
				YB[j][0] = yB[j];
			}
			std::wcout << std::wstring(L" values of matrix XA  :  ") << std::endl;
			print(X);
			std::wcout << std::wstring(L" values of matrix YA  :  ") << std::endl;
			print(Y);

			Matrix *mX = new Matrix(X);
			Matrix *mY = new Matrix(Y);
			Matrix *mXB = new Matrix(XB);
			Matrix *mYB = new Matrix(YB);

			Matrix *mB = regressParam(mX, mY);

			std::wcout << std::wstring(L" values of coefficient b  :  ") << std::endl;
			mB->print();
			Matrix *mYmod;
			mYmod = MatrixMathematics::multiply(mX, mB);
			std::wcout << std::wstring(L" values of matrix mYmod  :  ") << std::endl;
			mYmod->print();

			std::vector<std::vector<int>> models = setOfModels(m);
			std::wcout << std::wstring(L"models with least squares criterion and regularization criterion:") << std::endl;

			std::vector<double> CriterionUnbiesedness = getCriterionUnbiesedness(models, mX, mY, mXB, mYB); // критерій мінімум зсуву

			auto send_per_node = models.size() / numworkers;
			auto offset = 0;
			int* minIds = new int[numworkers];
			for (int i = 1; i <= numworkers; i++)
			{
				
				auto size = i >= send_per_node ? send_per_node : i;
				int** models_array = new int*[size];
				for (size_t i = 0; i < size; i++)
				{
					models_array[i] = new int[models[i].size()];
					for (size_t j = 0; j < models[i].size(); j++)
						models_array[i][j] = models[i].at(j);
				}

				MPI_Send(&size, 1, MPI_INT, i, TAG, MPI_COMM_WORLD);
				MPI_Send(&CriterionUnbiesedness.front() + offset, size, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
				MPI_Send(&offset, 1, MPI_INT, i, TAG, MPI_COMM_WORLD);
				offset += size;
				MPI_Recv(&minIds[i - 1], 1, MPI_INT, i, TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			}

			/*for (auto model : models)
			{
				for (auto mm : model)
				{
					std::wcout << mm << std::wstring(L"\t");
				}
				std::wcout << getCriterion(model, mX, mY) << std::wstring(L"  ") << getCriterionReg(model, mX, mY, mXB, mYB) << std::wstring(L"  \t ") << MGUA::getCriterionUnbiesedness(model, mX, mY, mXB, mYB) << std::endl;
			}*/
			std::vector<int> modelOpt = models[minIds[0]];
			double min = CriterionUnbiesedness[minIds[0]];
			for (int j = 1; j < numworkers; j++)
			{
				if (CriterionUnbiesedness[minIds[j]] < min)
				{
					min = CriterionUnbiesedness[minIds[j]];
					modelOpt = models[minIds[j]];
				}
			}

			std::vector<std::vector<double>> XAB = RectangularVectors::ReturnRectangularDoubleVector(X.size() + XB.size(), m + 1);
			for (int j = 0; j < X.size(); j++)
			{
				for (int i = 0; i <= m; i++)
				{
					XAB[j][i] = X[j][i];
				}
			}
			for (int j = X.size(); j < X.size() + XB.size(); j++)
			{
				for (int i = 0; i <= m; i++)
				{
					XAB[j][i] = XB[j - X.size()][i];
				}
			}

			std::vector<std::vector<double>> YAB = RectangularVectors::ReturnRectangularDoubleVector(X.size() + XB.size(), 1);
			for (int j = 0; j < X.size(); j++)
			{
				YAB[j][0] = yA[j];
			}
			for (int j = X.size(); j < X.size() + XB.size(); j++)
			{
				YAB[j][0] = yB[j - X.size()];
			}

			Matrix tempVar(XAB);
			Matrix *XofModel = subMatrix(modelOpt, &tempVar);
			Matrix tempVar2(YAB);
			mB = regressParam(XofModel, &tempVar2);

			std::wcout << std::wstring(L" values of coefficient mB  :  ") << std::endl;
			mB->print();
			std::vector<std::wstring> nameX(m + 1);
			nameX[0] = L"";
			for (int j = 1; j <= m; j++)
			{
				nameX[j] = std::wstring(L"x") + to_wstring((long double)j);
			}

			std::wcout << std::wstring(L"\n We have such results:");
			std::wcout << std::wstring(L"\n f(x) = ");
			int k = 0;
			for (int j = 0; j <= m; j++)
			{
				if (modelOpt[j] == 1)
				{
					std::wcout << std::abs(mB->getValues()[k][0]) << L" " << nameX[j] << std::endl;
					//wprintf(L"%5.4f %ls ", , nameX[j]);
					k++;
					if (j < mB->getNrows() - 1)
					{
						if (mB->getValues()[k][0] >= 0)
						{
							std::wcout << std::wstring(L" + ");
						}
						else
						{
							std::wcout << std::wstring(L" - ");
						}
					}
				}
			}
			wprintf(L"\n Criterion = %5.10f", min);
			if (min < 0.05)
			{
				std::wcout << std::wstring(L"\n Congratulations! You have a quality model ");
			}

		auto t2 = MPI_Wtime();
		std::cout << "elapsed: " << t2 - t1;
					std::wcout << std::endl;
		}
		else
		{
			int size, offset = 0;
			std::vector<double> criteria;
			MPI_Recv(&size, 1, MPI_INT, MASTER, TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			criteria.resize(size);
			MPI_Recv(&criteria.front(), size, MPI_DOUBLE, MASTER, TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(&offset, 1, MPI_INT, MASTER, TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			//MPI_Recv(&rowLength, 1, MPI_INT, MASTER, TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			
			int min = criteria.at(0);
			int minId = 0;
			for (size_t i = 1; i < criteria.size(); i++)
			{
				if (min > criteria.at(i)) {
					min = criteria.at(i);
					minId = i;
				}
			}

			minId += offset;
			MPI_Send(&minId, 1, MPI_INT, MASTER, TAG, MPI_COMM_WORLD);
		}
		MPI_Finalize();
	}



	std::vector<std::vector<double>> MGUA::getData()
	{
		std::default_random_engine generator;
		std::normal_distribution<double> distribution;

		std::vector<std::vector<double>> data = RectangularVectors::ReturnRectangularDoubleVector(7, 3);
		double delta = 0.0001;
		data[0][0] = -2;
		data[0][1] = -1;
		data[0][2] = 1.0 - 2.0 * data[0][1] + delta * distribution(generator);
		data[1][0] = -1;
		data[1][1] = 0;
		data[1][2] = 1.0 - 2.0 * data[1][1] + delta * distribution(generator);
		data[2][0] = 0;
		data[2][1] = 5;
		data[2][2] = 1.0 - 2.0 * data[2][1] + delta * distribution(generator);
		data[3][0] = 1;
		data[3][1] = 2;
		data[3][2] = 1.0 - 2.0 * data[3][1] + delta * distribution(generator);
		data[4][0] = 2;
		data[4][1] = 3;
		data[4][2] = 1.0 - 2.0 * data[4][1] + delta * distribution(generator);
		data[5][0] = 3;
		data[5][1] = -2;
		data[5][2] = 1.0 - 2.0 * data[5][1] + delta * distribution(generator);
		data[6][0] = 4;
		data[6][1] = 5;
		data[6][2] = 1.0 - 2.0 * data[6][1] + delta * distribution(generator);

		return data;
	}

	std::vector<std::vector<double>> MGUA::readData(const std::wstring &name)
	{
		std::ifstream f(name);

		int rows, cols;
		f >> cols;
		f >> rows;
		cols += 3;

		std::vector<std::vector<double>> result(rows);
		for (size_t i = 0; i < rows; i++)
		{
			int delCos = 365;
			int delSin = 365;
			int delCof = 2;
			double x, y;

			f >> x;
			f >> y;
			result[i] = std::vector<double>(cols + 1);
			for (size_t j = 0; j < cols; j++)
			{
				if (j % 2 == 0) {
					
					result[i][j] = sin(2 * M_PI * x / delSin);
					delSin = delSin / delCof;
					delCof++;
				}
				else
				{
					result[i][j] = cos(2 * M_PI * x / delCos);
					delCos = delCos / delCof;
				}
				result[i][cols] = y;
			}
		}

		f.close();
		return result;
	}

	std::vector<double> MGUA::getCriterion(std::vector<std::vector<int>> &models, Matrix *mx, Matrix *my)
	{
		int n = models.size();

		std::vector<double> c(n);
		for (int j = 0; j < n; j++)
		{
			c[j] = getCriterion(models[j], mx, my);
		}

		return c;

	}

	std::vector<double> MGUA::getCriterionReg(std::vector<std::vector<int>> &models, Matrix *mx, Matrix *my, Matrix *mxB, Matrix *myB)
	{
		int n = models.size();

		std::vector<double> c(n);
		for (int j = 0; j < n; j++)
		{
			c[j] = getCriterionReg(models[j], mx, my, mxB, myB);
		}

		return c;

	}

	std::vector<double> MGUA::getCriterionUnbiesedness(std::vector<std::vector<int>> &models, Matrix *mx, Matrix *my, Matrix *mxB, Matrix *myB)
	{
		int n = models.size();

		std::vector<double> c(n);
		for (int j = 0; j < n; j++)
		{
			c[j] = getCriterionUnbiesedness(models[j], mx, my, mxB, myB);
		}

		return c;

	}

	double MGUA::getSquaresOfY(Matrix *mx, Matrix *my)
	{
		double c = 0;

		for (int j = 0; j < my->getNrows(); j++)
		{
			c += (my->getValueAt(j, 0)) * (my->getValueAt(j, 0));
		}
		return c;
	}

	double MGUA::getCriterion(std::vector<int> &model, Matrix *mx, Matrix *my)
	{
		double c = 0;
		Matrix *XofModel = subMatrix(model, mx);
		Matrix *mB = regressParam(XofModel, my);
		Matrix *mYmod = MatrixMathematics::multiply(XofModel, mB);
		for (int j = 0; j < my->getNrows(); j++)
		{
			c += (my->getValueAt(j, 0) - mYmod->getValueAt(j, 0)) * (my->getValueAt(j, 0) - mYmod->getValueAt(j, 0));
		}
		c = c / getSquaresOfY(mx, my); //нормалізація критерію
		return c;
	}

	double MGUA::getCriterionReg(std::vector<int> &model, Matrix *mx, Matrix *my, Matrix *mxB, Matrix *myB)
	{
		double c = 0;
		Matrix *XofModel = subMatrix(model, mx);
		Matrix *mB = regressParam(XofModel, my);
		Matrix *mYmodOnB = MatrixMathematics::multiply(subMatrix(model, mxB), mB);

		for (int j = 0; j < myB->getNrows(); j++)
		{
			c += (myB->getValueAt(j, 0) - mYmodOnB->getValueAt(j, 0)) * (myB->getValueAt(j, 0) - mYmodOnB->getValueAt(j, 0));
		}
		c = c / getSquaresOfY(mx, my); //нормалізація критерію
		return c;
	}

	double MGUA::getCriterionUnbiesedness(std::vector<int> &model, Matrix *mx, Matrix *my, Matrix *mxB, Matrix *myB)
	{
		double c = 0;
		//коефіцієнт екстраполяції =1
		Matrix *XofModel = subMatrix(model, mx);
		Matrix *mBforA = regressParam(XofModel, my);

		Matrix *XBofModel = subMatrix(model, mxB);
		Matrix *mBforB = regressParam(XBofModel, myB);

		Matrix *mYmodOnAforA = MatrixMathematics::multiply(subMatrix(model, mx), mBforA);
		Matrix *mYmodOnAforB = MatrixMathematics::multiply(subMatrix(model, mxB), mBforA);
		Matrix *mYmodOnBforA = MatrixMathematics::multiply(subMatrix(model, mx), mBforB);
		Matrix *mYmodOnBforB = MatrixMathematics::multiply(subMatrix(model, mxB), mBforB);

		for (int j = 0; j < my->getNrows(); j++)
		{
			c += (mYmodOnAforA->getValueAt(j, 0) - mYmodOnBforA->getValueAt(j, 0)) * (mYmodOnAforA->getValueAt(j, 0) - mYmodOnBforA->getValueAt(j, 0));
		}

		for (int j = 0; j < myB->getNrows(); j++)
		{
			c += (mYmodOnAforB->getValueAt(j, 0) - mYmodOnBforB->getValueAt(j, 0)) * (mYmodOnAforB->getValueAt(j, 0) - mYmodOnBforB->getValueAt(j, 0));
		}
		c = c / getSquaresOfY(mx, my); //нормалізація критерію
		return c;
	}

	double MGUA::getCriterionUnbiesedness(std::vector<int> &model, Matrix *mx, Matrix *my, Matrix *mxB, Matrix *myB, Matrix *mxC)
	{
		double c = 0;
		//коефіцієнт екстраполяції =1
		Matrix *XofModel = subMatrix(model, mx);
		Matrix *mBforA = regressParam(XofModel, my);

		Matrix *XBofModel = subMatrix(model, mxB);
		Matrix *mBforB = regressParam(XBofModel, myB);

		Matrix *mYmodOnAforA = MatrixMathematics::multiply(subMatrix(model, mx), mBforA);
		Matrix *mYmodOnAforB = MatrixMathematics::multiply(subMatrix(model, mxB), mBforA);
		Matrix *mYmodOnBforA = MatrixMathematics::multiply(subMatrix(model, mx), mBforB);
		Matrix *mYmodOnBforB = MatrixMathematics::multiply(subMatrix(model, mxB), mBforB);
		Matrix *mYmodOnAforC = MatrixMathematics::multiply(subMatrix(model, mxC), mBforA);
		Matrix *mYmodOnBforC = MatrixMathematics::multiply(subMatrix(model, mxC), mBforB);

		for (int j = 0; j < mYmodOnAforA->getNrows(); j++)
		{
			c += (mYmodOnAforA->getValueAt(j, 0) - mYmodOnBforA->getValueAt(j, 0)) * (mYmodOnAforA->getValueAt(j, 0) - mYmodOnBforA->getValueAt(j, 0));
		}

		for (int j = 0; j < mYmodOnAforB->getNrows(); j++)
		{
			c += (mYmodOnAforB->getValueAt(j, 0) - mYmodOnBforB->getValueAt(j, 0)) * (mYmodOnAforB->getValueAt(j, 0) - mYmodOnBforB->getValueAt(j, 0));
		}
		for (int j = 0; j < mYmodOnAforC->getNrows(); j++)
		{
			c += (mYmodOnAforC->getValueAt(j, 0) - mYmodOnBforC->getValueAt(j, 0)) * (mYmodOnAforC->getValueAt(j, 0) - mYmodOnBforC->getValueAt(j, 0));
		}
		double alfa = 1 + mxC->getNrows() / (my->getNrows() + myB->getNrows());
		c = c / (getSquaresOfY(mx, my) * alfa); //нормалізація критерію
		return c;
	}

	Matrix *MGUA::regressParam(Matrix *x, Matrix *y)
	{
		Matrix *mB;
		mB = MatrixMathematics::inverse(MatrixMathematics::multiply(MatrixMathematics::transpose(x), x));
		mB = MatrixMathematics::multiply(mB, MatrixMathematics::multiply(MatrixMathematics::transpose(x), y));

		return mB;
	}

	void MGUA::print(std::vector<std::vector<double>> &x)
	{
		for (int i = 0; i < x.size(); i++)
		{
			for (int j = 0; j < x[0].size(); j++)
			{
				std::wcout << x[i][j] << std::wstring(L"\t");
			}
			std::wcout << std::endl;
		}

	}

	void MGUA::print(std::vector<double> &x)
	{
		for (int j = 0; j < x.size(); j++)
		{
			std::wcout << x[j] << std::endl;
		}
		std::wcout << std::endl;
	}

	std::vector<std::vector<double>> MGUA::dataX(int n, int m)
	{
		std::default_random_engine generator;
		std::normal_distribution<double> distribution;

		//ORIGINAL LINE: double[][] x = new double[n][m];
		std::vector<std::vector<double>> x = RectangularVectors::ReturnRectangularDoubleVector(n, m);
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < m; i++)
			{
				x[j][i] = ceil(-10 + 20 * distribution(generator));
			}
		}

		return x;
	}

	std::vector<double> MGUA::dataY(std::vector<std::vector<double>> &dataX)
	{
		std::vector<double> y(dataX.size());


		for (int j = 0; j < dataX.size(); j++)
		{
			y[j] = 1.0;
			for (int i = 0; i < dataX[j].size(); i++)
			{
				y[j] += 1.0 * dataX[j][i];
			}
			std::default_random_engine generator;
			std::normal_distribution<double> distribution;
			y[j] += distribution(generator);
		}

		return y;
	}

	std::vector<double> MGUA::dataY(std::vector<std::vector<double>> &dataX, std::vector<double> &b)
	{
		std::vector<double> y(dataX.size());


		for (int j = 0; j < dataX.size(); j++)
		{
			y[j] = b[0];
			for (int i = 0; i < dataX[j].size(); i++)
			{
				y[j] += b[i + 1] * dataX[j][i];
			}
			std::default_random_engine generator;
			std::normal_distribution<double> distribution;
			y[j] += distribution(generator);

		}

		return y;
	}

	Matrix *MGUA::subMatrix(std::vector<int> &model, Matrix *X)
	{
		int cols = 0;
		for (int j=0;j< model.size();j++)
		{
			if (model[j] == 1)
			{
				cols++;
			}
		}
		Matrix *subX = new Matrix(X->getNrows(), cols);
		int col = 0;
		for (int j = 0; j < model.size(); j++)
		{
			if (model[j] == 1)
			{
				for (int i = 0; i < X->getNrows(); i++)
				{
					subX->setValueAt(i, col, X->getValueAt(i, j));
				}
				col++;
			}

		}
		return subX;
	}

	std::vector<std::vector<int>> MGUA::setOfModels(int q)
	{
		int min = pow2(q) + 1;
		int max = 0;
		for (int i = 0; i <= q; i++)
		{
			max += pow2(i);
		}

		std::vector<std::vector<int>> models = RectangularVectors::ReturnRectangularIntVector(pow2(q) - 1, q);
		int i = 0;
		for (int j = min; j <= max; j++)
		{
			models[i] = convertIntToBinary(q, j);
			i++;
		}

		return models;

	}

	std::vector<int> MGUA::convertIntToBinary(int q, int r)
	{

		std::vector<int> w(q + 1);
		int k = 0;
		for (int j = q; j >= 0; j--)
		{
			if (r < pow2(j))
			{
				w[k] = 0;
				k++;
			}
			else
			{
				w[k] = r / pow2(j);
				r = r % pow2(j);
				k++;
			}
		}

		return w;
	}

	int MGUA::pow2(int n)
	{
		int a = 1;
		if (n < 0)
		{
			return -1;
		}
		if (n == 0)
		{
			return 1;
		}
		for (int j = 1; j <= n; j++)
		{
			a *= 2;
		}
		return a;
	}
}
