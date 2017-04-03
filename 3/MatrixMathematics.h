#pragma once

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:
namespace matrix { class Matrix; }

namespace matrix
{

	class MatrixMathematics
	{

		/// <summary>
		/// This class a matrix utility class and cannot be instantiated.
		/// </summary>
	private:
		MatrixMathematics();
		/// <summary>
		/// Transpose of a matrix - Swap the columns with rows </summary>
		/// <param name="matrix">
		/// @return </param>
	public:
		static Matrix *transpose(Matrix *matrix);

		/// <summary>
		/// Inverse of a matrix - A-1 * A = I where I is the identity matrix
		/// A matrix that have inverse is called non-singular or invertible. If the matrix does not have inverse it is called singular.
		/// For a singular matrix the values of the inverted matrix are either NAN or Infinity
		/// Only square matrices have inverse and the following method will throw exception if the matrix is not square. </summary>
		/// <param name="matrix">
		/// @return </param>
		/// <exception cref="NoSquareException"> </exception>

		static Matrix *inverse(Matrix *matrix);

		/// <summary>
		/// Determinant of a square matrix
		/// The following function find the determinant in a recursively. </summary>
		/// <param name="matrix">
		/// @return </param>
		/// <exception cref="NoSquareException"> </exception>
		static double determinant(Matrix *matrix);

		/// <summary>
		/// Determine the sign; i.e. even numbers have sign + and odds - </summary>
		/// <param name="i">
		/// @return </param>
	private:
		static int changeSign(int i);
		/// <summary>
		/// Creates a submatrix excluding the given row and column </summary>
		/// <param name="matrix"> </param>
		/// <param name="excluding_row"> </param>
		/// <param name="excluding_col">
		/// @return </param>
	public:
		static Matrix *createSubMatrix(Matrix *matrix, int excluding_row, int excluding_col);

		/// <summary>
		/// The cofactor of a matrix </summary>
		/// <param name="matrix">
		/// @return </param>
		/// <exception cref="NoSquareException"> </exception>
		static Matrix *cofactor(Matrix *matrix);

		/// <summary>
		/// Adds two matrices of the same dimension </summary>
		/// <param name="matrix1"> </param>
		/// <param name="matrix2">
		/// @return </param>
		/// <exception cref="IllegalDimensionException"> </exception>
		static Matrix *add(Matrix *matrix1, Matrix *matrix2);

		/// <summary>
		/// subtract two matrices using the above addition method. Matrices should be the same dimension. </summary>
		/// <param name="matrix1"> </param>
		/// <param name="matrix2">
		/// @return </param>
		/// <exception cref="IllegalDimensionException"> </exception>
		static Matrix *subtract(Matrix *matrix1, Matrix *matrix2);

		/// <summary>
		/// Multiply two matrices
		/// </summary>
		/// <param name="matrix1"> </param>
		/// <param name="matrix2"> </param>
		/// <returns>  </returns>
		static Matrix *multiply(Matrix *matrix1, Matrix *matrix2);
	};

}
