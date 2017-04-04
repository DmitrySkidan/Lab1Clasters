//----------------------------------------------------------------------------------------
//	Copyright © 2007 - 2017 Tangible Software Solutions Inc.
//	This class can be used by anyone provided that the copyright notice remains intact.
//
//	This class includes methods to convert multidimensional arrays to C++ vectors.
//----------------------------------------------------------------------------------------
class RectangularVectors
{
    using namespace std;
public:
    static vector<vector<double> > ReturnRectangularDoubleVector(int size1, int size2)
    {
        vector<vector<double> > newVector(size1);
        for (int vector1 = 0; vector1 < size1; vector1++)
        {
            newVector[vector1] = vector<double>(size2);
        }

        return newVector;
    }

    static vector<vector<int> > ReturnRectangularIntVector(int size1, int size2)
    {
        vector<vector<int> > newVector(size1);
        for (int vector1 = 0; vector1 < size1; vector1++)
        {
            newVector[vector1] = vector<int>(size2);
        }

        return newVector;
    }
};
