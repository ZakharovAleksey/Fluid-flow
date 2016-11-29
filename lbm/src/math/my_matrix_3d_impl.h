#include"my_matrix_3d.h"


template<typename T>
Matrix3D<T>::Matrix3D(const int rows, const int colls, const int height) : rows_(rows), colls_(colls), height_(height)
{
	body3d_.resize(height_);
	for (auto & curLayer : body3d_)
		curLayer.Resize(rows_, colls_);

}

template<typename T>
inline void Matrix3D<T>::Set(const int x, const int y, const int z, T const value)
{
	body3d_.at(z)(1, 0) = value;
}


template<typename T1>
std::ostream & operator<<(std::ostream & os, Matrix3D<T1> const & matrix)
{
	using std::endl;

	int curZ = matrix.height_ - 1;
	for (auto curLayer : matrix.body3d_)
	{
		os << "Z = " << curZ-- << endl;
		os << curLayer << std::endl;
	}
	return os;
}