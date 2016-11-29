#include"my_matrix_3d.h"


template<typename T>
Matrix3D<T>::Matrix3D(const int rows, const int colls, const int height) : countY_(rows), countX_(colls), countZ_(height)
{
	body3d_.resize(countZ_);
	for (auto & curLayer : body3d_)
		curLayer.Resize(countY_, countX_);

}

template<typename T>
inline long double Matrix3D<T>::GetSum() const
{
	long double sum{ 0.0 };
	std::vector<Matrix2D<T>>::const_iterator curLayerIt = body3d_.begin();

#pragma omp parallel
	{
		while (curLayerIt != body3d_.end())
			#pragma omp atomic
			sum += curLayerIt++->GetSum();
	}
	return sum;
}

template<typename T>
inline std::vector<T> Matrix3D<T>::GetRow(int const y) const
{
	
}

template<typename T>
inline T const Matrix3D<T>::operator()(const int y, const int x, const int z) const
{
	return body3d_.at(z)(y, x);
}

template<typename T>
inline T const Matrix3D<T>::At(const int y, const int x, const int z) const
{
	return body3d_.at(z)(y, x);
}

template<typename T>
inline T & Matrix3D<T>::operator()(const int y, const int x, const int z)
{
	return body3d_.at(z)(y, x);
}

template<typename T>
inline void Matrix3D<T>::Set(const int y, const int x, const int z, T const value)
{
	body3d_.at(z)(y, x) = value;
}

template<typename T>
inline int Matrix3D<T>::CountY() const
{
	return countY_;
}

template<typename T>
inline int Matrix3D<T>::CountX() const
{
	return countX_;
}

template<typename T>
inline int Matrix3D<T>::CountZ() const
{
	return countZ_;
}


template<typename T1>
std::ostream & operator<<(std::ostream & os, Matrix3D<T1> const & matrix)
{
	using std::endl;

	int curZ = matrix.countZ_ - 1;
	for (auto curLayer : matrix.body3d_)
	{
		os << "Z = " << curZ-- << endl;
		os << curLayer << std::endl;
	}
	return os;
}