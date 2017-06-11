#pragma once

#ifndef IMMERSED_BODY_H
#define IMMERSED_BODY_H


#include<memory>
#include<windows.h>

#include <fstream> // file streams
#include <sstream> // string streams



//! Point in 2D space
struct Point
{
	double x_;
	double y_;

	Point() : x_(0.0), y_(0.0) {}
	Point(double x, double y) : x_(x), y_(y) {}
};

//! Single node of immersed in fluid body
struct IBNode
{
	//! Current node position
	Point cur_pos_;
	//! Reference node position
	Point ref_pos_;

	//! Node velocity along x-axis
	double vx_;
	//! Node velocity along y-axis
	double vy_;
	//! Elastic force, acting on the node along x-axis
	double Fx_;
	//! Elastic force, acting on the node along y-axis
	double Fy_;

	IBNode() : cur_pos_(), ref_pos_(), vx_(0.0), vy_(0.0), Fx_(0.0), Fy_(0.0) {}
};

//! Immersed in body fluid
class ImmersedBody
{
public:

	ImmersedBody(int domainX, int domainY, int nodesNumber, Point center, double radius);
	~ImmersedBody() {}


	void CalculateForces();

	void SpreadForces(Matrix2D<double> &fx, Matrix2D<double> & fy);

	void SpreadVelocity(Fluid & fluid);



	void UpdatePosition();

	//! Writes data about boundary of immersed body to *.txt file
	void WriteBodyFormToTxt(const int time);

	//! Writes data about boundary of immersed body to *.vtk file
	void WriteBodyFormToVtk(std::string file_path, const int time);

private:

	//! Calculate strain forces for all nodes of immersed body
	void CalculateStrainForces();

	void CalculateBendingForces();


private:

	int domain_x_;
	int domain_y_;
	//! Number of Lagragian nodes
	int nodes_num;
	//! Stiffness modulus 
	const double stiffness_ = 0.1;
	//! Bending modulus
	const double bending_ = 0.001;

	//! Radius of body
	double radius_;
	//! Center position of body
	Point center_;

	//! Nodes of Lagragian greed
	std::vector<IBNode> body_;

};

#endif // !IMMERSED_BODY_H