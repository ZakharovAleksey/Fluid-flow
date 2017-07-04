#pragma once

#ifndef IMMERSED_BODY_H
#define IMMERSED_BODY_H

#include<memory>
#include<windows.h>
#include<math.h> // for ceil

#include <fstream> // file streams
#include <sstream> // string streams

#include"..\..\modeling_area\fluid.h"
#include"..\..\modeling_area\medium.h"

# define M_PI 3.14159265358979323846  /* pi */
#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code

//! Point in 2D space
struct Point
{
	double y_;
	double x_;

	Point() : x_(0.0), y_(0.0) {}
	Point(double y, double x) : y_(y), x_(x) {}
};

//! Type of immersed body node
enum class IBNodeType
{
	STATIC = 0,	// fixed position in space
	MOVING = 1,	// position could change
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
	//! Type of node (static/moving)
	IBNodeType type_;

	IBNode() : cur_pos_(), ref_pos_(), vx_(0.0), vy_(0.0), Fx_(0.0), Fy_(0.0), type_(IBNodeType::STATIC) {}
};

//! Abstract class of immersed in fluid body
class ImmersedBody
{
	friend class Microphone;
public:

	ImmersedBody();
	ImmersedBody(int domainX, int domainY, int nodesNumber);
	virtual ~ImmersedBody() = 0 {}

	//! Performs calculation of elastic forces, acting between nodes of immersed boundaries
	void CalculateForces();
	//! Performs calculation of additional fluid forces 'fx' and 'fy' based on values of elastic forces
	void SpreadForces(Matrix2D<double> &fx, Matrix2D<double> & fy);
	//! Performs calculation of additional velocities based on values of immersed body nodes velocities
	void SpreadVelocity(Fluid & fluid);
	//! Update immersed body current position
	void UpdatePosition();


	//! Writes data about boundary of immersed body to *.txt file
	void WriteBodyFormToTxt(const int time, const int body_id);
	//! Writes data about boundary of immersed body to *.vtk file
	void WriteBodyFormToVtk(std::string file_path, const int body_id, const int time);

	//! First bad version of RBC-Wall interaction
	friend void Interaction(ImmersedBody* moving_body, ImmersedBody* static_body)
	{
		double rc = 1;
		double kr = 0.005;

		for (int i = 0; i < moving_body->body_.size(); ++i)
		{
			double move_pos_x = moving_body->body_.at(i).cur_pos_.x_;
			double move_pos_y = moving_body->body_.at(i).cur_pos_.y_;

			double static_pos_x = static_body->body_.at(0).cur_pos_.x_;
			double static_pos_y = static_body->body_.at(0).cur_pos_.y_;

			double min_distance = SQ(move_pos_x - static_pos_x) + SQ(move_pos_y - static_pos_y);

			for (int j = 1; j < static_body->body_.size(); ++j)
			{
				double static_pos_x = static_body->body_.at(j).cur_pos_.x_;
				double static_pos_y = static_body->body_.at(j).cur_pos_.y_;

				double cur_distance = SQ(move_pos_x - static_pos_x) + SQ(move_pos_y - static_pos_y);

				if (cur_distance < min_distance)
					min_distance = cur_distance;
			}
			
			min_distance = sqrt(min_distance);

			for (int j = 0; j < static_body->body_.size(); ++j)
			{
				double x = move_pos_x - static_body->body_.at(j).cur_pos_.x_;
				double y = move_pos_y - static_body->body_.at(j).cur_pos_.y_;

				if (SQ(x) + SQ(y) < rc)
				{
					moving_body->body_.at(i).Fx_ += kr * x / abs(pow(min_distance, 3));
					moving_body->body_.at(i).Fy_ += kr * y / abs(pow(min_distance, 3));
				}
			}

		}
	}

protected:

	virtual double GetArcLen() { return 0; };

	//! Performs calculation of strain forces for all nodes of immersed body
	void CalculateStrainForces();
	//! Performs calculation of bending forces for all nodes of immersed body
	void CalculateBendingForces();

protected:

	int domain_x_;
	int domain_y_;
	//! Number of Lagragian nodes
	int nodes_num;
	//! Stiffness modulus 
	const double stiffness_ = 0.1; // 0.1
	//! Bending modulus
	const double bending_ = 0.001; // 0.001

	//! Radius of body
	double radius_;
	//! Center position of body
	Point center_;

	//! Nodes of Lagragian greed
	std::vector<IBNode> body_;

};


// >> Good implementation of immersed bodies


//! Immersed in fluid RBC (Red Blood Cell)
class ImmersedRBC : public ImmersedBody
{
public:
	ImmersedRBC(int domainX, int domainY, int nodesNumber, Point center, double radius);

protected:
	double GetArcLen() override { return 2.0 * M_PI * radius_ / nodes_num; }

private:
	Point center_;
	double radius_;
};

//! A part of circle in range [startAngle, finishAngle] with center point included, or full circle if [0, 2 pi]
class ImmersedCircle : public ImmersedBody
{
public:
	ImmersedCircle(int domainX, int domainY, int nodesNumber, Point center, double radius, double startAngle, double finishAngle);

protected:
	// ÏÅÐÅÑ×ÈÒÀÒÜ Â ÇÀÂÈÑÈÌÎÑÒÈ ÎÒ ÓÃËÀ
	double GetArcLen() override { return 2.0 * M_PI * radius_ / nodes_num; }

private:
	Point center_;
	double radius_;
};


class ImmersedRectangle : public ImmersedBody
{
public:
	ImmersedRectangle(int domainX, int domainY, int nodesNumber, Point rightTop, double width, double height);

protected:
	// ÏÅÐÅÑ×ÈÒÀÒÜ Â ÇÀÂÈÑÈÌÎÑÒÈ ÎÒ ÓÃËÀ
	double GetArcLen() override { return 2.0 * (width_ + height_) / nodes_num; }

private:
	Point rigth_top_;
	double width_;
	double height_;
};

//
////! Immersed in fluid tromb
//class ImmersedBottomTromb : public ImmersedBody
//{
//public:
//	ImmersedBottomTromb(int domainX, int domainY, int nodesNumber, Point center, double radius);
//};

//
////! Immersed in fluid tromb
//class ImmersedTopTromb : public ImmersedBody
//{
//public:
//	ImmersedTopTromb(int domainX, int domainY, int nodesNumber, Point center, double radius);
//};
//
//
////! Immersed in fluid tromb
//class ImmersedTopRect : public ImmersedBody
//{
//public:
//	ImmersedTopRect(int domainX, int domainY, int nodesNumber, Point center, double width, double height);
//};
//
//
////! Immersed in fluid tromb
//class ImmersedBottomRect : public ImmersedBody
//{
//public:
//	ImmersedBottomRect(int domainX, int domainY, int nodesNumber, Point center, double width, double height);
//};
//



#pragma region Microphone

class BloodFlowMicrophone
{
public:
	BloodFlowMicrophone(Point ref_point) : receive_point_(ref_point) 
	{
		assert(ref_point.x_ >= 0);
		assert(ref_point.y_ >= 0);
	}

	~BloodFlowMicrophone() {}

	void AddMeasurePoint(const Point point)
	{
		measure_points_.push_back(point);
	}

	void TakeOffParameters(const int time, const MacroscopicParam<double> & physVal, const std::string physValName)
	{
		std::string outputFileName = "Data\\" + physValName + "_measure.txt";
		std::ofstream outFile;

		// If file open in first time create appropriate header
		if (time == 0)
		{
			outFile.open(outputFileName);
			if (outFile.is_open())
			{
				for (int i = 0; i < measure_points_.size(); ++i)
					outFile << "t" << std::to_string(i) << "\tval" << std::to_string(i) << "\t";
				outFile << "\n";
			}
			else
			{
				outFile.close();
				std::cout << "Eror! Could not open file " << outputFileName << " for writing data!\n";
			}
		}
		else
		{
			// Writes time (lag time is taken into account) and value of measuring physical value
			outFile.open(outputFileName, std::ios_base::app);
		}

		for (auto point : measure_points_)
		{
			int reciveTime = time + int(LagTime(point, receive_point_));
			outFile << reciveTime << "\t" << physVal(int(point.y_), int(point.x_)) << "\t";
		}
		outFile << "\n";
		
		outFile.close();
	}

private:
	//! Calculates lag time for signal coming from point of measure to recieve point (microphone position)
	double LagTime(const Point & first, const Point & second)
	{
		// Speed of sound in LBM is 1/sqrt(3), so below we divide path on speed of sound
		return sqrt(3 * SQ(first.x_ - second.x_) + SQ(first.y_ - second.y_));
	}

	//! Stores all point from which microphone receives incoming signal 
	// with chosen macroscopic physical value in this point
	std::vector<Point> measure_points_;
	//! This pint is the real position of microphone, which receive all incoming
	//! signals from the points of measure
	Point receive_point_;
};




#pragma endregion


#endif // !IMMERSED_BODY_H