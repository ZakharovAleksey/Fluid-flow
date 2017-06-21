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
	// œ≈–≈—◊»“¿“‹ ¬ «¿¬»—»ÃŒ—“» Œ“ ”√À¿
	double GetArcLen() override { return 2.0 * M_PI * radius_ / nodes_num; }

private:
	Point center_;
	double radius_;
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


class Microphone
{
public:

	Microphone()
	{
		CreateDataFolder("Data\\ib_lbm_data\\spectrum_data");
	}

	~Microphone() {}

	void PerformMeasurements(const int iter, const std::vector<ImmersedBody*> & im_bodies, const MacroscopicParam<double> & physVal, std::string physValName)
	{
		std::cout << "Perform Measurements\n";

		std::ofstream output_file;
		std::string fileName = "Data\\ib_lbm_data\\spectrum_data\\spectrum_" + physValName + ".txt";
		(iter == 0 ) ? (output_file.open(fileName)) : output_file.open(fileName, std::ios_base::app);

		if (iter == 0)
		{
			output_file << "time\t";
			for (int i = 0; i < im_bodies.size(); ++i)
			{
				if (i != im_bodies.size() - 1)
					output_file << "body" << std::to_string(i) << "\t";
				else
					output_file << "body" << std::to_string(i) << "\n";
			}
		}
		

		if (output_file.is_open())
		{
			std::vector<double> res;

			for (auto imBody : im_bodies)
			{
				double totalValue = 0.0;

				for (auto node : imBody->body_)
				{
					std::vector<std::pair<int, int>> ids;

					double curX = node.cur_pos_.x_;
					double curY = node.cur_pos_.y_;

					// Check if position is integer values (in this case it match with eulerian grid node)
					bool isYinteger = (curY == ceil(curY)) ? true : false;
					bool isXinteger = (curX == ceil(curX)) ? true : false;

					if (!isYinteger && !isXinteger)
					{
						ids.push_back(std::make_pair(floor(curY), floor(curX)));
						ids.push_back(std::make_pair(floor(curY), ceil(curX)));
						ids.push_back(std::make_pair(ceil(curY), floor(curX)));
						ids.push_back(std::make_pair(ceil(curY), ceil(curX)));
					}
					else if (isYinteger && !isXinteger)
					{
						ids.push_back(std::make_pair(curY, floor(curX)));
						ids.push_back(std::make_pair(curY, ceil(curX)));
					}
					else if (!isYinteger && isXinteger)
					{
						ids.push_back(std::make_pair(floor(curY), curX));
						ids.push_back(std::make_pair(ceil(curY), curX));
					}
					else
						ids.push_back(std::make_pair(curY, curX));

					// Calculate average value for all nodes
					for (auto id : ids)
					{
						totalValue += physVal(id.first, id.second);
					}
					totalValue /= ids.size();
				}

				res.push_back(totalValue);
			}
			
			output_file << iter << "\t";
			for (auto val : res)
				output_file << val << "\t";
			output_file << "\n";
			res.clear();
		}
		else
		{
			std::cout << "Could not open file Data\\ib_lbm_data\\spectrum_data\\spectrum.txt for spectrum writing. \n";
		}

		output_file.close();
	}


private:
	void CreateDataFolder(std::string folder_name) const
	{
		// Get path to current directory
		char buffer[MAX_PATH];
		GetModuleFileName(NULL, buffer, MAX_PATH);

		std::string::size_type pos = std::string(buffer).find_last_of("\\/");
		std::string path = std::string(buffer).substr(0, pos);
		path = path.substr(0, path.size() - 6) + "\\" + folder_name;

		char *cstr = new char[path.length() + 1];
		strcpy(cstr, path.c_str());

		// Create folder if not exist yet
		if (GetFileAttributes(cstr) == INVALID_FILE_ATTRIBUTES)
			CreateDirectory(cstr, NULL);
	}

};



#endif // !IMMERSED_BODY_H