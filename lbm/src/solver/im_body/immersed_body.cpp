#include"immersed_body.h"

#define 	M_PI   3.14159265358979323846
#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code


ImmersedBody::ImmersedBody() : domain_x_(0), domain_y_(0), nodes_num(0)
{
	body_.resize(nodes_num, IBNode());
}

ImmersedBody::ImmersedBody(int domainX, int domainY, int nodesNumber) : domain_x_(domainX), domain_y_(domainY), nodes_num(nodesNumber)
{
	body_.resize(nodes_num, IBNode());
}

void ImmersedBody::CalculateForces()
{
	// Set initial values of force equal to zero
	for (auto & node : body_)
	{
		node.Fx_ = 0.0;
		node.Fy_ = 0.0;
	}

	const double arcLen = GetArcLen(); // 2.0 * M_PI * radius_ / nodes_num;

	CalculateStrainForces();
	CalculateBendingForces();

	for (int i = 0; i < nodes_num; ++i)
	{
		if (body_.at(i).type_ == IBNodeType::STATIC)
		{
			body_.at(i).Fx_ = -stiffness_ * (body_.at(i).cur_pos_.x_ - body_.at(i).ref_pos_.x_) * arcLen;
			body_.at(i).Fy_ = -stiffness_ * (body_.at(i).cur_pos_.y_ - body_.at(i).ref_pos_.y_) * arcLen;
		}
	}
}

void ImmersedBody::SpreadForces(Matrix2D<double>& fx, Matrix2D<double>& fy)
{
	//fx.FillWith(0.0);
	//fy.FillWith(0.0);

	for (int i = 0; i < nodes_num; ++i)
	{
		int x_int = (int)(body_.at(i).cur_pos_.x_ - 0.5 + domain_x_) - domain_x_;
		int y_int = (int)(body_.at(i).cur_pos_.y_ + 0.5);

		for (int X = x_int; X <= x_int + 1; ++X) {
			for (int Y = y_int; Y <= y_int + 1; ++Y) {

				// Compute distance between object node and fluid lattice node.

				const double dist_x = body_.at(i).cur_pos_.x_ - 0.5 - X;
				const double dist_y = body_.at(i).cur_pos_.y_ + 0.5 - Y;

				// Compute interpolation weights for x- and y-direction based on the distance.

				const double weight_x = 1 - abs(dist_x);
				const double weight_y = 1 - abs(dist_y);

				// Compute lattice force.

				fx(Y, (X + domain_x_) % domain_x_) += body_.at(i).Fx_ * weight_x * weight_y;
				fy(Y, (X + domain_x_) % domain_x_) += body_.at(i).Fy_ * weight_x * weight_y;
			}
		}

	}
}

void ImmersedBody::SpreadVelocity(Fluid & fluid)
{
	for (int i = 0; i < nodes_num; ++i)
	{
		// Reset node velocity first since '+=' is used.
		body_.at(i).vx_ = 0.0;
		body_.at(i).vy_ = 0.0;

		// Identify the lowest fluid lattice node in interpolation range (see spreading).

		int x_int = (int)(body_.at(i).cur_pos_.x_ - 0.5 + domain_x_) - domain_x_;
		int y_int = (int)(body_.at(i).cur_pos_.y_ + 0.5);

		// Run over all neighboring fluid nodes.
		// In the case of the two-point interpolation, it is 2x2 fluid nodes.

		for (int X = x_int; X <= x_int + 1; ++X) {
			for (int Y = y_int; Y <= y_int + 1; ++Y) {

				// Compute distance between object node and fluid lattice node.

				const double dist_x = body_.at(i).cur_pos_.x_ - 0.5 - X;
				const double dist_y = body_.at(i).cur_pos_.y_ + 0.5 - Y;

				// Compute interpolation weights for x- and y-direction based on the distance.

				const double weight_x = 1 - abs(dist_x);
				const double weight_y = 1 - abs(dist_y);

				// Compute node velocities.

				body_.at(i).vx_ += (fluid.vx_(Y, (X + domain_x_) % domain_x_) * weight_x * weight_y);
				body_.at(i).vy_ += (fluid.vy_(Y, (X + domain_x_) % domain_x_) * weight_x * weight_y);
			}
		}
	}
}

void ImmersedBody::UpdatePosition()
{
	// Reset center position
	//center_.x_ = 0.0;
	//center_.y_ = 0.0;

	// Update node and center positions

	for (int i = 0; i < nodes_num; ++i)
	{
		body_.at(i).cur_pos_.x_ += body_.at(i).vx_;
		body_.at(i).cur_pos_.y_ += body_.at(i).vy_;

		//center_.x_ += body_.at(i).cur_pos_.x_ / nodes_num;
		//center_.y_ += body_.at(i).cur_pos_.y_ / nodes_num;
	}

	/// Check for periodicity along the x-axis

	/*if (center_.x_ < 0)
	{
	center_.x_ += domain_x_;

	for (int n = 0; n < nodes_num; ++n)
	{
	body_.at(n).cur_pos_.x_ += domain_x_;
	}
	}
	else if (center_.x_ >= domain_x_)
	{
	center_.x_ -= domain_x_;

	for (int n = 0; n < nodes_num; ++n)
	{
	body_.at(n).cur_pos_.x_ -= domain_x_;
	}
	}*/


}

void ImmersedBody::WriteBodyFormToTxt(const int time, const int body_id)
{
	std::string file_name = "Data/ib_lbm_data/body_form_txt/body_form" + std::to_string(body_id) + "_t" + std::to_string(time) + ".txt";

	std::ofstream output_file;
	output_file.open(file_name);

	if (output_file.is_open())
	{
		for (auto node : body_)
		{
			output_file << node.cur_pos_.x_ << " " << node.cur_pos_.y_ << std::endl;
		}

		// Add first point to data in file to display closed boundary
		output_file << body_.begin()->cur_pos_.x_ << " " << body_.begin()->cur_pos_.y_;

	}
	else
		std::cout << "Error! Could not open file " << file_name << " to write form of immersed body. \n";

	output_file.close();
}

void ImmersedBody::WriteBodyFormToVtk(std::string file_path, const int body_id, const int time)
{
	std::string file_name = file_path + "\\body_form" + std::to_string(body_id) + "_t" + std::to_string(time) + ".vtk";

	std::ofstream output_file;
	output_file.open(file_name);

	if (output_file.is_open())
	{
		// Write VTK header
		output_file << "# vtk DataFile Version 3.0\n";
		output_file << "particle_state\n";
		output_file << "ASCII\n";
		output_file << "DATASET POLYDATA\n";

		// Write node positions
		output_file << "POINTS " << nodes_num << " float\n";

		for (auto node : body_)
			output_file << node.cur_pos_.x_ << " " << node.cur_pos_.y_ << " 0\n";

		// Write lines between neighboring nodes
		output_file << "LINES " << nodes_num << " " << 3 * nodes_num << "\n";

		for (int i = 0; i < nodes_num; ++i)
			output_file << "2 " << i << " " << (i + 1) % nodes_num << "\n";

		// Write vertices
		output_file << "VERTICES 1 " << nodes_num + 1 << "\n";
		output_file << nodes_num << " ";

		for (int i = 0; i < nodes_num; ++i)
			output_file << i << " ";
	}
	else
		std::cout << "Error! Could not open file " << file_name << " to write form of immersed body. \n";

	output_file.close();
}

void ImmersedBody::CalculateStrainForces()
{
	for (int i = 0; i < nodes_num; ++i)
	{
		if (body_.at(i).type_ == IBNodeType::MOVING)
		{
			const double distance = SQ(body_.at(i).cur_pos_.x_ - body_.at((i + 1) % nodes_num).cur_pos_.x_) + SQ(body_.at(i).cur_pos_.y_ - body_.at((i + 1) % nodes_num).cur_pos_.y_);
			const double distance_ref = SQ(body_.at(i).ref_pos_.x_ - body_.at((i + 1) % nodes_num).ref_pos_.x_) + SQ(body_.at(i).ref_pos_.y_ - body_.at((i + 1) % nodes_num).ref_pos_.y_);

			const double fx = stiffness_ * (distance - distance_ref) * (body_.at(i).cur_pos_.x_ - body_.at((i + 1) % nodes_num).cur_pos_.x_);
			const double fy = stiffness_ * (distance - distance_ref) * (body_.at(i).cur_pos_.y_ - body_.at((i + 1) % nodes_num).cur_pos_.y_);

			// Signs of forces are chosen to satisfy third Newton law
			body_.at(i).Fx_ += -fx;
			body_.at(i).Fy_ += -fy;

			body_.at((i + 1) % nodes_num).Fx_ += fx;
			body_.at((i + 1) % nodes_num).Fy_ += fy;
		}
	}
}

void ImmersedBody::CalculateBendingForces()
{
	for (int i = 0; i < nodes_num; ++i)
	{
		if (body_.at(i).type_ == IBNodeType::MOVING)
		{
			int prevId = (i - 1 + nodes_num) % nodes_num;
			int nextId = (i + 1) % nodes_num;

			const double x_l = body_.at(prevId).cur_pos_.x_;
			const double y_l = body_.at(prevId).cur_pos_.y_;
			const double x_m = body_.at(i).cur_pos_.x_;
			const double y_m = body_.at(i).cur_pos_.y_;
			const double x_r = body_.at(nextId).cur_pos_.x_;
			const double y_r = body_.at(nextId).cur_pos_.y_;

			const double x_l_ref = body_.at(prevId).ref_pos_.x_;
			const double y_l_ref = body_.at(prevId).ref_pos_.y_;
			const double x_m_ref = body_.at(i).ref_pos_.x_;
			const double y_m_ref = body_.at(i).ref_pos_.y_;
			const double x_r_ref = body_.at(nextId).ref_pos_.x_;
			const double y_r_ref = body_.at(nextId).ref_pos_.y_;


			// x-координата вектора, соединяющая l и r
			const double tang_x_ref = x_r_ref - x_l_ref;
			// y-координата вектора, соединяющая l и r
			const double tang_y_ref = y_r_ref - y_l_ref;
			double normal_x_ref;
			double normal_y_ref;

			// Тут просто задем нормаль так чтобы скалярное произведение вектора нормали на 
			// вектор разности l и r были перпендикулярны : а разность модулей  для того чтобы 
			// навпраление было всегда от вне
			if (abs(tang_x_ref) < abs(tang_y_ref)) {
				normal_x_ref = 1;
				normal_y_ref = -tang_x_ref / tang_y_ref;
			}
			else {
				normal_y_ref = 1;
				normal_x_ref = -tang_y_ref / tang_x_ref;
			}

			// То же самое для обычныых
			const double tang_x = x_r - x_l;
			const double tang_y = y_r - y_l;
			double normal_x;
			double normal_y;

			if (abs(tang_x) < abs(tang_y)) {
				normal_x = 1;
				normal_y = -tang_x / tang_y;
			}
			else {
				normal_y = 1;
				normal_x = -tang_y / tang_x;
			}

			// Просто нормализация вектора
			const double normal_length_ref = sqrt(SQ(normal_x_ref) + SQ(normal_y_ref));
			normal_x_ref /= normal_length_ref;
			normal_y_ref /= normal_length_ref;

			if (normal_x_ref * tang_y_ref - normal_y_ref * tang_x_ref > 0) {
				normal_x_ref *= -1;
				normal_y_ref *= -1;
			}

			const double normal_length = sqrt(SQ(normal_x) + SQ(normal_y));
			normal_x /= normal_length;
			normal_y /= normal_length;

			if (normal_x * tang_y - normal_y * tang_x > 0) {
				normal_x *= -1;
				normal_y *= -1;
			}
			// Angle calculations

			// Скалярное произведение векторов деленное на длину векторов
			double angle_ref_cos = (x_l_ref - x_m_ref) * (x_m_ref - x_r_ref) + (y_l_ref - y_m_ref) * (y_m_ref - y_r_ref);
			angle_ref_cos /= (sqrt(SQ(x_l_ref - x_m_ref) + SQ(y_l_ref - y_m_ref)) * sqrt(SQ(x_m_ref - x_r_ref) + SQ(y_m_ref - y_r_ref)));

			// Addition if because of surrounding cos = 1.000000002
			if (angle_ref_cos > 1.0)
				angle_ref_cos = 1.0;
			else if (angle_ref_cos < -1.0)
				angle_ref_cos = -1.0;


			double angle_ref = acos(angle_ref_cos);

			const double convex_x_ref = (x_l_ref + x_r_ref) / 2 - x_m_ref;
			const double convex_y_ref = (y_l_ref + y_r_ref) / 2 - y_m_ref;

			if (convex_x_ref * normal_x_ref + convex_y_ref * normal_y_ref > 0) {
				angle_ref *= -1;
			}

			double angle_cos = (x_l - x_m) * (x_m - x_r) + (y_l - y_m) * (y_m - y_r);
			angle_cos /= (sqrt(SQ(x_l - x_m) + SQ(y_l - y_m)) * sqrt(SQ(x_m - x_r) + SQ(y_m - y_r)));

			// Addition if because of surrounding cos = 1.000000002
			if (angle_cos > 1.0)
				angle_cos = 1.0;
			else if (angle_cos < -1.0)
				angle_cos = -1.0;


			double angle = acos(angle_cos);

			const double convex_x = (x_l + x_r) / 2 - x_m;
			const double convex_y = (y_l + y_r) / 2 - y_m;

			if (convex_x * normal_x + convex_y * normal_y > 0) {
				angle *= -1;
			}
			// force calculations 
			const double force_mag = bending_ * (angle - angle_ref);
			const double length_l = abs(tang_x * (x_m - x_l) + tang_y * (y_m - y_l));
			const double length_r = abs(tang_x * (x_m - x_r) + tang_y * (y_m - y_r));

			body_.at(prevId).Fx_ += normal_x * force_mag * length_l / (length_l + length_r);
			body_.at(prevId).Fy_ += normal_y * force_mag * length_l / (length_l + length_r);

			body_.at(i).Fx_ += -normal_x * force_mag;
			body_.at(i).Fy_ += -normal_y * force_mag;

			body_.at(nextId).Fx_ += normal_x * force_mag * length_r / (length_l + length_r);
			body_.at(nextId).Fy_ += normal_y * force_mag * length_r / (length_l + length_r);
		}
	}
}




ImmersedRBC::ImmersedRBC(int domainX, int domainY, int nodesNumber, Point center, double radius) : ImmersedBody(domainX, domainY, nodesNumber), center_(center), radius_(radius_)
{
	for (int id = 0; id < nodes_num; ++id)
	{
		body_.at(id).type_ = IBNodeType::MOVING;

		// Parametrization of the RBC shape in 2D
		body_.at(id).cur_pos_.y_ = center.y_ + radius * sin(2. * M_PI * (double)id / nodes_num);
		body_.at(id).ref_pos_.y_ = center.y_ + radius * sin(2. * M_PI * (double)id / nodes_num);
		body_.at(id).cur_pos_.x_ = radius * cos(2. * M_PI * (double)id / nodes_num);

		if (body_.at(id).cur_pos_.x_ > 0)
		{
			body_.at(id).cur_pos_.x_ = center.x_ + sqrt(1 - SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius)) * (0.207 + 2.00 * SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius) - 1.12 * SQ(SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius))) * radius / 2;
			body_.at(id).ref_pos_.x_ = center.x_ + sqrt(1 - SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius)) * (0.207 + 2.00 * SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius) - 1.12 * SQ(SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius))) * radius / 2;
		}
		else 
		{
			body_.at(id).cur_pos_.x_ = center.x_ - sqrt(1 - SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius)) * (0.207 + 2.00 * SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius) - 1.12 * SQ(SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius))) * radius / 2;
			body_.at(id).ref_pos_.x_ = center.x_ - sqrt(1 - SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius)) * (0.207 + 2.00 * SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius) - 1.12 * SQ(SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius))) * radius / 2;
		}
	}
}

ImmersedCircle::ImmersedCircle(int domainX, int domainY, int nodesNumber, Point center, double radius, double startAngle, double finishAngle) : ImmersedBody(domainX, domainY, nodesNumber), center_(center), radius_(radius)
{
	assert(startAngle < finishAngle);
	assert(finishAngle - startAngle <= 2.0 * M_PI);

	// True if user input to plot fill circle, or interval [0, 2 * pi]
	bool isFullCircle = finishAngle - startAngle == 2.0 * M_PI;

	// Obtain one step rotation angle and number of points on circle (in not full - 1 node, because it need to center node add)
	double angleStep = (!isFullCircle) ? (finishAngle - startAngle) / (nodes_num - 2) : (finishAngle - startAngle) / (nodes_num - 1);
	double circleNumb = (!isFullCircle) ? nodes_num - 1 : nodes_num;
	
	// Fill circle
	for (int id = 0; id < circleNumb; ++id)
	{
		body_.at(id).type_ = IBNodeType::MOVING;

		// Parametrization of the circle shape in 2D
		body_.at(id).cur_pos_.x_ = center_.x_ + radius_ * cos(startAngle + (double)id * angleStep);
		body_.at(id).ref_pos_.x_ = body_.at(id).cur_pos_.x_;

		body_.at(id).cur_pos_.y_ = center_.y_ + radius_ * sin(startAngle + (double)id * angleStep);
		body_.at(id).ref_pos_.y_ = body_.at(id).cur_pos_.y_;
	}

	// Add additional center node if user input not full circle
	if (!isFullCircle)
	{
		body_.at(nodesNumber - 1).type_ = IBNodeType::STATIC;

		body_.at(nodesNumber - 1).cur_pos_.x_ = center.x_;
		body_.at(nodesNumber - 1).ref_pos_.x_ = center.x_;

		body_.at(nodesNumber - 1).cur_pos_.y_ = center.y_;
		body_.at(nodesNumber - 1).ref_pos_.y_ = center.y_;
	}

	body_.at(0).type_ = IBNodeType::STATIC;
	body_.at(nodes_num - 2).type_ = IBNodeType::STATIC;
	body_.at(nodes_num - 1).type_ = IBNodeType::STATIC;
}


ImmersedRectangle::ImmersedRectangle(int domainX, int domainY, int nodesNumber, Point rightTop, double width, double height) : ImmersedBody(domainX, domainY, nodesNumber), rigth_top_(rightTop), width_(width), height_(height)
{
	int pointToSide = int(nodes_num / 4.0);
	const double xStep = width_ / pointToSide;
	const double yStep = height_ / pointToSide;

	for (int y = 0; y < pointToSide; ++y)
	{
		body_.at(y).cur_pos_.x_ = rigth_top_.x_;
		body_.at(y).cur_pos_.y_ = rigth_top_.y_ - y * yStep;

		body_.at(y).ref_pos_.x_ = rigth_top_.x_;
		body_.at(y).ref_pos_.y_ = rigth_top_.y_ - y * yStep;
	}

	for (int x = 0; x < pointToSide; ++x)
	{
		body_.at(pointToSide + x).cur_pos_.x_ = rigth_top_.x_ + x * xStep;
		body_.at(pointToSide + x).cur_pos_.y_ = rigth_top_.y_ - height_;

		body_.at(pointToSide + x).ref_pos_.x_ = rigth_top_.x_ + x * xStep;
		body_.at(pointToSide + x).ref_pos_.y_ = rigth_top_.y_ - height_;
	}

	for (int y = 0; y < pointToSide; ++y)
	{
		body_.at(2 * pointToSide + y).cur_pos_.x_ = rigth_top_.x_ + width_;
		body_.at(2 * pointToSide + y).cur_pos_.y_ = rigth_top_.y_ - height_ + y * yStep;

		body_.at(2 * pointToSide + y).ref_pos_.x_ = rigth_top_.x_ + width_;
		body_.at(2 * pointToSide + y).ref_pos_.y_ = rigth_top_.y_ - height_ + y * yStep;
	}

	for (int x = 0; x < pointToSide; ++x)
	{
		body_.at(3 * pointToSide + x).cur_pos_.x_ = rigth_top_.x_ + width_ - x * xStep;
		body_.at(3 * pointToSide + x).cur_pos_.y_ = rigth_top_.y_;

		body_.at(3 * pointToSide + x).ref_pos_.x_ = rigth_top_.x_ + width_ - x * xStep;
		body_.at(3 * pointToSide + x).ref_pos_.y_ = rigth_top_.y_;
	}

	for (auto & node : body_)
		node.type_ = IBNodeType::STATIC;
}

//ImmersedBottomTromb::ImmersedBottomTromb(int domainX, int domainY, int nodesNumber, Point center, double radius) : ImmersedBody(domainX, domainY, nodesNumber, center, radius)
//{
//	// Parametrization of the circle shape in 2D (half of circle)
//	double h = radius_ * 4.0 / nodes_num;
//
//	for (int id = 0; id < nodes_num / 2; ++id)
//	{
//		body_.at(id).type_ = IBNodeType::STATIC;
//
//		body_.at(id).cur_pos_.x_ = center_.x_ + radius_ - id * h;
//		body_.at(id).ref_pos_.x_ = body_.at(id).cur_pos_.x_;
//
//		body_.at(id).cur_pos_.y_ = center_.y_;
//		body_.at(id).ref_pos_.y_ = body_.at(id).cur_pos_.y_;
//	}
//
//	for (int id = nodes_num / 2; id < nodes_num; ++id)
//	{
//		body_.at(id).type_ = IBNodeType::MOVING; // Moveing
//
//		body_.at(id).cur_pos_.x_ = center_.x_ + radius_ * cos(2.0 * M_PI * (double)id / nodes_num);
//		body_.at(id).ref_pos_.x_ = body_.at(id).cur_pos_.x_;
//
//		body_.at(id).cur_pos_.y_ = center_.y_ - radius_ * sin(2.0 * M_PI * (double)id / nodes_num);
//		body_.at(id).ref_pos_.y_ = body_.at(id).cur_pos_.y_;
//	}
//}
//
//ImmersedTopTromb::ImmersedTopTromb(int domainX, int domainY, int nodesNumber, Point center, double radius) : ImmersedBody(domainX, domainY, nodesNumber, center, radius)
//{
//	// Parametrization of the circle shape in 2D (half of circle)
//	double h = radius_ * 4.0 / nodes_num;
//
//	for (int id = 0; id < nodes_num / 2; ++id)
//	{
//		body_.at(id).type_ = IBNodeType::STATIC;
//
//		body_.at(id).cur_pos_.x_ = center_.x_ - radius_ + id * h;
//		body_.at(id).ref_pos_.x_ = body_.at(id).cur_pos_.x_;
//
//		body_.at(id).cur_pos_.y_ = center_.y_;
//		body_.at(id).ref_pos_.y_ = body_.at(id).cur_pos_.y_;
//	}
//
//	for (int id = nodes_num / 2; id < nodes_num; ++id)
//	{
//		body_.at(id).type_ = IBNodeType::MOVING; // Moving
//
//		body_.at(id).cur_pos_.x_ = center_.x_ - radius_ * cos(2.0 * M_PI * (double)id / nodes_num);
//		body_.at(id).ref_pos_.x_ = body_.at(id).cur_pos_.x_;
//
//		body_.at(id).cur_pos_.y_ = center_.y_ + radius_ * sin(2.0 * M_PI * (double)id / nodes_num);
//		body_.at(id).ref_pos_.y_ = body_.at(id).cur_pos_.y_;
//	}
//
//}
//
//ImmersedTopRect::ImmersedTopRect(int domainX, int domainY, int nodesNumber, Point center, double  width, double height) : ImmersedBody(domainX, domainY, nodesNumber, center, (width + height) / M_PI)
//{
//	double x_start = center.x_ - width / 2;
//	double y_start = center.y_ + height / 2;
//
//	double x_step = width / nodesNumber * 4.0;
//	double y_step = height / nodesNumber * 4.0;
//
//	// TOP
//	int j = 0;
//
//	for (int i = 0; i < nodesNumber / 4; ++i)
//	{
//		body_.at(i).type_ = IBNodeType::STATIC;
//
//		body_.at(i).cur_pos_.x_ = x_start + j * x_step;
//		body_.at(i).ref_pos_.x_ = x_start + j * x_step;
//
//		body_.at(i).cur_pos_.y_ = y_start;
//		body_.at(i).ref_pos_.y_ = y_start;
//		j++;
//	}
//
//	// RIGHT
//	j = 0;
//
//	for (int i = nodesNumber / 4; i < nodesNumber / 2; ++i)
//	{
//		body_.at(i).type_ = IBNodeType::MOVING;
//
//		body_.at(i).cur_pos_.x_ = x_start + width;
//		body_.at(i).ref_pos_.x_ = x_start + width;
//
//		body_.at(i).cur_pos_.y_ = y_start - j * y_step;
//		body_.at(i).ref_pos_.y_ = y_start - j * y_step;
//		j++;
//	}
//	
//	// BOTTOM
//	j = 0;
//
//	for (int i = nodesNumber / 2; i <  3 * nodesNumber / 4; ++i)
//	{
//		body_.at(i).type_ = IBNodeType::MOVING;
//
//		body_.at(i).cur_pos_.x_ = x_start + width - j * x_step;
//		body_.at(i).ref_pos_.x_ = x_start + width - j * x_step;
//
//		body_.at(i).cur_pos_.y_ = y_start - height;
//		body_.at(i).ref_pos_.y_ = y_start - height;
//		j++;
//	}
//
//	// LEFT
//	j = 0;
//
//	for (int i = 3 * nodesNumber / 4; i < nodesNumber; ++i)
//	{
//		body_.at(i).type_ = IBNodeType::MOVING;
//
//		body_.at(i).cur_pos_.x_ = x_start;
//		body_.at(i).ref_pos_.x_ = x_start;
//
//		body_.at(i).cur_pos_.y_ = y_start - height + j * y_step;
//		body_.at(i).ref_pos_.y_ = y_start - height + j * y_step;
//		j++;
//	}
//}
//
//ImmersedBottomRect::ImmersedBottomRect(int domainX, int domainY, int nodesNumber, Point center, double width, double height) : ImmersedBody(domainX, domainY, nodesNumber, center, (width + height) / M_PI)
//{
//	double x_start = center.x_ - width / 2;
//	double y_start = center.y_ + height / 2;
//
//	double x_step = width / nodesNumber * 4.0;
//	double y_step = height / nodesNumber * 4.0;
//
//	// TOP
//	int j = 0;
//
//	for (int i = 0; i < nodesNumber / 4; ++i)
//	{
//		body_.at(i).type_ = IBNodeType::MOVING;
//
//		body_.at(i).cur_pos_.x_ = x_start + j * x_step;
//		body_.at(i).ref_pos_.x_ = x_start + j * x_step;
//
//		body_.at(i).cur_pos_.y_ = y_start;
//		body_.at(i).ref_pos_.y_ = y_start;
//		j++;
//	}
//
//	// RIGHT
//	j = 0;
//
//	for (int i = nodesNumber / 4; i < nodesNumber / 2; ++i)
//	{
//		body_.at(i).type_ = IBNodeType::MOVING;
//
//		body_.at(i).cur_pos_.x_ = x_start + width;
//		body_.at(i).ref_pos_.x_ = x_start + width;
//
//		body_.at(i).cur_pos_.y_ = y_start - j * y_step;
//		body_.at(i).ref_pos_.y_ = y_start - j * y_step;
//		j++;
//	}
//
//	// BOTTOM
//	j = 0;
//
//	for (int i = nodesNumber / 2; i < 3 * nodesNumber / 4; ++i)
//	{
//		body_.at(i).type_ = IBNodeType::STATIC;
//
//		body_.at(i).cur_pos_.x_ = x_start + width - j * x_step;
//		body_.at(i).ref_pos_.x_ = x_start + width - j * x_step;
//
//		body_.at(i).cur_pos_.y_ = y_start - height;
//		body_.at(i).ref_pos_.y_ = y_start - height;
//		j++;
//	}
//
//	// LEFT
//	j = 0;
//
//	for (int i = 3 * nodesNumber / 4; i < nodesNumber; ++i)
//	{
//		body_.at(i).type_ = IBNodeType::MOVING;
//
//		body_.at(i).cur_pos_.x_ = x_start;
//		body_.at(i).ref_pos_.x_ = x_start;
//
//		body_.at(i).cur_pos_.y_ = y_start - height + j * y_step;
//		body_.at(i).ref_pos_.y_ = y_start - height + j * y_step;
//		j++;
//	}
//}
