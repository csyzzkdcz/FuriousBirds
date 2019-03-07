#include "RigidBodyTemplate.h"
#include <iostream>
#include <igl/readOBJ.h>
#include <Eigen/Dense>
#include <fstream>
#include <map>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

RigidBodyTemplate::RigidBodyTemplate(const std::string &meshFilename, double scale) : volume_(0), radius_(0)
{
    inertiaTensor_.setZero();

    igl::readOBJ(meshFilename, V, F);

    V *= scale;

    initialize();
}

RigidBodyTemplate::~RigidBodyTemplate()
{    
}

void RigidBodyTemplate::initialize()
{
  // TODO compute quantities such as volume, center of mass, and intertia tensor
  // translate body so its center of mass is the world origin

  // Compute all numerical integrals
  // For time efficiency, compute once for all
  // The code might not have good readability, but you will only iterate faces once
	volume_ = 0.0;
	inertiaTensor_.setZero();
	originalCenterOfMass_.setZero();
	
	/*
	Inertia tensor is the integral of the following
	[(x,y,z)]_x^T[(x,y,z)]_x = 
	    y^2+z^2		-xy		-xz
		-xy		x^2+z^2     -yz
		-xz         -yz x^2+y^2

	*/

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::Vector3d v1, v2, v3;
		v1 = V.row(F(i, 0));
		v2 = V.row(F(i, 1));
		v3 = V.row(F(i, 2));
		Eigen::Vector3d n;
		n = (v2 - v1).cross(v3 - v1);
	
		volume_ = volume_ + ((v1 + v2 + v3).dot(n)) / 18.0;  // norm(n)= 2area(v1,v2,v3)

		for (int j = 0; j < 3; j++)
		{
			originalCenterOfMass_[j] += (n[j] * (v1[j] * v1[j] + v2[j] * v2[j] + v3[j] * v3[j] + v1[j] * v2[j] + v1[j] * v3[j] + v2[j] * v3[j]) / 24.0);
		}

	} 

	originalCenterOfMass_ /= volume_;

	std::cout << "Volume: " << volume_ << std::endl;
	std::cout << "Center of Mass: " << originalCenterOfMass_.transpose() << std::endl;

	// Shift the geometry
	V = V.rowwise() - originalCenterOfMass_.transpose();
}		
