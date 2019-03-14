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
  
  // Compute all numerical integrals
  // For time efficiency, compute once for all
  // The code might not have good readability, but you will only iterate faces once
	

	/*
	Inertia tensor is the integral of the following
	[(x,y,z)]_x^T[(x,y,z)]_x =
		y^2+z^2		-xy		-xz
		-xy		x^2+z^2     -yz
		-xz         -yz x^2+y^2

	*/
	volume_ = 0.0;
	inertiaTensor_.setZero();
	originalCenterOfMass_.setZero();
	

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::Vector3d v1, v2, v3;
		v1 = V.row(F(i, 0));
		v2 = V.row(F(i, 1));
		v3 = V.row(F(i, 2));
		Eigen::Vector3d n;
		// norm(n)= 2area(v1,v2,v3)
		n = (v2 - v1).cross(v3 - v1);

		volume_ = volume_ + ((v1 + v2 + v3).dot(n)) / 18.0;
		for (int j = 0; j < 3; j++)
		{
			originalCenterOfMass_[j] += (n[j] * (v1[j] * v1[j] + v2[j] * v2[j] + v3[j] * v3[j] + v1[j] * v2[j] + v1[j] * v3[j] + v2[j] * v3[j]) / 24.0);
			
		}
	}

	if (volume_ < 1e-14)
		std::cout << "Warning: Volume goes wrong!\n";
	originalCenterOfMass_ /= volume_;
	// Shift the geometry
	V = V.rowwise() - originalCenterOfMass_.transpose();

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::Vector3d v1, v2, v3;
		v1 = V.row(F(i, 0));
		v2 = V.row(F(i, 1));
		v3 = V.row(F(i, 2));
		Eigen::Vector3d n;
		// norm(n)= 2area(v1,v2,v3)
		n = (v2 - v1).cross(v3 - v1);
		double squareTerms[3];

		// compute Squared Term
		for (int j = 0; j < 3; j++)
		{n = (v2 - v1).cross(v3 - v1);
			squareTerms[j] = n[j] * (v1[j] * v1[j] * v1[j] + v2[j] * v2[j] * v2[j] + v3[j] * v3[j] * v3[j] +
				v1[j] * v1[j] * v2[j] + v1[j] * v1[j] * v3[j] + v2[j] * v2[j] * v1[j] + v2[j] * v2[j] * v3[j]
				+ v3[j] * v3[j] * v1[j] + v3[j] * v3[j] * v2[j] + v1[j] * v2[j] * v3[j]) / 60.0;
		}
		// Compute cross term : 0 yz 1 xz 2 xy

		Matrix3d auxV;
		auxV.row(0) = V.row(F(i, 0));
		auxV.row(1) = V.row(F(i, 1));
		auxV.row(2) = V.row(F(i, 2));
		double crossTerms[3];
		double homogeniousTerm = 0.0;
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				for (int r = 0; r < 3; r++)
				{
					int count = ((p == q) + (q == r) + (r == p));
					double factor = 1.0;
					if (count == 3)
						factor = 6.0;
					if (count == 1)
						factor = 2.0;
					homogeniousTerm += (factor * auxV(p, 0) * auxV(q, 1) * auxV(r, 2));
				}
			}
		}
		homogeniousTerm /= 120.0;
		for(int l = 0; l < 3; l++)
		{
			crossTerms[l] = n[l] * homogeniousTerm;
		}

		// Assemble Interia Tensor
		inertiaTensor_(0, 0) += (squareTerms[1] + squareTerms[2]);
		inertiaTensor_(1, 1) += (squareTerms[0] + squareTerms[2]);
		inertiaTensor_(2, 2) += (squareTerms[0] + squareTerms[1]);

		inertiaTensor_(0, 1) -= crossTerms[2];
		inertiaTensor_(1, 0) -= crossTerms[2];

		inertiaTensor_(0, 2) -= crossTerms[1];
		inertiaTensor_(2, 0) -= crossTerms[1];

		inertiaTensor_(1, 2) -= crossTerms[0];
		inertiaTensor_(2, 1) -= crossTerms[0];
			
	} 


	std::cout << "Volume: " << volume_ << std::endl;
	std::cout << "Center of Mass: " << originalCenterOfMass_.transpose() << std::endl;
	std::cout << "Inertia Tensor: " << std::endl << inertiaTensor_ << std::endl;

	
}		
