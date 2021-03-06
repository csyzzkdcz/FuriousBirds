#include "BirdsHook.h"
#include "RigidBodyTemplate.h"
#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include "NewtonSolver.h"
#include "igl/opengl/glfw/imgui/ImGuiHelpers.h"

using namespace Eigen;

void BirdsHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputText("Filename", sceneFile_);
    }
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputFloat("Timestep", &params_.timeStep, 0, 0, 3);
        ImGui::DragFloat("Newton Tolerance", &params_.NewtonTolerance, 0.01, 1e-16, 1e-1, "%.3e", 10);
        ImGui::InputInt("Newton Max Iters", &params_.NewtonMaxIters);
    }
    if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Gravity Enabled", &params_.gravityEnabled);
        ImGui::InputFloat("Gravity G", &params_.gravityG, 0, 0, 3);
    }
}

void BirdsHook::updateRenderGeometry()
{
    int totverts = 0;
    int totfaces = 0;
    for (RigidBodyInstance *rbi : bodies_)
    {
        totverts += rbi->getTemplate().getVerts().rows();
        totfaces += rbi->getTemplate().getFaces().rows();
    }
    renderQ.resize(totverts, 3);
    renderF.resize(totfaces, 3);
    int voffset = 0;
    int foffset = 0;
    for (RigidBodyInstance *rbi : bodies_)
    {
        int nverts = rbi->getTemplate().getVerts().rows();
        for (int i = 0; i < nverts; i++)
            renderQ.row(voffset + i) = (rbi->c + VectorMath::rotationMatrix(rbi->theta)*rbi->getTemplate().getVerts().row(i).transpose()).transpose();
        int nfaces = rbi->getTemplate().getFaces().rows();
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                renderF(foffset + i, j) = rbi->getTemplate().getFaces()(i, j) + voffset;
            }
        }
        voffset += nverts;
        foffset += nfaces;
    }
}


void BirdsHook::initSimulation()
{
    time_ = 0;    
    loadScene();
    updateRenderGeometry();
}

void BirdsHook::tick()
{    
}

bool BirdsHook::simulateOneStep()
{   
    time_ += params_.timeStep;
    Eigen::VectorXd c, cvel, theta, w;
    buildConfiguration(c, cvel, theta, w);
    timeIntegration(c, cvel, theta, w);
    unbuildConfiguration(c, cvel, theta, w);
	/* Check every states
	std::cout << "\nCurrent States:\n\n";
	for (int i = 0; i < (int)bodies_.size(); i++)
	{
		RigidBodyInstance &body = *bodies_[i];
		std::cout << "Body " << i+1 << std::endl;
		std::cout << "position:" << body.c.transpose() << std::endl;
		std::cout << "axisangle:" << body.theta.transpose() << std::endl;
		std::cout << "velocity:" << body.cvel.transpose() << std::endl;
		std::cout << "angular velocity:" << body.w.transpose() << std::endl;

	}
	*/
    return false;
}

void BirdsHook::timeIntegration(Eigen::VectorXd &c, Eigen::VectorXd &cvel, Eigen::VectorXd &theta, Eigen::VectorXd &w)
{
    c += params_.timeStep * cvel;
    int nbodies =  bodies_.size();
    VectorMath vecOp;
	VectorXd previousTheta = theta;
    for(int i = 0; i<nbodies; i++)
    {
        Eigen::Matrix3d R = vecOp.rotationMatrix(theta.segment(3*i, 3)) * vecOp.rotationMatrix(params_.timeStep * w.segment(3*i, 3));
		theta.segment(3 * i, 3) = vecOp.axisAngle(R);
    }
    Eigen::VectorXd force = Eigen::VectorXd::Zero(3*nbodies);
    processGravityFieldForce(c, force);
    for(int i = 0; i<nbodies; i++)
    {
        cvel.segment(3*i, 3) += params_.timeStep * bodies_[i]->massInv * force.segment(3*i, 3);
    }
	Eigen::VectorXd prevw = w;
    bool isSuccess = newtonSolver(w, [this, w, previousTheta](Eigen::VectorXd curw, Eigen::VectorXd &F, Eigen::SparseMatrix<double> *gradF)
                                  {
                                      this->computeValueAndGrad(curw, w, previousTheta, &F, gradF);
                                  }, params_.NewtonMaxIters, params_.NewtonTolerance);
	if(isSuccess)
    {
        std::cout<<"Newton Solver succeeded!!"<<std::endl;
    }
}

void BirdsHook::buildConfiguration(Eigen::VectorXd &c, Eigen::VectorXd &cvel, Eigen::VectorXd &theta, Eigen::VectorXd &w)
{
    int nbodies =  bodies_.size();
    c.resize(3 * nbodies);
    cvel.resize(3 * nbodies);
    theta.resize(3 * nbodies);
    w.resize(3 * nbodies);
    
    for(int i=0;i<nbodies;i++)
    {
        c.segment(3*i, 3) = bodies_[i]->c;
        cvel.segment(3*i, 3) = bodies_[i]->cvel;
        theta.segment(3*i, 3) = bodies_[i]->theta;
        w.segment(3*i, 3) = bodies_[i]->w;
    }
}

void BirdsHook::unbuildConfiguration(Eigen::VectorXd c, Eigen::VectorXd cvel, Eigen::VectorXd theta, Eigen::VectorXd w)
{
    int nbodies =  bodies_.size();
    
    for(int i=0;i<nbodies;i++)
    {
        bodies_[i]->c = c.segment(3*i, 3);
        bodies_[i]->cvel = cvel.segment(3*i, 3);
        bodies_[i]->theta = theta.segment(3*i, 3);
        bodies_[i]->w = w.segment(3*i, 3);
    }
    
}

void BirdsHook::loadScene()
{
    for (RigidBodyInstance *rbi : bodies_)
        delete rbi;
    for (RigidBodyTemplate *rbt : templates_)
        delete rbt;
    bodies_.clear();
    templates_.clear();

    std::string prefix;
    std::string scenefname = std::string("scenes/") + sceneFile_;
    std::ifstream ifs(scenefname);
    if (!ifs)
    {
        // run from the build directory?
        prefix = "../";
        scenefname = prefix + scenefname;        
        ifs.open(scenefname);
        if(!ifs)
        {
            prefix = "../../";   // For xcode && VS build/release
            scenefname = prefix + std::string("scenes/") + sceneFile_;
            ifs.open(scenefname);
            if(!ifs)
                return;
        }
            
    }
        

    int nbodies;
    ifs >> nbodies;
    for (int body = 0; body < nbodies; body++)
    {
        std::string meshname;
        ifs >> meshname;
        meshname = prefix + std::string("meshes/") + meshname;
        double scale;
        ifs >> scale;
        RigidBodyTemplate *rbt = new RigidBodyTemplate(meshname, scale);
        double rho;
        ifs >> rho;
        Eigen::Vector3d c, theta, cvel, w;
        for (int i = 0; i < 3; i++)
            ifs >> c[i];
        for (int i = 0; i < 3; i++)
            ifs >> theta[i];
        for (int i = 0; i < 3; i++)
            ifs >> cvel[i];
        for (int i = 0; i < 3; i++)
            ifs >> w[i];
        RigidBodyInstance *rbi = new RigidBodyInstance(*rbt, c, theta, cvel, w, rho);
        templates_.push_back(rbt);
        bodies_.push_back(rbi);
    }
}



/////////////////////////////////////////////////////////////////////////////////////
/////
///// Gravity Force Computation 
/////
/////////////////////////////////////////////////////////////////////////////////////

void BirdsHook::processGravityFieldForce(const VectorXd &c, VectorXd &F)
{
	// I did not set F to be zero in this code.
	int nbodies = (int)bodies_.size();
	for (int i = 0; i < nbodies; i++)
	{
		for (int j = i + 1; j < nbodies; j++)
		{
			Eigen::Vector3d c1 = c.segment<3>(3 * i);
			Eigen::Vector3d c2 = c.segment<3>(3 * j);
			double dist = (c1 - c2).norm();
			double constPart = params_.gravityG * bodies_[i]->mass * bodies_[j]->mass/pow(dist,3);

			// F = -dV
			F.segment<3>(3 * i) -= (constPart * (c1 - c2));
			F.segment<3>(3 * j) -= (constPart * (c2 - c1));
		}
	}
}

void BirdsHook::processGravityFieldHessian(const Eigen::VectorXd & c, std::vector<Eigen::Triplet<double>>& HTriplet)
{
	int nbodies = (int)bodies_.size();
	for (int i = 0; i < nbodies; i++)
	{
		for (int j = i + 1; j < nbodies; j++)
		{
			Eigen::Vector3d c1 = c.segment<3>(3 * i);
			Eigen::Vector3d c2 = c.segment<3>(3 * j);
			double dist = (c1 - c2).norm();
			double constPart = params_.gravityG * bodies_[i]->mass * bodies_[j]->mass / pow(dist, 3);

			// H = -dF
			Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
			Eigen::Matrix3d localH = I * constPart - 3 * constPart * (c1 - c2)*(c1 - c2).transpose() / dist / dist;
			
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					HTriplet.push_back(Eigen::Triplet<double>(3 * i + k, 3 * i + l, localH(k, l)));
					HTriplet.push_back(Eigen::Triplet<double>(3 * j + k, 3 * j + l, localH(k, l)));

					HTriplet.push_back(Eigen::Triplet<double>(3 * i + k, 3 * j + l, -localH(k, l)));
					HTriplet.push_back(Eigen::Triplet<double>(3 * j + k, 3 * i + l, -localH(k, l)));
				}
			}
		}
	}
}



void BirdsHook::computeValueAndGrad(Eigen::VectorXd curw, Eigen::VectorXd prevw, 
	Eigen::VectorXd prevtheta, Eigen::VectorXd *f, Eigen::SparseMatrix<double> *df)
{
    // Right now, in our case, d_{\theta} V(q^{i+1}) = 0
    int nbodies =  bodies_.size();
    VectorMath vecOp;
    std::vector<Eigen::Triplet<double> > triplet;
    if(f != NULL)
        f->resize(3*nbodies);
    for(int i=0;i<nbodies;i++)
    {
		Eigen::Matrix3d MI = bodies_[i]->intertiaTensor;
		Eigen::Vector3d currentw = curw.segment(3 * i, 3);
		Eigen::Vector3d previousw = prevw.segment(3 * i, 3);
		Eigen::Vector3d previoustheta = prevtheta.segment(3 * i, 3);
		Eigen::Vector3d fval = (currentw.transpose()* MI * vecOp.TMatrix(-params_.timeStep * currentw).inverse() * vecOp.TMatrix(previoustheta) -
			previousw.transpose() * MI * vecOp.TMatrix(params_.timeStep * previousw).inverse() * vecOp.TMatrix(previoustheta)).transpose();
		if (f != NULL)
			f->segment(3 * i, 3) = fval;
		// Just give up for derivative of T matrix term since it is a small term.
		Eigen::Matrix3d localDeriv = (MI * vecOp.TMatrix(-params_.timeStep * currentw).inverse() * vecOp.TMatrix(previoustheta)).transpose();
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
            {
                triplet.push_back(Eigen::Triplet<double>(3*i + j, 3*i + k, localDeriv(j,k)));
            }
    }
    if(df != NULL)
    {
        df->resize(3*nbodies, 3*nbodies);
        df->setFromTriplets(triplet.begin(), triplet.end());
    }
    
}

/////////////////////////////////////////////////////////////////////////////////////
/////
///// Test If Computation Makes Sense
/////
/////////////////////////////////////////////////////////////////////////////////////

void BirdsHook::testValueAndGrad()
{
    Eigen::VectorXd c, cvel, theta, w;
    buildConfiguration(c, cvel, theta, w);
    
    Eigen::VectorXd prevw = w;
	Eigen::VectorXd prevtheta = theta;
    timeIntegration(c, cvel, theta, w);
    Eigen::VectorXd f;
    Eigen::SparseMatrix<double> df;
    
    computeValueAndGrad(w, prevw, prevtheta, &f, &df);
    
    Eigen::VectorXd dir = Eigen::VectorXd::Random(w.size());
    dir.normalized();
    
    for(int i = 4; i< 14; i++)
    {
        double eps = powf(10, -i);
        Eigen::VectorXd updatedw = w + eps * dir;
        Eigen::VectorXd updatedf;
        computeValueAndGrad(updatedw, prevw, prevtheta, &updatedf, NULL);
        
        std::cout<<"EPS is: "<<eps<<std::endl;
        std::cout<<"The norm of directional derivative is: "<<(df*dir).norm()<<" The norm of finite difference is: "<< ((updatedf - f)/eps).norm() <<std::endl;
        std::cout<<"The error is: "<<( df*dir - (updatedf - f)/eps).norm()<<std::endl;
    }
}

void BirdsHook::testGravityFieldForce()
{
	Eigen::VectorXd c, cvel, theta, w;
	buildConfiguration(c, cvel, theta, w);

	Eigen::VectorXd f;
	Eigen::SparseMatrix<double> H;
	std::vector<Eigen::Triplet<double>> HCoef;

	f.resize(c.size());
	f.setZero();
	processGravityFieldForce(c, f);
	processGravityFieldHessian(c, HCoef);

	H.resize(c.size(), c.size());
	H.setFromTriplets(HCoef.begin(), HCoef.end());
	
	Eigen::VectorXd dir = Eigen::VectorXd::Random(c.size());
	dir.normalized();

	for (int i = 4; i< 14; i++)
	{
		double eps = powf(10, -i);
		Eigen::VectorXd updatedc = c + eps * dir;
		Eigen::VectorXd updatedf;
		updatedf.resize(c.size());
		updatedf.setZero();
		processGravityFieldForce(updatedc, updatedf);

		std::cout << "EPS is: " << eps << std::endl;
		std::cout << "The norm of Hessian is: " << (H*dir).norm() << " The norm of finite difference is: " << ((updatedf - f) / eps).norm() << std::endl;
		std::cout << "The error is: " << (H*dir + (updatedf - f) / eps).norm() << std::endl;
	}
}
