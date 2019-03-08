#include "BirdsHook.h"
#include "RigidBodyTemplate.h"
#include "RigidBodyInstance.h"
#include "VectorMath.h"
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

    // TODO: rigid body dynamics

    return false;
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
            return;
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
	// I simply assume there will be configurations of current position (c in notation), which is a 3n-dim vector
	// If not, the program will compute gravity force based on c in each of element in queue variable bodies_
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
			F.segment<3>(3 * i) -= constPart * (c1 - c2);
			F.segment<3>(3 * j) -= constPart * (c2 - c1);
		}
	}
}

void BirdsHook::processGravityFieldForce(VectorXd &F)
{
	// Delete this if you have configuration binding
	int nbodies = (int)bodies_.size();
	for (int i = 0; i < nbodies; i++)
	{
		for (int j = i + 1; j < nbodies; j++)
		{
			Eigen::Vector3d c1 = bodies_[i]->c;
			Eigen::Vector3d c2 = bodies_[j]->c;
			double dist = (c1 - c2).norm();
			double constPart = params_.gravityG * bodies_[i]->mass * bodies_[j]->mass / pow(dist, 3);

			// F = -dV
			F.segment<3>(3 * i) -= constPart * (c1 - c2);
			F.segment<3>(3 * j) -= constPart * (c2 - c1);
		}
	}
}

