#include "PhysicsHook.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/ViewerData.h>
#include <deque>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include "SimParameters.h"

class RigidBodyTemplate;
class RigidBodyInstance;

class BirdsHook : public PhysicsHook
{
public:
    BirdsHook() : PhysicsHook(), sceneFile_("box.scn") {}

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);

    virtual void initSimulation();

    virtual void mouseClicked(double x, double y, int button)
    {
    }

    virtual void updateRenderGeometry();

    virtual void tick();

    virtual bool simulateOneStep();

    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
    {
        viewer.data().clear();
        viewer.data().set_mesh(renderQ, renderF);
    }
    
    void testValueAndGrad();
    
private:
    void loadScene();

	void processGravityFieldForce(const Eigen::VectorXd &c, Eigen::VectorXd &F);
	void processGravityFieldForce(Eigen::VectorXd &F);
    
    void buildConfiguration(Eigen::VectorXd &c, Eigen::VectorXd &cvel, Eigen::VectorXd &theta, Eigen::VectorXd &w);
    void unbuildConfiguration(Eigen::VectorXd c, Eigen::VectorXd cvel, Eigen::VectorXd theta, Eigen::VectorXd w);
    
    void timeIntegration(Eigen::VectorXd &c, Eigen::VectorXd &cvel, Eigen::VectorXd &theta, Eigen::VectorXd &w);

    void computeValueAndGrad(Eigen::VectorXd curw, Eigen::VectorXd prevw, Eigen::VectorXd *f, Eigen::SparseMatrix<double> *df);  // This function is used to update w
    
    
    double time_;
    SimParameters params_;
    std::string sceneFile_;

    std::vector<RigidBodyTemplate *> templates_;
    std::vector<RigidBodyInstance *> bodies_;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;

};
