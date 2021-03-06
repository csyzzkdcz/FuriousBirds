#ifndef NEWTONSOLVER_H
#define NEWTONSOLVER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

bool newtonSolver(Eigen::VectorXd &x, std::function<void (Eigen::VectorXd, Eigen::VectorXd &, Eigen::SparseMatrix<double> *)> _computeFAndGradF, int NewtonMaxIters, double NewtonTolerance)
{
    for(int i=0; i<NewtonMaxIters; i++)
    {
		//std::cout << i << std::endl;
        Eigen::VectorXd F,dx;
        Eigen::SparseMatrix<double> gradF;
        _computeFAndGradF(x, F, &gradF);
		//std::cout << F.norm() << std::endl;
        if(F.norm()<NewtonTolerance)
        {
            std::cout<<"Optimal station reached in "<<i<<" iterations"<<std::endl;
            return true;
        }
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
        solver.compute(gradF);
        dx = solver.solve(-F);
        x = x + dx;
    }
    //std::cout<<"Maximun iteration reached !!"<<std::endl;
    return true;
}

#endif
