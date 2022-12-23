#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "bits/stdc++.h"
#include "integrate.hpp"
#include "../Eigen/Sparse"

typedef Eigen::SparseMatrix<double> SparseMatrix;

template <int N>
class Solver
{

    public:
        std::vector<std::function<double(double)>> functions;
        std::vector<std::function<double(double)>> derivatives;
        std::vector<std::pair<double, double>> domains;
        std::vector<double> weights;
        Integrator::Quadrature<N> gl;

        // given eps
        std::function<double(double)> epsilon;

        int n;
        double interval;
        double from;
        double to;

        Solver(std::function<double(double)> epsilon_, double from_, double to_, int n_, Integrator::Quadrature<N> gl_)
        {   
            gl = gl_;
            n = n_;
            epsilon = epsilon_;
            from = from_;
            to = to_;
            // dzivide (from, to) to n-1 intervals
            double interval = (to - from) / (double)(n - 1);

            std::cout << "Tworzenie funkcji testujacych, ich pochodnych i przedzialow.\n";

            // weight vector
            for(int i=0; i<n; i++){
                weights.push_back(1.0);
            }

            // vectors holding fuctions, derivatives, intervals
            for(int i=0; i<n; i++)
            {   
                // def functions
                functions.emplace_back([interval, i](double x){
                    if((i-1)*interval <= x && x < i*interval)
                        return (x - (i-1)*interval) / interval;

                    if(i*interval <= x && x <= (i + 1)*interval)
                        return ((i+1)*interval - x) / interval;

                    return 0.0;
                });

                // def derivatives
                derivatives.emplace_back([interval, i](double x){
                    if((i-1)*interval <= x && x < i*interval)
                        return 1.0 / interval;

                    if(i*interval < x && x <= (i + 1)*interval)
                        return -1.0 / interval;

                    return 0.0;
                });

                // def interval where f(x) != 0
                domains.emplace_back(std::make_pair((i-1)*interval, (i + 1)*interval));
            }
        }

        // helpers
        std::function<double(double)> functionAt(int i){
            return functions.at(i);
        }

        std::function<double(double)> derivativeAt(int i){
            return derivatives.at(i);
        }

        std::pair<double, double> domainAt(int i){
            return domains.at(i);
        }

        double weightAt(int i){
            return weights.at(i);
        }

        double B(int idx1, int idx2)
        {

            std::function<double(double)> a = this->derivativeAt(idx1);
            std::function<double(double)> b = this->derivativeAt(idx2);
            
            // helper for f(x) = v'(x)w'(x)
            std::function<double(double)> f = [a, b](double x){
                return a(x) * b(x);
            };

            std::pair<double, double> integrationDomain = std::make_pair(std::max({domainAt(idx1).first, domainAt(idx2).first, from}),
                                                                         std::min({domainAt(idx1).second, domainAt(idx2).second, to}));
            
            if(integrationDomain.first >= integrationDomain.second)
                return 0.0;

            return functionAt(idx1)(0.0)*functionAt(idx2)(0.0) - gl.integrate(integrationDomain.first, integrationDomain.second, f);
        }

        double L(int idx1)
        {
            std::function<double(double)> a = this->functionAt(idx1);
            std::function<double(double)> b = this->epsilon;
            
            // helper for f(x) =  v(x)/eps(x)
            std::function<double(double)> f = [a, b](double x){
                return a(x) / b(x);
            };

            std::pair<double, double> integrationDomain = std::make_pair(std::max({domainAt(idx1).first, from}),
                                                                         std::min({domainAt(idx1).second, to}));
            
            if(integrationDomain.first >= integrationDomain.second)
                return 0.0;

            return 5*functionAt(idx1)(0.0) - gl.integrate(integrationDomain.first, integrationDomain.second, f);
        }

        void solve()
        {
            // arr holding values before solving equation
            int m = n - 1;
            double Larr[m];
            double BdiagArrFirst[m];
            double BdiagArrSecond[m-1];

            // filling arrs
            std::cout << "Tworzenie macierzy B i wektora L.\n";
            for(int i=0; i<m; i++)
            {
                // L(v) - 2B(w,v)
                Larr[i] = L(i) - 2 * B(m, i);

                // diagonal B(v,v)
                BdiagArrFirst[i] = B(i, i);

                if(i >= 1)
                {
                    // above and below diagonal
                    BdiagArrSecond[i-1] = B(i-1, i);
                }
            }

            SparseMatrix Bmat(m, m);
            Eigen::VectorXd Lvec(m), Wvec;
            Bmat.reserve(Eigen::VectorXi::Constant(m, 3));

            // insert into B and L
            for(int i=0; i<m; i++)
            {
                Lvec[i] = Larr[i];
                Bmat.insert(i, i) = BdiagArrFirst[i];
                if(i >=1 )
                {
                    Bmat.insert(i-1, i) = BdiagArrSecond[i-1];
                    Bmat.insert(i, i-1) = BdiagArrSecond[i-1];
                }
            }

            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > matrixSolver;

            std::cout << "Dekompozycja B.\n";
            matrixSolver.compute(Bmat);

            if(matrixSolver.info()!=Eigen::Success) 
            {
                std::cout << "Decompozycja nieudana.\n";
                return;
            }

            std::cout << "Rozwiazywanie Bx = L.\n";
            Wvec = matrixSolver.solve(Lvec);
            if(matrixSolver.info()!=Eigen::Success) 
            {
                std::cout << "Rozwiazanie nieudane.\n";
                return;
            }

            // fill weight vector
            for(int i=0; i<m; i++)
                weights[i] = Wvec[i];
            
            weights[n - 1] = 2;
            std::cout << "Rozwiazanie udane.\n";
        }

        // evaluates testing sum of testing functions and weights at x
        double evaluateSolution(double x)
        {   
            double sum = 0;
            for(int i=0; i<n; i++)
            {
                if(domainAt(i).first <= x && x <= domainAt(i).second)
                {
                    sum += functionAt(i)(x) * weightAt(i);
                }
            }

            return sum;
        }
};

#endif