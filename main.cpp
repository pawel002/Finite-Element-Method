#include "utils/integrate.hpp"
#include "Eigen/Sparse"
#include "utils/solver.hpp"

#include <chrono>
#include <fstream>
#include "bits/stdc++.h"


// COMPILED USING g++ -I .\Eigen .\main.cpp -o main

int main()
{   
    int n;
    std::cout << "Podaj n - ilosc funkcji testujacych.\n";
    std::cin >> n; 
    double from = 0.0;
    double to = 3.0;
    Integrator::Quadrature<2> gl2;

    std::function<double(double)> eps = [](double x){
        if(0.0 <= x && x <= 1.0)
            return 10.0;
        
        if(1.0 < x && x <= 2.0)
            return 5.0;
        
        if(2.0 < x && x <= 3.0)
            return 1.0;

        return 0.0;
    };

    auto start = std::chrono::high_resolution_clock::now();

    // main solving function
    Solver solver = Solver(eps, from, to, n, gl2);
    solver.solve();

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Obliczanie zajelo " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds.\n";


    // Create and open a text file
    std::ofstream save(".\\saveplot\\plot.txt");
    double x = (to - from) / (double)(n - 1);
    for(int i=0; i<n; i++){
        save << from + x * i << " " << solver.evaluateSolution(from + x * i) << "\n";
    }
    save.close();

    system("python -u \".\\saveplot\\plot.py\"");
}