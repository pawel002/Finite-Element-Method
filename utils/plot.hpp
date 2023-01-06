#ifndef PLOT_HPP
#define PLOT_HPP

#include "bits/stdc++.h"
#include "mgl2/mgl.h"

void plot(double from, double to, int noPoints, std::vector<double> X, std::vector<double> Y, const char* output_file){
    mglGraph graph;
    graph.SetSize(1920, 1080);

    mglData graph_data_x(noPoints);
    mglData graph_data_y(noPoints);

    double minY = Y[0], maxY = Y[0];

    for(int i=0; i<X.size(); i++){
        graph_data_x.a[i] = X[i];
        graph_data_y.a[i] = Y[i];

        minY = std::min(minY, Y[i]);
        maxY = std::max(maxY, Y[i]);
    }
    
    std::string title = "FEM with " + std::to_string(noPoints) +  " Points";
    graph.Title(title.c_str());
    double margin = 0.1;
    graph.FillBackground(mglColor(255, 255, 255, 0.0));
    graph.LoadFont("adventor font");
    graph.SetFontSize(2.3);
    graph.Light(true);
    graph.SetRanges(from-margin, to+margin, -0.1, maxY+margin);
    graph.SetOrigin(0, 0);
    graph.Plot(graph_data_x, graph_data_y);
    graph.Axis();

    graph.WriteFrame(output_file);
}

#endif