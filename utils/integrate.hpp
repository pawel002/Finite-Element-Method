#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include "bits/stdc++.h"

namespace Integrator {

    template <int N>
    class Quadrature {
    public:
        enum {eDEGREE = N};

        template <typename Function>
        double integrate(double a, double b, Function f) {
            double p = (b - a) / 2;
            double q = (b + a) / 2;
            const LegendrePolynomial& legpoly = s_LegendrePolynomial;

            double sum = 0;
            for (int i = 1; i <= eDEGREE; ++i) {
                sum += legpoly.weight(i) * f(p * legpoly.root(i) + q);
            }
        
            return p * sum;
        }

    private:
        class LegendrePolynomial {
        public:
            LegendrePolynomial () {
                for (int i = 0; i <= eDEGREE; ++i) {
                    double dr = 1;

                    Evaluation eval(cos(M_PI * (i - 0.25) / (eDEGREE + 0.5)));
                    do {
                        dr = eval.v() / eval.d();
                        eval.evaluate(eval.x() - dr);
                    } while (fabs (dr) > 2e-16);

                    this->_r[i] = eval.x();
                    this->_w[i] = 2 / ((1 - eval.x() * eval.x()) * eval.d() * eval.d());
                }
            }

            double root(int i) const 
            { 
                return this->_r[i]; 
            }
            
            double weight(int i) const 
            {
                return this->_w[i]; 
            }
            
        private:
            double _r[eDEGREE + 1];
            double _w[eDEGREE + 1];

            class Evaluation {
            public:
                explicit Evaluation (double x) : _x(x), _v(1), _d(0) {
                    this->evaluate(x);
                }

                void evaluate(double x) {
                    this->_x = x;

                    double vsub1 = x;
                    double vsub2 = 1;
                    double f     = 1 / (x * x - 1);
                
                    for (int i = 2; i <= eDEGREE; ++i) {
                        this->_v = ((2 * i - 1) * x * vsub1 - (i - 1) * vsub2) / i;
                        this->_d = i * f * (x * this->_v - vsub1);

                        vsub2 = vsub1;
                        vsub1 = this->_v;
                    }
                }

                double v() const { return this->_v; }
                double d() const { return this->_d; }
                double x() const { return this->_x; }

            private:
                double _x;
                double _v;
                double _d;
            };
        };

        static LegendrePolynomial s_LegendrePolynomial;
    };

    template <int N>
    typename Quadrature<N>::LegendrePolynomial Quadrature<N>::s_LegendrePolynomial;
}


#endif