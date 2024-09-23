namespace thermalfist {

  namespace xMath {
// Computing Lambert W function (0-branch) using Halley's method
// The desired accuracy is within 10 * epsilon where epsilon is the machine precision
// Initial guess is z = 0
    template<typename T>
    T LambertW0(T z) {
      if (z < -(1 / M_E)) {
//        std::cout << "LambertW0(x): x < -(1 / M_E) = " << -(1 / M_E) << ". Returning W0(1/e)." << endl;
        z = -(1 / M_E);
      }

      T Wcur = 0, Wprev = 0;
      const int MAXITERS = 200;
      const double TOL = 10. * std::numeric_limits<double>::epsilon();
      for (int i = 0; i < MAXITERS; ++i) {
        T Wexp = std::exp(Wprev);
        T denom = (Wexp * (2. + Wprev * (2. + Wprev)) + (2. + Wprev) * z) / (2. * (1. + Wprev));

        Wcur = Wprev - ((Wprev * Wexp - z) / denom);

        if (std::abs(Wprev / Wcur - 1.) < TOL)
          break;

        Wprev = Wcur;

        if (i == MAXITERS - 1) {
//          std::cout << "LambertW0(x): Max iters reached for x = " << z << ".";
//          std::cout << "Relative error estimate: " << T(abs(Wprev/Wcur - 1.)) << endl;
        }

      }

      return Wcur;
    }

  }

}