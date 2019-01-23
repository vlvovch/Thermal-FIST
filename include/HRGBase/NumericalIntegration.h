/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H

/**
 * \file NumericalIntegration.h
 * 
 * Collection of Gauss-Legendre and Gauss-Laguerre quadratures
 * used in numerical integrations.
 */

#include <vector>

namespace thermalfist {

  /// \brief Contains various Gauss-Legendre and Gauss-Laguerre quadratures
  ///        used in numerical integrations.
  namespace NumericalIntegration {

    /**
     * \brief Nodes of the 5-point Gauss-Legendre quadrature.
     * 
     */
    const double coefficients_xleg5[5] = { -0.906179845939, -0.538469310106,
      0., 0.538469310106, 0.906179845939 };
    
    /**
     * \brief Weights of the 5-point Gauss-Legendre quadrature.
     * 
     */
    const double coefficients_wleg5[5] = { 0.236926885056, 0.478628670499, 0.568888888889, 0.478628670499,
      0.236926885056 };

    /**
     * \brief Nodes of the 10-point Gauss-Legendre quadrature.
     * 
     */
    const double coefficients_xleg10[10] = { -0.973906528517, -0.865063366689, -0.679409568299, -0.433395394129,
      -0.148874338982, 0.148874338982, 0.433395394129, 0.679409568299,
      0.865063366689, 0.973906528517 };
    
    /**
     * \brief Weights of the 10-point Gauss-Legendre quadrature.
     * 
     */
    const double coefficients_wleg10[10] = { 0.0666713443087, 0.149451349151, 0.219086362516, 0.269266719310,
      0.295524224715, 0.295524224715, 0.269266719310, 0.219086362516,
      0.149451349151, 0.0666713443087 };

    /**
     * \brief Nodes of the 32-point Gauss-Legendre quadrature.
     * 
     */
    const double coefficients_xleg32[32] = { -0.997263861849, -0.985611511545, -0.964762255588, -0.934906075938,
      -0.896321155766, -0.849367613733, -0.794483795968, -0.732182118740,
      -0.663044266930, -0.587715757241, -0.506899908932, -0.421351276131,
      -0.331868602282, -0.239287362252, -0.144471961583, -0.048307665688,
      0.048307665688, 0.144471961583, 0.239287362252, 0.331868602282,
      0.421351276131, 0.506899908932, 0.587715757241, 0.663044266930,
      0.732182118740, 0.794483795968, 0.849367613733, 0.896321155766,
      0.934906075938, 0.964762255588, 0.985611511545, 0.997263861849 };
    
    /**
     * \brief Weights of the 32-point Gauss-Legendre quadrature.
     * 
     */
    const double coefficients_wleg32[32] = { 0.00701861000947, 0.0162743947309, 0.0253920653093, 0.0342738629130,
      0.0428358980222, 0.0509980592624, 0.0586840934785, 0.0658222227764,
      0.0723457941088, 0.0781938957871, 0.0833119242269, 0.0876520930044,
      0.0911738786958, 0.0938443990808, 0.0956387200793, 0.0965400885147,
      0.0965400885147, 0.0956387200793, 0.0938443990808, 0.0911738786958,
      0.0876520930044, 0.0833119242269, 0.0781938957871, 0.0723457941088,
      0.0658222227764, 0.0586840934785, 0.0509980592624, 0.0428358980222,
      0.0342738629130, 0.0253920653093, 0.0162743947309, 0.00701861000947 };

    /**
     * \brief Nodes of the 40-point Gauss-Legendre quadrature.
     * 
     */
    const double coefficients_xleg40[40] = { -0.998237709711, -0.990726238699, -0.977259949984, -0.957916819214,
      -0.932812808279, -0.902098806969, -0.865959503212, -0.824612230833,
      -0.778305651427, -0.727318255190, -0.671956684614, -0.612553889668,
      -0.549467125095, -0.483075801686, -0.413779204372, -0.341994090826,
      -0.268152185007, -0.192697580701, -0.116084070675, -0.038772417506,
      0.038772417506, 0.116084070675, 0.192697580701, 0.268152185007,
      0.341994090826, 0.413779204372, 0.483075801686, 0.549467125095,
      0.612553889668, 0.671956684614, 0.727318255190, 0.778305651427,
      0.824612230833, 0.865959503212, 0.902098806969, 0.932812808279,
      0.957916819214, 0.977259949984, 0.990726238699, 0.998237709711 };
    
    /**
     * \brief Weights of the 40-point Gauss-Legendre quadrature.
     * 
     */
    const double coefficients_wleg40[40] = { 0.00452127709853, 0.0104982845312, 0.0164210583819, 0.0222458491942,
      0.0279370069800, 0.0334601952825, 0.0387821679745, 0.0438709081857,
      0.0486958076351, 0.0532278469839, 0.0574397690994, 0.0613062424929,
      0.0648040134566, 0.0679120458152, 0.0706116473913, 0.0728865823958,
      0.0747231690580, 0.0761103619006, 0.0770398181642, 0.0775059479784,
      0.0775059479784, 0.0770398181642, 0.0761103619006, 0.0747231690580,
      0.0728865823958, 0.0706116473913, 0.0679120458152, 0.0648040134566,
      0.0613062424929, 0.0574397690994, 0.0532278469839, 0.0486958076351,
      0.0438709081857, 0.0387821679745, 0.0334601952825, 0.0279370069800,
      0.0222458491942, 0.0164210583819, 0.0104982845312, 0.00452127709853 };

    /**
     * \brief Nodes of the 32-point Gauss-Laguerre quadrature.
     * 
     */
    const double coefficients_xlag32[32] = { 0.0444893658333, 0.234526109520, 0.576884629302, 1.07244875382,
      1.72240877644, 2.52833670643, 3.49221327302, 4.61645676975,
      5.90395850417, 7.35812673319, 8.98294092421, 10.7830186325,
      12.7636979867, 14.9311397555, 17.2924543367, 19.8558609403,
      22.6308890132, 25.6286360225, 28.8621018163, 32.3466291540,
      36.1004948058, 40.1457197715, 44.5092079958, 49.2243949873,
      54.3337213334, 59.8925091621, 65.9753772879, 72.6876280907,
      80.1874469779, 88.7353404179, 98.8295428683, 111.751398098 };
    
    /**
     * \brief Weights of the 32-point Gauss-Laguerre quadrature.
     * 
     */
    const double coefficients_wlag32[32] = { 0.114187105768, 0.266065216898, 0.418793137325, 0.572532846500,
      0.727648788381, 0.884536719340, 1.04361887589, 1.20534927415,
      1.37022133852, 1.53877725647, 1.71161935269, 1.88942406345,
      2.07295934025, 2.26310663400, 2.46088907249, 2.66750812640,
      2.88439209292, 3.11326132704, 3.35621769260, 3.61586985648,
      3.89551304495, 4.19939410471, 4.53311497853, 4.90427028761,
      5.32350097202, 5.80633321423, 6.37661467416, 7.07352658071,
      7.96769350930, 9.20504033128, 11.1630130908, 15.3901804153 };

    /**
     * \brief Nodes of the 32-point Gauss-Legendre quadrature in the interval [0,1].
     * 
     */
    const double coefficients_xleg32_zeroone[32] = { 0.00136806907526, 0.00719424422737, 0.0176188722062,
      0.0325469620311, 0.0518394221170, 0.0753161931337, 0.102758102016,
      0.133908940630, 0.168477866535, 0.206142121380, 0.246550045534,
      0.289324361935, 0.334065698859, 0.380356318874, 0.427764019209,
      0.475846167156, 0.524153832844, 0.572235980791, 0.619643681126,
      0.665934301141, 0.710675638065, 0.753449954466, 0.793857878620,
      0.831522133465, 0.866091059370, 0.897241897984, 0.924683806866,
      0.948160577883, 0.967453037969, 0.982381127794, 0.992805755773,
      0.998631930925 };
    
    /**
     * \brief Weights of the 32-point Gauss-Legendre quadrature in the interval [0,1].
     * 
     */
    const double coefficients_wleg32_zeroone[32] = { 0.00350930500474, 0.00813719736545, 0.0126960326546,
      0.0171369314565, 0.0214179490111, 0.0254990296312, 0.0293420467393,
      0.0329111113882, 0.0361728970544, 0.0390969478935, 0.0416559621135,
      0.0438260465022, 0.0455869393479, 0.0469221995404, 0.0478193600396,
      0.0482700442574, 0.0482700442574, 0.0478193600396, 0.0469221995404,
      0.0455869393479, 0.0438260465022, 0.0416559621135, 0.0390969478935,
      0.0361728970544, 0.0329111113882, 0.0293420467393, 0.0254990296312,
      0.0214179490111, 0.0171369314565, 0.0126960326546, 0.00813719736545,
      0.00350930500474 };

    /**
     * Integrates a function f(x,y) 
     * in range 0 < x < \infty, ay <= y <= by
     * using the combined 32-point
     * Gauss-Laguerre and the 32-point Gauss-Legendre quadrature.
     * 
     * \param func A point to a function to be integrate.
     * \param ay Left limit of integration for variable y.
     * \param by Right limit of integration for variable y.
     * \return Result of the integration.
     */
    double Integrate2DLaguerre32Legendre32(double(*func)(double, double), double ay, double by);

    /**
     * Populates the nodes and weights for integrating
     * a function f(x,y)
     * in the range 0 < x < \infty, ay <= y <= by
     * using the combined 32-point
     * Gauss-Laguerre and the 32-point Gauss-Legendre quadrature.
     * 
     * \param [in] ay Left limit of integration for variable y.
     * \param [in] by Right limit of integration for variable y.
     * \param [out] xlag Gauss-Laguerre nodes for the variable x. 
     * \param [out] wlag Gauss-Laguerre weights for the variable x. 
     * \param [out] xleg Gauss-Legendre nodes for the variable y. 
     * \param [out] wleg Gauss-Legendre weights for the variable y. 
     */
    void GetCoefs2DLaguerre32Legendre32(double ay, double by, std::vector<double> *xlag, std::vector<double> *wlag, std::vector<double> *xleg, std::vector<double> *wleg);

    /**
     * Populates the nodes and weights for integrating
     * a function f(x,y)
     * in the range ax < x < bx, ay <= y <= by
     * using two 32-point Gauss-Legendre quadratures.
     * 
     * \param [in] ax Left limit of integration for variable x.
     * \param [in] bx Right limit of integration for variable x.
     * \param [in] ay Left limit of integration for variable y.
     * \param [in] by Right limit of integration for variable y.
     * \param [out] xlag Gauss-Legendre nodes for the variable x. 
     * \param [out] wlag Gauss-Legendre weights for the variable x. 
     * \param [out] xleg Gauss-Legendre nodes for the variable y. 
     * \param [out] wleg Gauss-Legendre weights for the variable y. 
     */
    void GetCoefs2DLegendre32Legendre32(double ax, double bx, double ay, double by, std::vector<double> *xleg1, std::vector<double> *wleg1, std::vector<double> *xleg2, std::vector<double> *wleg2);

    /**
     * Populates the nodes and weights for integrating
     * a function f(x)
     * in the range ax < x < bx
     * using the 32-point Gauss-Legendre quadrature.
     * 
     * \param [in] a Left limit of integration.
     * \param [in] b Right limit of integration
     * \param [out] x Gauss-Legendre nodes.
     * \param [out] w Gauss-Legendre weights.
     */
    void GetCoefsIntegrateLegendre32(double a, double b, std::vector<double> *x, std::vector<double> *w);

    /**
     * Populates the nodes and weights for integrating
     * a function f(x)
     * in the range ax < x < bx
     * using the 10-point Gauss-Legendre quadrature.
     * 
     * \param [in] a Left limit of integration.
     * \param [in] b Right limit of integration
     * \param [out] x Gauss-Legendre nodes.
     * \param [out] w Gauss-Legendre weights.
     */
    void GetCoefsIntegrateLegendre10(double a, double b, std::vector<double> *x, std::vector<double> *w);

    /**
     * Populates the nodes and weights for integrating
     * a function f(x)
     * in the range ax < x < bx
     * using the 5-point Gauss-Legendre quadrature.
     * 
     * \param [in] a Left limit of integration.
     * \param [in] b Right limit of integration
     * \param [out] x Gauss-Legendre nodes.
     * \param [out] w Gauss-Legendre weights.
     */
    void GetCoefsIntegrateLegendre5(double a, double b, std::vector<double> *x, std::vector<double> *w);

    /**
     * Populates the nodes and weights for integrating
     * a function f(x)
     * in the range \param ax < x < \param bx
     * using the 40-point Gauss-Legendre quadrature.
     * 
     * \param [in] a Left limit of integration.
     * \param [in] b Right limit of integration
     * \param [out] x Gauss-Legendre nodes.
     * \param [out] w Gauss-Legendre weights.
     */
    void GetCoefsIntegrateLegendre40(double a, double b, std::vector<double> *x, std::vector<double> *w);

    /**
     * Populates the nodes and weights for integrating
     * a function f(x)
     * in the range 0 < x < \infty
     * using the 32-point Gauss-Laguerre quadrature.
     * 
     * \param [out] x Gauss-Legendre nodes.
     * \param [out] w Gauss-Legendre weights.
     */
    void GetCoefsIntegrateLaguerre32(std::vector<double> *x, std::vector<double> *w);

  }

} // namespace thermalfist

#endif
