#!/usr/bin/env python3
#
# [1] Smoothing splines:
#     https://docs.scipy.org/doc/scipy/tutorial/
#     interpolate/smoothing_splines.html
#
# [2] scipy.interpolate.splrep:
#     https://docs.scipy.org/doc/scipy/reference/generated/
#     scipy.interpolate.splrep.html#scipy.interpolate.splrep
#
# [3] Univariate spline in the B-spline basis:
#     https://docs.scipy.org/doc/scipy/reference/generated/
#     scipy.interpolate.BSpline.html#scipy.interpolate.BSpline
#

import numpy as np
from scipy.interpolate import splrep, BSpline
import matplotlib.pyplot as plt


def print_smspl(tck, x, y, w=None):
    print(f'degree of the spline is {tck[2]}')
    print(f'knots are {tck[0]} (N={len(tck[0])})')
    # number of B-spline coefficients are len(knots) - k - 1
    nb = len(tck[0]) - tck[2] - 1
    print(f'number of B-spline coefficients are (Nb={nb})')
    ys = BSpline(*tck)(x)
    ares = (y-ys)**2
    #  print(f'squared residuals point by point are\n{ares}')
    if w is None:
        print(f'sum of squared residuals is {sum(ares):.3f}')
    else:
        print(f'weighted sum of squared residuals is'
              f' {sum(ares*(w**2)):.3f}')


def data_label(data_name):
    t = data_name.split('_')
    if t[0] == 'Pim':
        t[0] = r'$\pi^-$'
    elif t[0] == 'Pip':
        t[0] = r'$\pi^+$'
    elif t[0] == 'Km':
        t[0] = r'$K^-$'
    elif t[0] == 'Kp':
        t[0] = r'$K^+$'
    return ' '.join(t)


def bes3_data(data_name, sys_err, tsc=0, show_pict=True):
    dat = []
    #  file = f'{data_name}.dat'
    file = f'Pict/TrkEff/{data_name}.dat'
    with open(file,'r') as fr:
        for txt in fr:
            t = [float(v) for v in txt[4:].strip(' []\n').split(',')]
            dat.append(t)

    x = np.array(dat[0])
    y = np.array(dat[1])
    err_y = np.array(dat[2])
    err_y = np.sqrt(err_y**2 + sys_err**2)

    # weights
    wy = 1./err_y

    tck = splrep(x, y, w=wy, s=0)

    if isinstance(tsc, float):
        sc = tsc
    else:
        sc = len(x) - np.sqrt(2*len(x))
        if tsc == 1:
            sc = len(x)
        elif tsc == 2:
            sc = len(x) + np.sqrt(2*len(x))
    tck_s = splrep(x, y, w=wy, s=sc)

    # chi^2 and Ndof
    #  print_smspl(tck, x, y, w=wy)
    #  print_smspl(tck_s, x, y, w=wy)
    ch2 = sum(((y-BSpline(*tck_s)(x))/err_y)**2)
    # number of points minus number of B-spline coefficients
    ndof = len(x) - (len(tck_s[0])-tck_s[2]-1)

    # test B-spline as root C++ function
    #  write_Bcpp(tck_s, f'Bsp_{data_name}', *dat)

    # plot figure: only one object in memory and auto cleaning
    fig = plt.figure(num=1, clear=True)
    Pmin = x[0] - (x[1]-x[0])/2
    Pmax = x[-1] + (x[1]-x[0])/2
    xnew = np.arange(Pmin, Pmax, 0.001)
    Rmin = 0.9
    if y[0] < Rmin:
        Rmin = 0.1 * np.floor((y[0]-0.01)*10)
    Rmax = 1.1
    if y[0] > Rmax:
        Rmax = 0.1 * np.ceil((y[0]+0.01)*10)
    plt.ylim([Rmin, Rmax])

    plt.errorbar(x, y, err_y, fmt='ko',
                 label=f'{data_label(data_name)}')
    plt.plot(xnew, BSpline(*tck)(xnew), '--', label='Cubic spline')
    plt.plot(xnew, BSpline(*tck_s)(xnew), '-',
             label=f'Smooth spline, s={sc:.1f}')
    plt.plot([], [], ' ', label=f'χ² / ndof = {ch2:.1f} / {ndof}')
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [3, 0, 1, 2]
    plt.legend([handles[idx] for idx in order],
               [labels[idx] for idx in order],
               #  title='Explanatory text',
               loc='best')
    plt.grid(True)
    plt.xlabel(r'$P_t$, GeV/$c$', loc='right')
    plt.ylabel(r'$\epsilon(\text{data})/\epsilon(\text{MC})$',
               loc='top')
    fig.savefig(f'rat_Bspl_{data_name}.pdf',
                format="pdf", bbox_inches="tight")
    if show_pict:
        plt.show()
    else:
        print(f'{data_name:9s}: s={sc:.1f} '\
              f'χ²/ndof= {ch2:.1f}/{ndof:2d} = {ch2/ndof:.2f}')

    return get_Fun(tck_s, f'Bsp_{data_name}')


# Write C++ function:
def get_Fun(tck, func_name):
    return f'''
double {func_name}(const double* xx, const double* p){{
   const std::vector<double> t {{
   {', '.join([str(e) for e in tck[0].tolist()])}
   }};
   const std::vector<double> c {{
   {', '.join([str(e) for e in tck[1].tolist()])}
   }};
   int k = {tck[2]};
   return bspline(xx[0],k,t,c);
}};

'''


def get_commonFun():
    return '''
double deBoor(int l, double x, int k,
      const std::vector<double>& t, const std::vector<double>& c) {
   // Evaluates S(x) = Sum c_i * B_i,k (x)

   // Arguments
   // ---------
   // l: Index of knot interval that contains x
   // x: Position
   // t: Array of knot positions
   // c: Array of control points
   // k: Degree of B-spline

   std::vector<double> d(k+1);
   for( int j = 0; j < k+1; ++j ) {
      d[j] = c[j+l-k];
   }

   for( int r = 1; r < k+1; ++r ) {
      for( int j = k; j > r-1; --j ) {
         double alpha = (x-t[j+l-k]) / (t[j+1+l-r] - t[j+l-k]);
         d[j] = (1.0 - alpha) * d[j-1] + alpha * d[j];
      }
   }

   return d[k];
}

double bspline(double x, int k,
      const std::vector<double>& t, const std::vector<double>& c) {
   int n = t.size() - k - 1;
   // find index l in range [k,n) what t[l] <= x < t[l+1]
   int l = k+1;
   while( x >= t[l] && l != n ) {
      l += 1;
   }
   return deBoor(l-1, x, k, t, c);
}

'''


def get_rootFun(func_name, *dat):
    # test root-script for debug
    fX = dat[0]
    fY = dat[1]
    dY = dat[2]
    return f'''
void test_{func_name}() {{
   gROOT->Reset();
   gStyle->SetOptStat(0);

   const std::vector<double> fX {{
   {', '.join([str(e) for e in fX])}
   }};
   const std::vector<double> fY {{
   {', '.join([str(e) for e in fY])}
   }};
   const std::vector<double> dY {{
   {', '.join([str(e) for e in dY])}
   }};
   int n = fX.size();
   const std::vector<double> dX(n);
   TGraphErrors* gr = new TGraphErrors(n,fX.data(),fY.data(),
         dX.data(),dY.data());

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,500);
   c1->cd();
   c1->SetGrid();

   double xmin = {fX[0]-(fX[1]-fX[0])/2};
   double xmax = {fX[-1]+(fX[1]-fX[0])/2};
   TF1* Fbspl = new TF1("Fbspl",{func_name},xmin,xmax,0);
   Fbspl->SetTitle("Smooth spline: {func_name}");
   Fbspl->SetLineWidth(2);

   gr->SetTitle(Fbspl->GetTitle());
   gr->Draw("AP");
   Fbspl->Draw("SAME");
   c1->Update();
}}
'''


def write_Bcpp(tck, func_name, *dat):
    with open(f'test_{func_name}.cc', 'w') as fw:
        fw.write(get_commonFun())
        fw.write(get_Fun(tck, func_name))
        fw.write(get_rootFun(func_name, *dat))


if __name__ == "__main__":
    # change the fontset so that
    # \varepsilon looks different than \epsilon
    plt.rcParams['mathtext.fontset'] = 'cm'

    # global params:
    plt.rcParams['legend.fontsize'] = 'large'
    plt.rcParams['axes.labelsize'] = 'large'

    # tests
    #  bes3_data('Kp_2021', 0.005, 2)
    #  bes3_data('Km_2021', 0.005, 0)
    #  bes3_data('Kp_2012', 0.005, 19.)
    #  bes3_data('Pip_2021', 0.005, 8.)
    #  bes3_data('Pim_2012', 0.005, 8.)
    #  bes3_data('Pim_2009', 0.005, 5.5)
    #  exit()

    # Correction functions
    show_pict = True
    with open(f'Bsp_K_Pi.cc', 'w') as fw:
        header = '// {{{1 Bsp_K... functions\n'
        fw.write(header)
        tsc = 0
        ret = bes3_data('Kp_2021', 0.005, tsc, show_pict)
        fw.write(ret)
        ret = bes3_data('Km_2021', 0.005, tsc, show_pict)
        fw.write(ret)
        ret = bes3_data('Kp_2012', 0.005, 19., show_pict)
        fw.write(ret)
        ret = bes3_data('Km_2012', 0.005, tsc, show_pict)
        fw.write(ret)
        ret = bes3_data('Kp_2009', 0.005, tsc, show_pict)
        fw.write(ret)
        ret = bes3_data('Km_2009', 0.005, tsc, show_pict)
        fw.write(ret)
        #
        tsc = 1
        header = '// {{{1 Bsp_Pi... functions\n'
        fw.write(header)
        ret = bes3_data('Pip_2021', 0.005, tsc, show_pict)
        fw.write(ret)
        ret = bes3_data('Pim_2021', 0.003, tsc, show_pict)
        fw.write(ret)
        ret = bes3_data('Pip_2012', 0.003, tsc, show_pict)
        fw.write(ret)
        ret = bes3_data('Pim_2012', 0.005, tsc, show_pict)
        fw.write(ret)
        ret = bes3_data('Pip_2009', 0.003, tsc, show_pict)
        fw.write(ret)
        ret = bes3_data('Pim_2009', 0.005, 5.5, show_pict)
        fw.write(ret)
