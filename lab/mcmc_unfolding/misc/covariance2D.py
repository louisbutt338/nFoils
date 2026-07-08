# Boilerplate imports
from numpy import array as ary; from numpy import log as ln
from numpy import pi, sqrt, arccos, arcsin
import numpy as np
from matplotlib import pyplot as plt
# Linear algebra functions
from numpy.linalg import inv, det, eig
from matplotlib.patches import Ellipse # for plotting ellipse
import matplotlib as mpl
from matplotlib.colorbar import ColorbarBase
from matplotlib import colors
# Stochastic distribution functions
import numpy.random as rd
#plt.rcParams['pgf.preamble'] = r'\usepackage{amsmath}'  # for \text command as well
#plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'  # for \text command
#plt.rcParams['text.usetex'] = True  # use TeX to render everything
#plt.rcParams['font.family'] = 'STIXGeneral'  # correct font!
multi_norm = rd.multivariate_normal

RED = [1.0,0.0,0.0]
ORANGE = [1.0,0.35,0.0]
GREEN = [0.0,1.0,0.0]
blue_green = GREEN+ary([0.0, 0.0, 0.2])
CHOSEN_COLOR = blue_green
WHITE = [1.0, 1.0, 1.0]

# program to generate a covariance matrix where the variance values are fixed at [2,2].
def generate_cov(off_diag):
    cov = ary([[2.0, off_diag], [off_diag, 1.0]])
    return cov
center = (2, 1)


def chi2_to_colour(values, max_chi2):
    """values: list of scalar (chi^2) values that we want to transform (independently of each other).
    max_chi2: maximum value of chi^2 before the colour saturates"""
    values = np.array(values)
    alpha = np.clip(sqrt(values)/sqrt(max_chi2), 0, 1)
    transparent_red = np.append(CHOSEN_COLOR, 0.0)
    template = np.tile(transparent_red, (len(values), 1))
    pixel_colours = template.T
    pixel_colours[-1] += alpha
    pixel_colours = pixel_colours.T
    return pixel_colours


if __name__=='__main__':
    #mpl.rcParams['text.usetex'] = True
    # covar_range = np.linspace(-1.3, 1.3, 3) # plotting the actual covariance circle
    (fig, (axes)) = plt.subplots(1, 3, gridspec_kw=dict(
        # width_ratios=[0.33, 0.33, 0.33 ,0.01]
        # height_ratios=[0.97, 0.03],
        ))
    for plot_index, cov_value in enumerate([-1.3, 0, 1.3]):
    # create a scatter of dots for plotting
        ax = axes[plot_index]
        rd.seed(0)
        # generate the covariance matrix, and evaluate the shape of the error ellipse
        cov = generate_cov(cov_value)
        print(cov)
        (minor, major), eig_mat = eig(inv(cov) * det(cov))
        print(major, minor)
        mat = inv(eig_mat)
        # ignore the arccos becuase it will always return a non-negative value.
        orientation = np.mean([arcsin(mat[0,1]), -arcsin(mat[1,0])])*np.sign(mat[0,0])
        def get_chi2(xvec_yvec_matrix):
            """xvec_yvec_matrix: shape (2, [N*M]) where the second axis is optional."""
            difference_in_N_space = ((xvec_yvec_matrix).T-center)
            return [(diff_vec @ inv(cov) @ diff_vec) for diff_vec in np.atleast_2d(difference_in_N_space)]
        # plotting
        ellipse = Ellipse(center, # centered at the origin
                    2*sqrt(major),
                    2*sqrt(minor),
                    angle=np.rad2deg(orientation),
                    fill=False,
                    #label=r"Area in which $\chi^2 \le 1$",
                    label=r"Area in which $\chi^2$ ≤ 1",
                    )
        # DUUUUDE I got the width=2*sqrt(major)/sqrt(det(inv(cov))) equation by trial and error LMAO
        ax.add_patch(ellipse)
        xrange, yrange = [-1.5+2, 1.5+2], [-1.5+1, 1.5+1]
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.set_xlabel(r"$x_1$")
        ax.set_ylabel(r"$x_2$")
        if plot_index==0:
            ax.set_title(r"mean=($2,1$)"" cov(x)= ( 2  -1.3, -1.3  1 ) ")
        elif plot_index==1:
            ax.set_title(r"mean=($2,1$) ""cov(x)= ( 2  0,  0  1 ) ")
        elif plot_index==2:
            ax.set_title(r"mean=($2,1$)"" cov(x)= ( 2  1.3 , 1.3  1) ")

        ax.set_aspect(1.0)
        sample_points = multi_norm(center, cov, size=5000)
        ax.scatter(*sample_points.T, marker='+', alpha=0.2, color='C0', label="samples") # scatter plot approach
        ax.legend()
    fig.set_size_inches(10,4.5)
    fig.set_tight_layout(True)
    plt.savefig("covariance2D.pdf")
    # plt.show()