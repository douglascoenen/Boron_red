"""Helper functions to declutter the main file."""
import numpy as np
import scipy 
import scipy.stats as stats
import matplotlib.pyplot as plt
import uncertainties as unc
from uncertainties import unumpy as unp


class Bunch(dict):
    """
    Create a subclass of dictionary that allows the access with "."
    
    Taken from Oscar's Branson cbsyst, itself simplified from Bunch library
    """
    def __init__(self, *args, **kwds):
       super(Bunch, self).__init__(*args, **kwds)
       self.__dict__ = self

def exp_fun(x, a, b, c): 
    """Exponential function with intercept for inaccuracy correction."""
    return a*x**b + c

def float_check(s):
    """
    Filter the .exp raw Neptune files.
    
    Checks if the first character of a string is a number. Which would indicate
    that it is a row of data. To use in conjoncture with filter()
    
    Parameters
    ----------
    string : str
        line of data.

    Returns
    -------
    True or False
    """
    try:
        # Try to convert the first character to a float.
        float(s[0])
        return True
    except ValueError:  # If the first character is not a number,
        return False   # it is part of the block header or empty line

def getloc(ar, val):
    """
    Find the closest value in an array.

    Parameters
    ----------
    ar : list or np.array
        array of values to find the value
    val : float or int
        value to look for in the array

    Returns
    -------
    int
        position of the closest value in the array

    """
    return np.argmin(np.abs(ar - val))

def predband(x, xd, yd, p, func, conf=0.95):
    """
    Calculate the prediction bands of a regression function.
    
    Parameters
    ----------
    x : array
        "x" axis requested points .
    xd : array
        x data to be fitted in the model.
    yd : array
        y data to compare the model to.
    p : array
        fitted parametres (a,b,c) .
    func : function
        function to fit the data in.
    conf : float, optional
        confidence interval. The default is 0.95.

    Returns
    -------
    lpb : array
        lower prediction band.
    upb : array
        upper prediction band.

    """
    alpha = 1.0 - conf    # significance
    N = xd.size          # data sample size
    var_n = len(p)  # number of parameters
    # Quantile of Student's t distribution for p=(1-alpha/2)
    q = scipy.stats.t.ppf(1.0 - alpha / 2.0, N - var_n)
    # Stdev of an individual measurement
    se = np.sqrt(1. / (N - var_n) *
                 np.sum((yd - func(xd, *p)) ** 2))
    # Auxiliary definitions
    sx = (x - xd.mean()) ** 2
    # sxd = np.sum((xd - xd.mean()) ** 2)
    sxd = (N - 1)*xd.std()**2
    # Predicted values (best-fit model)
    yp = func(x, *p)
    # Prediction band formula
    dy = q * se * np.sqrt(1.0 + (1.0/N) + (sx/sxd))
    # Upper and lower prediction bands.
    lpb, upb = yp - dy, yp + dy

    return lpb, upb

def mcerr(x, xsd, y, ysd, func):
    """Propagate uncertainties"""
    n = int(1e5)
    xm = np.random.normal(x, xsd, n)
    ym = np.random.normal(y, ysd, n)
    xy = func(xm, ym)
    return np.median(xy), np.std(xy)

def stdmask(df, sd):
    """Mask outliers of DataFrame above or below n SD"""
    mask = np.abs(scipy.stats.zscore(df)) > sd
    df = df.where(~mask)
    return df

def stdzrem(ar, sd):
    """Remove outliers from an array (2sd)."""
    arz = np.abs(scipy.stats.zscore(ar))
    return ar[arz < sd]

def stderem(are, sd):
        
    are.loc[abs(are - are.mean()) > are.std()*sd] = np.nan
    return are

def matsdrem(df, sd):
    """David Matlab way, remove whole row"""
    mask = np.abs(scipy.stats.zscore(df)) > sd
    ndf = df.copy()
    for i in range(len(df)):
        if True in mask.iloc[i, :]:
            ndf.iloc[i, :] = np.full(len(df.columns), np.nan)
    
    return ndf

def cooksdistance(x, y, mfun, popt):
    """
    Calculate the Cook's distance of a regression (linear and non-linear).
    
    This code was adapted from  the Yellowbrick library to make it 
    fit a non-linear equations. (https://github.com/DistrictDataLabs/yellowbrick)

    Parameters
    ----------
    x : np.array
        Values to fit into the model.
    y : np.array
        Real data.
    mfun : callable
        function to fit.
    popt : list
        parametres of the regression.

    Returns
    -------
    distance : np.array
        Array containing the cook's distance values for each data point.
    influence_threshold : float
        Threshold value after which a point is considered an outlier.
    r2 : float
        Correlation coefficient of the fit.

    """
    
    fy = mfun(x, *popt) # fitted y values

    mX = np.vstack(x)
    # Calculate leverage values
    # plt.plot(mX)
    # plt.plot(mX.T)
    a = np.dot(mX.T, mX)
    # print(mX, mX.T)
    # print(a)
    b = np.linalg.inv(np.dot(mX.T, mX))
    b2 = 1/ np.dot(mX.T, mX)
    
    # print(b, b2)
    c = np.dot(np.linalg.inv(np.dot(mX.T, mX)),mX.T)
    # print(np.shape(c))
    H = np.dot(mX, c)
    # print(np.shape(H))
    hii = np.diag(H)
    leverage = (mX * np.linalg.pinv(mX).T).flatten()
    # print(hii)
    

    # Compute the rank and the degrees of freedom of the regression
    rank = np.linalg.matrix_rank(mX)
    
    thr = 2 * (rank/mX.shape[0])  # Threshold for too high leverage (2*(p/n))
    plt.plot(hii)
    plt.axhline(thr)
    # print(thr)
    df = mX.shape[0] - rank
    # Compute the MSE from the residuals
    residuals = y - fy

    mse = np.dot(residuals, residuals) / df
    
    residuals_studentized = residuals / np.sqrt(mse) / np.sqrt(1 - leverage)
    
    # Calculate Cook's distance
    distance = residuals_studentized ** 2 / mX.shape[1]
    distance *= leverage / (1 - leverage)
    
    nd = residuals**2 / (rank*mse)
    nd *= leverage / (1-leverage)**2
    
    # print(distance - nd)

    # Compute the influence threshold rule of thumb
    # p_values = stats.f.sf(distance, mX.shape[1], df)
    # influence_threshold = 4 / (mX.shape[0]-len(popt))
    # Other threshold is the 50th percentile of the F function.
    
    pvals = stats.f.sf(nd, rank, df)
    # print(pvals)
    # influence_threshold = 4 / (len(distance) - rank - 1)
    influence_threshold = 4 / len(distance)  # Find source ?? 
    
    r2 = 1.0-(sum((y-mfun(x, *popt))**2) /
                              ((len(y)-1.0)*np.var(y, ddof=1)))
    
    return distance, influence_threshold, r2
    
def mancooks():
    # Try to make the normal manual version of the Cook's distance
    return None
def DFFITS():
    # Try to program DFFITS -> simples than Cooks distance
    # Cook's distance and DFFITS are the same conceptually
    return None


def an_plot(sig, i, n, sn):
    """Plot a single analysis."""
    fig, ax = plt.subplots()
    c1 = "blue"
    c2 = "red"   
    c3 = "magenta"
    
    rr = sig["11B"] / sig["10B"]
    leng = len(rr)
    ax.scatter(range(leng), rr, c=c1, edgecolor="k", s=20)
    ax.set_ylim(np.nanmin(rr.values)*0.95, np.nanmax(rr.values)*1.05)
    ax.set_xlabel("Integration cycles")
    ax.set_ylabel(r"raw $^{11}$B/$^{10}$B")
    # ax.set_ylabel(r"raw $^{138}$Ba/$^{135}$Ba")
    ax.tick_params(axis="y", colors=c1)
    ax.yaxis.label.set_color(c1)
    # Plot the raw 11B voltage to detect spikes
    ax2 = ax.twinx()
    ax2.plot(range(leng), sig["11B"], c=c2)
    # ax2.plot(range(len(sig["138Ba"])), sig["138Ba"], c=c2)
    yL = np.nanmax(sig["11B"].values)
    # yL = np.max(sig["138Ba"].values)
    ax2.set_ylim(0, yL*1.2)
    ax2.set_ylabel(r"$^{11}$B voltage (V)")
    # ax2.set_ylabel(r"$^{138}$Ba voltage (V)")
    # Nicer colour matching
    ax2.spines['left'].set_color(c1)
    ax2.spines["right"].set_color(c2)
    ax2.tick_params(axis="y", colors=c2)
    ax2.yaxis.label.set_color(c2)
    ax2.set_title("Click twice or press enter to use the whole analysis",
                  weight="bold")
    
    ax3 = ax.twinx()
    ax3.spines.right.set_position(("axes", 1.25))
    b11_ca10 = sig["11B"] / sig["10.035"]
    b11_ca9 = sig["11B"] / sig["9.979"]
    # ax3.plot(range(leng), b11_ca10, c=c3) 
    ax3.plot(range(leng), b11_ca9, c=c3, alpha=0.7)
    ax3.spines['left'].set_color(c1)
    ax3.spines["right"].set_color(c3)
    ax3.tick_params(axis="y", colors=c3)
    ax3.yaxis.label.set_color(c3)
    ax3.set_ylabel(r"$^{11}$B / 9.979")

    # Display the sample number and what is left to cut
    ax2.text(0, yL*0.05, f"Sample : {sn}")
    ax2.text(0, yL*1.1, fr"Analysis {i+1} of {n}")
    # plt.tight_layout()  # Magic to look nice
    
    return fig


def polynomfit(inp, *args):
    """Polynomial fitting function, from Jie script"""
    x = inp
    res = 0
    for order in range(len(args)):
        res += args[order] * x**order
    return res
    