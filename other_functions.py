#13/01/2025

#other useful functions to use with the EOS and EOS data
import numpy as np
import pandas as pd
from scipy.spatial import distance
import warnings
from scipy import optimize
from symbols import *

#find multiple roots of a function between a specific bracket (univariate function only)
def multi_root(f: callable = None, bracket = None, args: tuple = (), n: int = 1000) -> np.ndarray:
    """ Find all roots of f in `bracket`, given that resolution `n` covers the sign change.
    Fine-grained root finding is performed with `scipy.optimize.root_scalar`.
    Parameters
    ----------
    f: Callable
        Function to be evaluated
    bracket: Sequence of two floats
        Specifies interval within which roots are searched in log space.
    args: Iterable, optional
        Iterable passed to `f` for evaluation
    n: int
        Number of points sampled equidistantly from bracket to evaluate `f`.
        Resolution has to be high enough to cover sign changes of all roots but not finer than that.
        Actual roots are found using `scipy.optimize.root_scalar`.
    Returns
    -------
    roots: np.ndarray
        Array containing all unique roots that were found in `bracket`.
    """
    # Evaluate function in given bracket

    # x = np.linspace(*bracket, int(n))
    x = np.logspace(*bracket, int(n))

    y = f(x, *args)

    # Find where adjacent signs are not equal
    sign_changes = np.where(np.sign(y[:-1]) != np.sign(y[1:]))[0]

    # Find roots around sign changes
    root_finders = (
        optimize.root_scalar(
            method="bisect",
            f=f,
            args=args,
            bracket=(x[s], x[s+1]),
            xtol=1e-12, rtol=1e-12
        )
        for s in sign_changes
    )

    roots = np.array([
        r.root if r.converged else np.nan
        for r in root_finders
    ])


    if np.any(np.isnan(roots)):
        warnings.warn("Not all root finders converged for estimated brackets! Maybe increase resolution `n`.")
        roots = roots[~np.isnan(roots)]

    roots_unique = np.unique(roots)
    if len(roots_unique) != len(roots):
        warnings.warn("One root was found multiple times. "
                    "Try to increase or decrease resolution `n` to see if this warning disappears.")

    #check if the roots correspond to actual roots or a discontinuity
    #if inputting the root into the function results in a distance higher than 1 then assume its not a root and assign nan
    roots_unique[abs(f(roots_unique, *args)) > 1e-3] = np.nan 
    roots_unique = roots_unique[~np.isnan(roots_unique)]

    # sort the roots in ascending order
    roots_unique = np.sort(roots_unique)


    return roots_unique   

#choose the best root from a list by knowing the experimental value
def choose_root(roots: list, experimental_value=0):
    roots = np.array(roots)

    difference = np.abs(experimental_value-roots)

    try:
        best_index = np.nanargmin(difference)
        best_root = roots[best_index]

    except ValueError:
        best_root = np.nan

    return best_root

#see zandbox for trial and error
def calc_rho_inverse_distance(df_missing_rho, df_experimental, pascals=1e4, neighbours=3):
    #restrict experimental dataframe to Paper, P, T, rho and deep copy it
    df_experimental = df_experimental.copy(deep=True)[["Paper", "P", "T", "rho"]]

    #restrict the other dataframes to Paper, P,T and deep copy it
    df_missing_rho = df_missing_rho.copy(deep=True)[["Paper", "P", "T"]]

    #create og index columns
    df_missing_rho["index"] = df_missing_rho.index
    df_experimental["index"] = df_experimental.index
    
    #merge the two dataframes on pressure and temperature
    df = pd.merge(df_experimental, df_missing_rho, on=["Paper","P", "T"], how="outer")

    #sort by temperature and then by pressure
    df = df.sort_values(by=["T","P"])
    #delete pressure and temperature duplicates
    # df = df.drop_duplicates(subset=["P", "T"])

    #reset the indices
    df = df.reset_index()

    #create existance column
    df["exist"] = df["rho"].notna()

    #reduce the pressure by the pascals
    df["P"] = df["P"]/pascals

    #calculate distances between points
    distances = distance.cdist(df[["T","P"]], df[["T","P"]], "euclidean")

    #multiply the distances by the truth vector
    distances_true = (distances.T * np.array(df["exist"])).T

    #make the distances a dataframe
    df_distances = pd.DataFrame(data=distances_true)

    #get the point numbers for which you would like to linearly interpolate
    point_numbers = df.index[df["exist"]==False]

    #list to keep track of interpolation values
    rho_interp = []

    #for unknown point
    for i in point_numbers:
        distance_column = df_distances[i]

        #get the smallest neighbours nonzero values
        smallest_n = distance_column.loc[distance_column!=0].nsmallest(neighbours)

        #calculate the inverse distance density and add it to the list
        rho_interp +=[sum(df.loc[smallest_n.index, "rho"]*1/smallest_n.values)/sum(1/smallest_n.values)]

    #add the interpolated values to the dataframe
    df.loc[point_numbers,["rho"]] = rho_interp

    #return the pressure back to normal
    df["P"] = df["P"]*pascals

    #drop all the rows where the index_y is Nan
    df = df[df["index_y"].notna()]

    # set index to the index of the dataframe which is missing density
    df = df.set_index("index_y").sort_index().reset_index(drop=True)

    return df[["Paper", "P", "T", "rho"]]   #return the  dataframe but only Paper, P, T, rho

def calc_Psat(fluid, EOS, temperature_list: list):
    """
    Calculate the saturation pressure of a fluid at given temperatures.

    Parameters:
    fluid: The fluid object containing properties like critical and triple point temperatures.
    EOS: The equation of state object used for pressure and fugacity calculations.
    temperature_list (list): A list of temperatures at which to calculate saturation pressures.

    Returns:
    list: A list of calculated saturation pressures corresponding to the input temperatures.
    """
    
    pressure_equation = sp.utilities.lambdify((rho, T, P), P - EOS.pressure)
    fugacity_coefficient_function = sp.utilities.lambdify((rho, T), EOS.fugacity_coefficient)

    def get_pressure_guess(temperature):
        """
        Get an initial pressure guess for a specified temperature.

        Parameters:
        temperature: The temperature for which to estimate the pressure.

        Returns:
        float: An initial pressure guess.
        """
        pressure_guess = fluid.P_t# + (fluid.P_c - fluid.P_t) * (temperature - fluid.T_t) / (fluid.T_c - fluid.T_t)
        n_roots = 0
        while n_roots < 3:
            pressure_guess += 1e2
            density_roots = multi_root(pressure_equation, [1e-4, 1e4], args=(temperature, pressure_guess))
            n_roots = len(density_roots)
        return pressure_guess

    def calc_saturation_pressure(pressure, temperature):
        """
        Calculate the saturation pressure from temperature and initial pressure guess.

        Parameters:
        pressure: The initial pressure guess.
        temperature: The temperature for which to calculate saturation pressure.

        Returns:
        float: The result of the fugacity ratio equation.
        """
        density_roots = multi_root(pressure_equation, [1e-4, 1e4], args=(temperature, pressure))
        fugacities = fugacity_coefficient_function(density_roots, temperature)
        fugacities = fugacities[fugacities != max(fugacities)]
        eq = fugacities[0] / fugacities[1] - 1
        return eq

    P_sat = []
    pressure_guess = get_pressure_guess(temperature_list[0])

    for temperature in temperature_list:
        if temperature > fluid.T_c or temperature < fluid.T_t:
            print("Temperature is outside the critical and triple point bounds")
            break
        try:
            sol = optimize.root(calc_saturation_pressure, x0=[pressure_guess], args=(temperature,), method='broyden1', options={"xtol": 1e2})
            solution = sol.x
            pressure_guess = sol.x[0]
        except IndexError:
            print("IndexError: no solution has been found due to the root finding method reaching a point with less than two density roots")
            print("Assigning 0 to return value")
            solution = [0]
        P_sat += list(solution)

    return P_sat


