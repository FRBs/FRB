""" Utilities for TNS Transient cross-matching with FRBs
@author: Yuxin Dong
last edited: May 13, 2025 """

import numpy as np
from numpy.typing import ArrayLike
from scipy.stats import chi2
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib


def cov_matrix(a: float, b: float, theta: float) -> np.ndarray:
    """
    Calculate the covariance matrix for a 2D ellipse.

    Args:
        a (float): Semi-major axis of the ellipse
        b (float): Semi-minor axis of the ellipse
        theta (float): Position angle of the ellipse in degrees

    Returns:
        numpy.ndarray: 2x2 covariance matrix
    """

    # Convert theta to radians
    theta_rad = np.radians(theta + 90)  # Adjust angle for correct orientation
    # Rotation matrix for an ellipse with position angle
    R = np.array([[np.cos(theta_rad), -np.sin(theta_rad)], 
                  [np.sin(theta_rad), np.cos(theta_rad)]])
    # Diagonal matrix with square of semi-major and semi-minor axes
    D = np.diag([a**2, b**2])
    # Covariance matrix
    covariance_matrix = R @ D @ R.T
    return covariance_matrix


def mahalanobis_distance(point: ArrayLike, frbcenter: ArrayLike,
                         cov_matrix: ArrayLike) -> float:
    """
    Calculate the Mahalanobis distance between a given point and the center.
    effectively the Z-score in 1D, and it can be a proxy for how many sigma away you are from the mean.

    Args:
        point (array-like): The point for which to calculate the distance (should be [RA, Dec]) in degrees.
        frbcenter (array-like): The center point, usually the mean [RA, Dec], in degrees.
        cov_matrix (array-like): The covariance matrix of the data.

    Returns:
        float: The Mahalanobis distance.
    """

    # the frbcenter in this case is the 'mean'
    diff = np.array(point) - np.array(frbcenter)
    
    # correction term for spherical geometry
    correction = diff[0] * (1 / np.cos(np.radians(point[1])))
    diff_corrected = np.array([correction, diff[1]])
    inv_cov_matrix = np.linalg.inv(cov_matrix)
    md = np.sqrt(diff_corrected.T @ inv_cov_matrix @ diff_corrected)
    
    return md


def percentile(mahalanobis_distance: float, df: int = 2) -> float:
    """
    Convert a Mahalanobis distance into the cumulative probability of a chi-squared distribution.

    Args:
        mahalanobis_distance (float): The Mahalanobis distance.
        df (int, optional): Degrees of freedom of the chi-squared distribution.

    Returns:
        float: Probability that a point drawn from the distribution lies within
        ``mahalanobis_distance`` of the center.
    """
    p_value = chi2.cdf(mahalanobis_distance**2, df)
    return p_value


def gauss_contour(frbcenter: SkyCoord, cov_matrix: ArrayLike,
                  semi_major: float, transient_name: str,
                  transient_position: SkyCoord | None = None,
                  levels: list[float] = [0.68, 0.95, 0.99],
                  save_fig: bool = False) -> None:
    
    """
    Plot Gaussian contours around a given FRB center based on its covariance matrix.
    
    This function computes and visualizes the Gaussian contours (ellipses) 
    that represent the uncertainty in the position of the FRB as defined by 
    its covariance matrix. Then, it plots the position of the transient along with 
    the Mahalanobis distance from the FRB center.

    Args:
        frbcenter (astropy.coordinates.SkyCoord): The central position of the FRB.
        cov_matrix (array-like): The covariance matrix representing the uncertainties in the FRB's position.
        semi_major (float): The semi-major axis of the Gaussian ellipse in degrees.
            Currently not used.
        transient_name (str): The name of the transient source.
        transient_position (astropy.coordinates.SkyCoord, optional): Transient position.
        levels (list of float, optional): Confidence levels for the contours (default: [0.68, 0.95, 0.99]).
        save_fig (bool, optional): set flag to true to save the figure
            as ``<transient_name>_gaussian_map.png``.

    Returns:
        None
    """
    
    eigvals, eigvecs = np.linalg.eigh(cov_matrix)
    order = eigvals.argsort()[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    
    width, height = 2 * np.sqrt(eigvals)
    theta = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))
    threshold = chi2.ppf(0.9999, df=2)
    
    # 2 arcmins size, quite arbitrary so can be changed
    buffer = 2 * u.arcsec
    if transient_position is not None:
        size = (frbcenter.separation(transient_position)+buffer).to(u.arcmin)
    else:
        size = 10*u.arcmin #semi_major*60*30

    fig = plt.figure(figsize=(6, 6))
    ax = plt.axes(
    projection='astro zoom',
    center=frbcenter,
    radius=size)

    # Plot FRB center
    ax.plot(frbcenter.ra.deg, frbcenter.dec.deg, 
        transform=ax.get_transform('world'), color='black', marker='x', 
        markersize=5)
    
    cmap = matplotlib.colormaps.get_cmap('viridis')

    for i, level in enumerate(levels):
        chi_square_val = chi2.ppf(level, 2)
        color = cmap(i / len(levels))
        ellipse = matplotlib.patches.Ellipse(
            xy=(frbcenter.ra.deg, frbcenter.dec.deg),
            width=width * np.sqrt(chi_square_val),
            height=height * np.sqrt(chi_square_val),
            angle=theta,
            edgecolor=color,
            fc='None',
            lw=2,
            label=f'{int(level*100)}%',
            transform=ax.get_transform('world')
        )
        ax.add_patch(ellipse)
    
    if transient_position is not None:
        ax.plot(transient_position.ra.deg, transient_position.dec.deg, 
                transform=ax.get_transform('world'), 
                color='#fb5607', marker='o', label=transient_name)
        
        
        md = mahalanobis_distance([transient_position.ra.deg, transient_position.dec.deg],
                                  [frbcenter.ra.deg, frbcenter.dec.deg], cov_matrix)
        pt = percentile(md)
        ax.annotate(f'{(1-pt)*100:.3f}%', xy=[transient_position.ra.deg, transient_position.dec.deg],
                    xytext=(12, 12), textcoords='offset points')     
        
    ax.set_xlabel('RA (J2000)', fontsize=14)
    ax.set_ylabel('DEC (J2000)', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(False)
    plt.tight_layout()
    if save_fig:
        plt.savefig(f'{transient_name}_gaussian_map.png')
