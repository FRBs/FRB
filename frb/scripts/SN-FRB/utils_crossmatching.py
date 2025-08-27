""" Utilities for TNS Transient cross-matching with FRBs
@author: Yuxin Dong
last edited: May 13, 2025 """

import numpy as np
from scipy.stats import chi2
from astropy import units as u
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def cov_matrix(a, b, theta):

    """
    Calculate the covariance matrix for a 2D ellipse.
    
    Parameters:
    a (float): Semi-major axis of the ellipse
    b (float): Semi-minor axis of the ellipse
    theta (float): Position angle of the ellipse in degrees
    
    Returns:
    np.ndarray: 2x2 covariance matrix
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


def mahalanobis_distance(point, frbcenter, cov_matrix):
    """
    Calculate the Mahalanobis distance between a given point and the center.
    effectively the Z-score in 1D, and it can be a proxy for how many sigma away you are from the mean.

    Parameters:
    point (array-like): The point for which to calculate the distance (should be [RA, Dec]).
    frbcenter (array-like): The center point, usually the mean [RA, Dec].
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


def percentile(mahalanobis_distance, df=2):
    p_value = chi2.cdf(mahalanobis_distance**2, df)
    return p_value


def gauss_contour(frbcenter, cov_matrix, semi_major, transient_name,
                  transient_position=None, levels=[0.68, 0.95, 0.99]):
    
    """
    Plot Gaussian contours around a given FRB center based on its covariance matrix.
    
    This function computes and visualizes the Gaussian contours (ellipses) 
    that represent the uncertainty in the position of the FRB as defined by 
    its covariance matrix. Then, it plots the position of the transient along with 
    the Mahalanobis distance from the FRB center.
    
    Parameters:
    frbcenter (Astropy SkyCoord): The central position of the FRB.
    cov_matrix (array-like): The covariance matrix representing the uncertainties in the FRB's position.
    semi_major (float): The semi-major axis of the Gaussian ellipse in degrees.
    transient_name (str): The name of the transient source.
    transient_position (Astropy SkyCoord): Transient position.
    levels: List of confidence levels for the contours (default: [0.68, 0.95, 0.99]).
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
        ellipse = Ellipse(
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
    #plt.savefig(f'{transient_name}_gaussian_map.png')
