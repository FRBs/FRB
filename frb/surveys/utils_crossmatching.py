""" Utilities for TNS Transient cross-matching with FRBs
@author: Yuxin Dong
last edited: May 13, 2025 """

import numpy as np
import pandas as pd
import json
import time
import requests
import os
import sys
from collections import OrderedDict

from astropy.coordinates import SkyCoord
import ligo.skymap.plot # KEEP: needed for projections in matplotlib
from astropy import units as u

from scipy.stats import chi2
from astropy.io import ascii

from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import matplotlib




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