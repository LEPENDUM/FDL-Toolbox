# FDL toolbox

## Introduction

This is the FDL toolbox providing matlab functions for Light Fields processing in the Fourier Disparity Layer representation.

The toolbox includes:

1. Calibration of sub-aperture image input : Determines view positions on the camera plane and a set of disparity values.

2. FDL construction algorithm.

3. Interactive Light Field Rendering application from a FDL model.

## Requirements

The toolbox has been developped using Matlab 2017a.

Matlab's Parallel Processing Toolbox and a CUDA GPU are required.

## Demo

A demo of the complete processing chain is provided in ./Demo/TestFDL.m