% Contents of the Complex Function Explorer (former Phase Plot) package
%
% Complex Function Explorer - Visualization of complex functions
% Version 1.1, February 1, 2014
% Copyright(c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)
%
% For background information on phase plots see 
% E. Wegert, Visual Complex Functions, Birkhaeuser-Springer 2012
%
% -----------------------------------------------------------------------------
%
% Demonstrations
%
% PPDemo        - demonstration of phase plots 
% StartMe       - the same as PPDemo
%
% -----------------------------------------------------------------------------
%
% User interface
%
% FEXGUI        - call of graphical user interface
%
% -----------------------------------------------------------------------------
%
% Main routines 
%
% PhasePlot     - phase plot of complex function
% PPOnSphere    - phase plot of function on Riemann sphere
% PPOnSphLow    - stereographic projection of phase plot from lower hemisphere
% PPOnSphUp     - stereographic projection of phase plot from upper hemisphere
%
% -----------------------------------------------------------------------------
%
% auxiliary functions
%
% MyFunction    - example of matlab routine (singular inner function) 
%
% colscheme     - algorithms for converting complex numbers to colors 
% fldrect       - discrete field of complex numbers, domain of function 
% unitcirc      - vector with discrete points on complex unit circle
% zdomain       - the same as fldrect (former version)
% zplane        - discrete field of points in the entire complex plane
% zplanePP      - the same as zplane (to avoid confusion with existing m-file)
% 
% -----------------------------------------------------------------------------
%
% auxiliary functions, internal use
%
% brightenRGB   - modification of color scheme 
% gui_phaseplot - graphical user interface 
% pltphase      - phase plot
% stereoP2S     - stereographic projection from plane to sphere
% stereoS2P     - stereographic projection from sphere to plane
% unitcirc      - uniformly distributed points on the complex unit circle
%
% -----------------------------------------------------------------------------
%
% directories
%
% Examples      - contains various examples for use with PPGUI
% Figures       - default directory for saving figures (will be created ) 
%
% -----------------------------------------------------------------------------
%
% Contents      - this file

