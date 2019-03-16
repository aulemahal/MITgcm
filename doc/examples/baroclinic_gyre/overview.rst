Overview
--------

This example experiment demonstrates using the MITgcm to simulate a
baroclinic, wind-forced, ocean gyre circulation. The experiment is a
numerical rendition of the gyre circulation problem similar to the
problems described analytically by Stommel in 1966
:cite:`Stommel66` and numerically in Holland et. al
:cite:`Holland75`.
In this experiment the model is configured to represent a mid-latitude
enclosed sector of fluid on a sphere,
:math:`60^{\circ} \times 60^{\circ}` in lateral extent. The fluid is
:math:`2` km deep and is forced by a constant in time zonal wind
stress, :math:`\tau_{\lambda}`, that varies sinusoidally in the
north-south direction. Topologically the simulated domain is a sector
on a sphere and the coriolis parameter, :math:`f`, is defined
according to latitude, :math:`\varphi`

.. math::
   :label: eq_eg_fourlayer_fcori
   
   f(\varphi) = 2 \Omega \sin( \varphi )

with the rotation rate, :math:`\Omega` set to :math:`\frac{2 \pi}{86400{\rm s}}`.
The sinusoidal wind-stress variations are defined according to

.. math::
   :label: eq_eg_fourlayer_taux

   \tau_{\lambda}(\varphi) = \tau_{0}\sin(\pi \frac{\varphi}{L_{\varphi}})

where :math:`L_{\varphi}` is the lateral domain extent (:math:`60^{\circ}`) and :math:`\tau_0` is set to :math:`0.1N m^{-2}`.

:numref:`fig_eg_fourlayer_simulation_config` summarizes the configuration simulated. In contrast to the example in
section :numref:`sec_eg_baro`, the current experiment
simulates a spherical polar domain. As indicated by the axes in the
lower left of the figure the model code works internally in a locally
orthogonal coordinate :math:`(x,y,z)`. For this experiment description
the local orthogonal model coordinate :math:`(x,y,z)` is synonymous
with the coordinates :math:`(\lambda,\varphi,r)` shown in :numref:`sphere_coor`
The experiment has four levels in the vertical, each of equal
thickness, :math:`\Delta z = 500` m. Initially the fluid is stratified
with a reference potential temperature profile,
:math:`\theta_{250}=20^{\circ}` C, :math:`\theta_{750}=10^{\circ}` C,
:math:`\theta_{1250}=8^{\circ}` C, :math:`\theta_{1750}=6^{\circ}` C.
The equation of state used in this experiment is linear

.. math::
   :label: eq_eg_fourlayer_linear1_eos

   \rho = \rho_{0} ( 1 - \alpha_{\theta}\theta^{'} )

which is implemented in the model as a density anomaly equation

.. math::
   :label: eq_eg_fourlayer_linear1_eos_pert

   \rho^{'} = -\rho_{0}\alpha_{\theta}\theta^{'}

with :math:`\rho_{0}=999.8\,{\rm kg\,m}^{-3}` and :math:`\alpha_{\theta}=2\times10^{-4}\,{\rm degrees}^{-1}`. Integrated forward in this configuration the model state variable theta is
equivalent to either in-situ temperature, :math:`T`, or potential
temperature, :math:`\theta`. For consistency with later examples, in
which the equation of state is non-linear, we use :math:`\theta` to
represent temperature here. This is the quantity that is carried in the
model core equations.

  .. figure:: baroclinic_gyre/baroclinic_simulation_config.*
      :width: 80%
      :align: center
      :alt: baroclinic gyre configuration
      :name: fig_eg_fourlayer_simulation_config

      Schematic of simulation domain and wind-stress forcing function  for the four-layer gyre numerical experiment. The domain is enclosed by solid walls at :math:`0^{\circ}` E, :math:`60^{\circ}` E, :math:`0^{\circ}` N and :math:`60^{\circ}` N. An initial stratification is imposed by setting the potential temperature, :math:`\theta`, in each layer. The vertical spacing, :math:`\Delta z`, is constant and equal to 500 m.


Equations solved
----------------

For this problem the implicit free surface, HPE (see section :numref:`hydrostatic_quasihydrostatic_forms`)
form of the equations described in Marshall et. al
:cite:`marshall:97a` are employed. The flow is
three-dimensional with just temperature, :math:`\theta`, as an active
tracer. The equation of state is linear. A horizontal Laplacian operator
:math:`\nabla_{h}^2` provides viscous dissipation and provides a
diffusive sub-grid scale closure for the temperature equation. A
wind-stress momentum forcing is added to the momentum equation for the
zonal flow, :math:`u`. Other terms in the model are explicitly switched
off for this experiment configuration (see section :numref:`eg_baroclinic_code_config`). This yields an active set of
equations solved in this configuration, written in spherical polar
coordinates as follows

.. math::
   :label: eq_eg_fourlayer_model_equations
   \frac{Du}{Dt} - fv + \frac{1}{\rho}\frac{\partial p^{\prime}}{\partial \lambda} - 
     A_{h}\nabla_{h}^2u - A_{z}\frac{\partial^{2}u}{\partial z^{2}} 
   = \mathcal{F}_{\lambda}
   \\
   \frac{Dv}{Dt} + fu + \frac{1}{\rho}\frac{\partial p^{\prime}}{\partial \varphi} - 
     A_{h}\nabla_{h}^2v - A_{z}\frac{\partial^{2}v}{\partial z^{2}} 
   = 0

.. math::
   :label: eq_eg_fourlayer_example_continuity
   \frac{\partial \eta}{\partial t} + \frac{\partial H \widehat{u}}{\partial \lambda} +
   \frac{\partial H \widehat{v}}{\partial \varphi}
   = 0
   
.. math::
   :label: eq_eg_fourlayer_theta
   \frac{D\theta}{Dt} - K_{h}\nabla_{h}^2\theta  - K_{z}\frac{\partial^{2}\theta}{\partial z^{2}} 
   = 0
   \\
   p^{\prime}  = g\rho_{0} \eta + \int^{0}_{-z}\rho^{\prime} dz
   \\
   \rho^{\prime} = - \alpha_{\theta}\rho_{0}\theta^{\prime}
   \\
   {\cal F}_{\lambda} |_{s} = \frac{\tau_{\lambda}}{\rho_{0}\Delta z_{s}}
   \\
   {\cal F}_{\lambda} |_{i} = 0\end{aligned}

where :math:`u` and :math:`v` are the components of the horizontal flow
vector :math:`\vec{u}` on the sphere
(:math:`u=\dot{\lambda},v=\dot{\varphi}`). The terms
:math:`H\widehat{u}` and :math:`H\widehat{v}` are the components of the
vertical integral term given in equation :eq:`free-surface` and explained in more detail in
section :numref:`press_meth_linear`.
However, for the problem presented here, the continuity relation
(equation :eq:`eq_eg_fourlayer_example_continuity`
differs from the general form given in section :numref:`press_meth_linear`,
equation :eq:`linear-free-surface=P-E`,
because the source terms :math:`{\cal P}-{\cal E}+{\cal R}` are all
:math:`0`.

The pressure field, :math:`p^{\prime}`, is separated into a barotropic
part due to variations in sea-surface height, :math:`\eta`, and a
hydrostatic part due to variations in density, :math:`\rho^{\prime}`,
integrated through the water column.

The suffices :math:`{s},{i}` indicate surface layer and the interior of
the domain. The windstress forcing, :math:`{\cal F}_{\lambda}`, is
applied in the surface layer by a source term in the zonal momentum
equation. In the ocean interior this term is zero.

In the momentum equations lateral and vertical boundary conditions for
the :math:`\nabla_{h}^{2}` and
:math:`\frac{\partial^{2}}{\partial z^{2}}` operators are specified when
the numerical simulation is run - see section :numref:`eg_fourl_code_config`
For temperature the boundary condition is “zero-flux” e.g.
:math:`\frac{\partial \theta}{\partial \varphi}= \frac{\partial \theta}{\partial \lambda}=\frac{\partial \theta}{\partial z}=0`.

Discrete Numerical Configuration
--------------------------------

The domain is discretised with a uniform grid spacing in latitude and
longitude :math:`\Delta \lambda=\Delta \varphi=1^{\circ}`, so that there
are sixty grid cells in the zonal and meridional directions. Vertically
the model is configured with four layers with constant depth,
:math:`\Delta z`, of :math:`500` m. The internal, locally orthogonal,
model coordinate variables :math:`x` and :math:`y` are initialized from
the values of :math:`\lambda`, :math:`\varphi`, :math:`\Delta \lambda`
and :math:`\Delta \varphi` in radians according to

.. math::

   \begin{aligned}
   x=r\cos(\varphi)\lambda,~\Delta x & = &r\cos(\varphi)\Delta \lambda \\
   y=r\varphi,~\Delta y &= &r\Delta \varphi\end{aligned}

The procedure for generating a set of internal grid variables from a
spherical polar grid specification is discussed in section
:numref:`spatial_discrete_horizontal_grid`.

As described in :ref:`tracer_eqns`, the time evolution of potential temperature, :math:`\theta`, 
(equation :eq:`eq_eg_fourlayer_theta` is evaluated
prognostically. The centered second-order scheme with Adams-Bashforth
time stepping described in section :numref:`sub_tracer_eqns_abII` is used
to step forward the temperature equation. Prognostic terms in the
momentum equations are solved using flux form as described in section
:numref:`flux-form_momentum_equations`.
The pressure forces that drive the fluid motions, (
:math:`\frac{\partial p^{'}}{\partial \lambda}` and
:math:`\frac{\partial p^{'}}{\partial \varphi}`), are found by summing
pressure due to surface elevation :math:`\eta` and the hydrostatic
pressure. The hydrostatic part of the pressure is diagnosed explicitly
by integrating density. The sea-surface height, :math:`\eta`, is
diagnosed using an implicit scheme. The pressure field solution method
is described in sections
:numref:`press_meth_linear`
and
:numref:`finding_the_pressure_field`.

Numerical Stability Criteria
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Laplacian viscosity coefficient, :math:`A_{h}`, is set to
:math:`400 m s^{-1}`. This value is chosen to yield a Munk layer width,

.. math::
   :label: eq_eg_fourlayer_munk_layer
   M_{w} = \pi ( \frac { A_{h} }{ \beta } )^{\frac{1}{3}}

of :math:`\approx 100`\ km. This is greater than the model resolution
in mid-latitudes
:math:`\Delta x=r \cos(\varphi) \Delta \lambda \approx 80~{\rm km}` at
:math:`\varphi=45^{\circ}`, ensuring that the frictional boundary
layer is well resolved.
The model is stepped forward with a time step
:math:`\delta t=1200`\ secs. With this time step the stability
parameter to the horizontal Laplacian friction

.. math::
   :label: eq_eg_fourlayer_laplacian_stability
   S_{l} = 4 \frac{A_{h} \delta t}{{\Delta x}^2}

evaluates to 0.012, which is well below the 0.3 upper limit for
stability for this term under ABII time-stepping.
The vertical dissipation coefficient, :math:`A_{z}`, is set to
:math:`1\times10^{-2} {\rm m}^2{\rm s}^{-1}`. The associated stability
limit

.. math::
   :label: eq_eg_fourlayer_laplacian_stability_z
   S_{l} = 4 \frac{A_{z} \delta t}{{\Delta z}^2}

evaluates to :math:`4.8 \times 10^{-5}` which is again well below the
upper limit. The values of :math:`A_{h}` and :math:`A_{z}` are also
used for the horizontal (:math:`K_{h}`) and vertical (:math:`K_{z}`)
diffusion coefficients for temperature respectively.
The numerical stability for inertial oscillations

.. math::
   :label: eq_eg_fourlayer_inertial_stability
   S_{i} = f^{2} {\delta t}^2

evaluates to :math:`0.0144`, which is well below the :math:`0.5` upper
limit for stability.
The advective CFL for a extreme maximum horizontal flow speed of
:math:`| \vec{u} | = 2 ms^{-1}`

.. math::
   :label: eq_eg_fourlayer_cfl_stability
   C_{a} = \frac{| \vec{u} | \delta t}{ \Delta x}

| evaluates to :math:`5 \times 10^{-2}`. This is well below the
  stability limit of 0.5.
| The stability parameter for internal gravity waves propagating at
  :math:`2~{\rm m}~{\rm s}^{-1}`

.. math::
   :label: eq_eg_fourlayer_igw_stability
   S_{c} = \frac{c_{g} \delta t}{ \Delta x}

evaluates to :math:`\approx 5 \times 10^{-2}`. This is well below the
linear stability limit of 0.25.

.. _eg_baroclinic_code_config:

Code Configuration
------------------

The model configuration for this experiment resides under the directory
:filelink:`verification/tutorial_baroclinic_gyre/`. The experiment files 

-  :code:`input/data`

-  :code:`input/data.pkg`

-  :code:`input/eedata,`

-  :code:`input/windx.sin_y,`

-  :code:`input/topog.box,`

-  :code:`code/CPP_EEOPTIONS.h`

-  :code:`code/CPP_OPTIONS.h,`

-  :code:`code/SIZE.h.`

contain the code customisations and parameter settings for this
experiment. Below we describe the customisations to these files
associated with this experiment.

File input/data
~~~~~~~~~~~~~~~

This file, reproduced completely below, specifies the main parameters
for the experiment. The parameters that are significant for this
configuration are

-  Line 4,

   ::

       tRef=20.,10.,8.,6., 

   this line sets the initial and reference values of potential
   temperature at each model level in units of
   :math:`^{\circ}\mathrm{C}`. The entries are ordered from surface to
   depth. For each depth level the initial and reference profiles will
   be uniform in :math:`x` and :math:`y`. The values specified here are
   read into the variable in the model code, by procedure

-  Line 6,

   ::

       viscAz=1.E-2, 

   this line sets the vertical Laplacian dissipation coefficient to
   :math:`1
   \times 10^{-2} {\rm m^{2}s^{-1}}`. Boundary conditions for this
   operator are specified later. The variable is read in the routine and
   is copied into model general vertical coordinate variable At each
   time step, the viscous term contribution to the momentum equations is
   calculated in routine

-  Line 7,

   ::

      viscAh=4.E2,

   this line sets the horizontal laplacian frictional dissipation
   coefficient to :math:`1 \times 10^{-2} {\rm m^{2}s^{-1}}`. Boundary
   conditions for this operator are specified later. The variable is
   read in the routine and applied in routine .

-  Line 8,

   ::

      no_slip_sides=.FALSE.

   this line selects a free-slip lateral boundary condition for the
   horizontal laplacian friction operator e.g. :math:`\frac{\partial
       u}{\partial y}`\ =0 along boundaries in :math:`y` and
   :math:`\frac{\partial
       v}{\partial x}`\ =0 along boundaries in :math:`x`. The variable
   is read in the routine and the boundary condition is evaluated in
   routine

-  Lines 9,

   ::

      no_slip_bottom=.TRUE.

   this line selects a no-slip boundary condition for bottom boundary
   condition in the vertical laplacian friction operator e.g.
   :math:`u=v=0` at :math:`z=-H`, where :math:`H` is the local depth of
   the domain. The variable is read in the routine and is applied in the
   routine .

-  Line 10,

   ::

      diffKhT=4.E2,

   this line sets the horizontal diffusion coefficient for temperature
   to :math:`400\,{\rm m^{2}s^{-1}}`. The boundary condition on this
   operator is
   :math:`\frac{\partial}{\partial x}=\frac{\partial}{\partial y}=0` at
   all boundaries. The variable is read in the routine and used in
   routine .

-  Line 11,

   ::

      diffKzT=1.E-2,

   this line sets the vertical diffusion coefficient for temperature to
   :math:`10^{-2}\,{\rm m^{2}s^{-1}}`. The boundary condition on this
   operator is :math:`\frac{\partial}{\partial z}` = 0 on all
   boundaries. The variable is read in the routine . It is copied into
   model general vertical coordinate variable which is used in routine .

-  Line 13,

   ::

      tAlpha=2.E-4,

   This line sets the thermal expansion coefficient for the fluid to
   :math:`2
     \times 10^{-4}\,{\rm degrees}^{-1}` The variable is read in the
   routine . The routine makes use of tAlpha.

-  Line 18,

   ::

      eosType='LINEAR'

   This line selects the linear form of the equation of state. The
   variable is read in the routine . The values of eosType sets which
   formula in routine FIND_RHO is used to calculate density.

-  Line 40,

   ::

      usingSphericalPolarGrid=.TRUE.,

   This line requests that the simulation be performed in a spherical
   polar coordinate system. It affects the interpretation of grid input
   parameters, for example delX and delY and causes the grid generation
   routines to initialize an internal grid based on spherical polar
   geometry. The variable is read in the routine . When set to .TRUE.
   the settings of delX and delY are taken to be in degrees. These
   values are used in the routine

-  Line 41,

   ::

      ygOrigin=0.,

   This line sets the southern boundary of the modeled domain to
   :math:`0^{\circ}` latitude. This value affects both the generation of
   the locally orthogonal grid that the model uses internally and
   affects the initialization of the coriolis force. Note - it is not
   required to set a longitude boundary, since the absolute longitude
   does not alter the kernel equation discretisation. The variable is
   read in the routine and is used in routine

-  Line 42,

   ::

      delX=60*1.,

   This line sets the horizontal grid spacing between each y-coordinate
   line in the discrete grid to :math:`1^{\circ}` in longitude. The
   variable is read in the routine and is used in routine

-  Line 43,

   ::

      delY=60*1.,

   This line sets the horizontal grid spacing between each y-coordinate
   line in the discrete grid to :math:`1^{\circ}` in latitude. The
   variable is read in the routine and is used in routine

-  Line 44,

   ::

      delZ=500.,500.,500.,500.,

   This line sets the vertical grid spacing between each z-coordinate
   line in the discrete grid to :math:`500\,{\rm m}`, so that the total
   model depth is :math:`2\,{\rm km}`. The variable is read in the
   routine . It is copied into the internal model coordinate variable
   which is used in routine

-  Line 47,

   ::

      bathyFile='topog.box'

   This line specifies the name of the file from which the domain
   bathymetry is read. This file is a two-dimensional (:math:`x,y`) map
   of depths. This file is assumed to contain 64-bit binary numbers
   giving the depth of the model at each grid cell, ordered with the x
   coordinate varying fastest. The points are ordered from low
   coordinate to high coordinate for both axes. The units and
   orientation of the depths in this file are the same as used in the
   MITgcm code. In this experiment, a depth of :math:`0m` indicates a
   solid wall and a depth of :math:`-2000m` indicates open ocean. The
   matlab program input/gendata.m shows an example of how to generate a
   bathymetry file. The variable is read in the routine . The bathymetry
   file is read in the routine

-  Line 50,

   ::

      zonalWindFile='windx.sin_y'

   This line specifies the name of the file from which the x-direction
   (zonal) surface wind stress is read. This file is also a
   two-dimensional (:math:`x,y`) map and is enumerated and formatted in
   the same manner as the bathymetry file. The matlab program
   input/gendata.m includes example code to generate a valid
   zonalWindFile file. The variable is read in the routine . The
   wind-stress file is read in the routine

other lines in the file input/data are standard values.

<PRE>

</PRE>

File input/data.pkg
~~~~~~~~~~~~~~~~~~~

This file uses standard default values and does not contain
customisations for this experiment.

File input/eedata
~~~~~~~~~~~~~~~~~

This file uses standard default values and does not contain
customisations for this experiment.

File input/windx.sin_y
~~~~~~~~~~~~~~~~~~~~~~

The input/windx.sin_y file specifies a two-dimensional (:math:`x,y`) map
of wind stress ,\ :math:`\tau_{x}`, values. The units used are
:math:`Nm^{-2}` (the default for MITgcm). Although :math:`\tau_{x}` is
only a function of latitude, :math:`y`, in this experiment this file
must still define a complete two-dimensional map in order to be
compatible with the standard code for loading forcing fields in MITgcm
(routine EXTERNAL_FIELDS_LOAD. The included matlab program
input/gendata.m gives a complete code for creating the input/windx.sin_y
file.

File input/topog.box
~~~~~~~~~~~~~~~~~~~~

The input/topog.box file specifies a two-dimensional (:math:`x,y`) map
of depth values. For this experiment values are either :math:`0~{\rm m}`
or :math:`-2000\,{\rm m}`, corresponding respectively to a wall or to
deep ocean. The file contains a raw binary stream of data that is
enumerated in the same way as standard MITgcm two-dimensional,
horizontal arrays. The included matlab program input/gendata.m gives a
complete code for creating the input/topog.box file.

File code/SIZE.h
~~~~~~~~~~~~~~~~

Two lines are customized in this file for the current experiment

-  Line 39,

   ::

       sNx=60, 

   this line sets the lateral domain extent in grid points for the axis
   aligned with the x-coordinate.

-  Line 40,

   ::

       sNy=60, 

   this line sets the lateral domain extent in grid points for the axis
   aligned with the y-coordinate.

-  Line 49,

   ::

       Nr=4,   

   this line sets the vertical domain extent in grid points.

File code/CPP_OPTIONS.h
~~~~~~~~~~~~~~~~~~~~~~~

This file uses standard default values and does not contain
customisations for this experiment.

File code/CPP_EEOPTIONS.h
~~~~~~~~~~~~~~~~~~~~~~~~~

This file uses standard default values and does not contain
customisations for this experiment.

Other Files 
~~~~~~~~~~~~

Other files relevant to this experiment are

-  model/src/ini_cori.F. This file initializes the model coriolis
   variables fCorU and fCorV.

-  model/src/ini_spherical_polar_grid.F This file initializes the model
   grid discretisation variables dxF, dyF, dxG, dyG, dxC, dyC.

-  model/src/ini_parms.F.

Running The Example
-------------------

Code Download
~~~~~~~~~~~~~

In order to run the examples you must first download the code
distribution. Instructions for downloading the code can be found in
section `[sec:obtainingCode] <#sec:obtainingCode>`__.

Experiment Location
~~~~~~~~~~~~~~~~~~~

This example experiments is located under the release sub-directory

Running the Experiment
~~~~~~~~~~~~~~~~~~~~~~

To run the experiment

#. Set the current directory to input/

   ::

      % cd input

#. Verify that current directory is now correct

   ::

      % pwd

   You should see a response on the screen ending in

   verification/exp2/input

#. Run the genmake script to create the experiment Makefile

   ::

      % ../../../tools/genmake -mods=../code

#. Create a list of header file dependencies in Makefile

   ::

      % make depend

#. Build the executable file.

   ::

      % make

#. Run the mitgcmuv executable

   ::

      % ./mitgcmuv
