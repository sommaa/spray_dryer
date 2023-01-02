<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/sommaa/spray_dryer">
    <img src="https://user-images.githubusercontent.com/120776791/210235026-1aae8c7f-3931-4496-b92d-35e1187e60de.png" alt="Logo" width="450" height="240">
  </a>
      <br />
    <h1 align="center">Spray Dryer</h3>
    <h4 align="center">Mathematical modeling of a spray dryer in matlab and OpenFoam (only flow-dynamic), modified from https://super.chem.polimi.it/ script, detailed informations in pdf folder.</h4>
</div>

## 1-D model
Simulation details [Matlab-new](./spray_matlab.m):
- non-dimensional numbers calculated inside the integral function;
- constrained evaluation of mass and particle diameter;
- weighted evaluation of density;

Simulation details [Matlab-corrected](./spray_matlab_corr.m):
- non-dimensional numbers calculated inside the integral function;
- constrained evaluation of mass and particle diameter;
- weighted evaluation of density;
- surface confined axial jet correlation from [Hydraulic characteristics of turbulent circular jets under surface confinement](https://doi.org/10.1080/09715010.2013.876725), see [jet folder](./jet) for numerical representation.
## 3D-model
Simulation details [RAS](./RAS):
- 2Mln mesh cells;
- k-omega-SST RAS model;
- inlet Q=20Kg/s;
- outlet Pout = 0;

Simulation details [RAS + DPFoam](./PART):
- 2Mln mesh cells;
- k-omega-SST RAS model;
- inlet Q=20Kg/s;
- outlet Pout = 0;
- 8000 particles;
- Drag Force + gravity + boundaries rebound/escape;

## Clone it!
```
$ git clone https://github.com/sommaa/spray_dryer.git
```
# Results:
- plot [Matlab-new](./spray_matlab.m):

![new](https://user-images.githubusercontent.com/120776791/210239070-f98a9f94-06b4-4d32-8390-32ecc7b09e4d.png)

- plot [Matlab-corrected](./spray_dryer_corrected) for D-ratio = 1.5

![corr](https://user-images.githubusercontent.com/120776791/210239108-a7332387-e0ab-49cf-8ee4-2ab084810a67.png)

- 3D flow in OpenFoam:

![bitmap](https://user-images.githubusercontent.com/120776791/210233646-6c381613-8675-4e05-a793-e01691cfe480.png)

- matrix particle plot

![matrix_view](https://user-images.githubusercontent.com/120776791/210233664-ca56476a-c940-4ee9-b692-c8fd118b6ad2.png)

- turbulent jet in OpenLB:

![jet](https://user-images.githubusercontent.com/120776791/210234167-184b5964-2845-4d13-9c66-8c63ec6dc062.png)

