<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/sommaa/spray_dryer">
    <img src="./image/blender_ext.png" alt="Logo" width="450" height="240">
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

Simulation details [Matlab-corrected](./spray_dryer_corr):
- non-dimensional numbers calculated inside the integral function;
- constrained evaluation of mass and particle diameter;
- weighted evaluation of density;
- surface confined axial jet correlation from [Hydraulic characteristics of turbulent circular jets under surface confinement](https://doi.org/10.1080/09715010.2013.876725)
## 3D-model
Simulation details [RAS](./spray_RAS):
- 2Mln mesh cells;
- k-omega-SST RAS model;
- inlet Q=20Kg/s;
- outlet Pout = 0;

Simulation details [RAS + DPFoam](./spray_part):
- 2Mln mesh cells;
- k-omega-SST RAS model;
- inlet Q=20Kg/s;
- outlet Pout = 0;
- 8000 particles;
- Drag Force + gravity + boundaries rebound/escape;

# Results:
- plot [Matlab-new](./spray_matlab.m):

![SI](https://user-images.githubusercontent.com/120776791/209564200-ea86d385-22e0-42ba-93f6-3a1c8d237bb8.svg)

- plot [Matlab-corrected](./spray_dryer_corrected) for D-ratio = 1.5

![1 5](https://user-images.githubusercontent.com/120776791/209564053-b41b8e1c-f7ac-46d1-89bb-2a8bc1c4743d.svg)

- 3D flow in OpenFoam:

![spray_time](https://user-images.githubusercontent.com/120776791/208314061-622f8a5d-346f-4b20-85b9-ec6a555f7d24.png)

![spray_gly](https://user-images.githubusercontent.com/120776791/208314063-88906cc5-339f-478f-a1a1-0c2c6802754f.png)

![spray_tracer](https://user-images.githubusercontent.com/120776791/208314248-bc7ea047-9b46-487f-a056-4fa3db7421a4.png)

https://user-images.githubusercontent.com/120776791/208930344-1c460252-828a-4064-b2e2-bc88c6fd3eec.mp4

https://user-images.githubusercontent.com/120776791/208931821-2a182562-6b0d-4896-b570-e6eeabedb428.mp4

- matrix particle plot

![matrix_view](https://user-images.githubusercontent.com/120776791/209577374-08aca8b9-f289-4600-b0a1-6500061462bc.png)

- mesh

![spray_mesh](https://user-images.githubusercontent.com/120776791/208314557-cd86c945-0d5a-4090-be0d-0090b2b051e3.png)
