# pocketTHERM

pocketTHERM is a suite of web-based calculators that can be used to model organic Rankine cycle (ORC) and heat pump systems. The tool is primarily intended as an educational tool that can be run directly from a web-browser without needing any form of installation.

---

### Features

pocketTHERM contains calculators to model ORC and heat pump systems. These models combine functions to predict the thermophysical properties of common working fluids with functions that predict the performance of the thermodynamic system.

For a defined a heat source and heat sink stream, the user can vary the working fluid and key cycle design variables and see the impact of these variables on system performance in real-time. The tool enables both single-point and parametric studies to be completed.  

Alongside this, pocketTHERM includes interactive self-paced lessons that enable users to develop their understanding of these technologies.

---

### Implementation

pocketTHERM is written in HTML and Python. All the underlying calculations are computed using custom Python modules. The calculators and lessons each have their own webpage, and these are coupled to an intermediate Python script that interfaces between HTML and the underlying Python modules. PyScript is used to call the Python scripts directly from HTML.

More information on PyScript: [https://pyscript.net/](https://pyscript.net/)

---

### Pre-requisites

There are no pre-requisites. pocketTHERM runs directly from the web-browser.

---

### How to get started
Go the pocketTHERM webpage: [https://pockettherm.github.io/](https://pockettherm.github.io/)

---

### Developers
**Martin T. White**, Senior Lecturer in Mechanical Engineering, University of Sussex, [martin.white@sussex.ac.uk](mailto:martin.white@sussex.ac.uk)
