<img  height="24" src="https://raw.githubusercontent.com/richard-weinhold/pomato/main/docs/_static/graphics/pomato_logo_small.png"> POMATO - Power Market Tool <img  height="24" src="https://raw.githubusercontent.com/richard-weinhold/pomato/main/docs/_static/graphics/pomato_logo_small.png">
=========================================================================================================================================================

![MarketModel](https://github.com/richard-weinhold/MarketModel/workflows/MarketModel/badge.svg)

Overview
--------

The Power Market Tool [(POMATO)](https://github.com/richard-weinhold/pomato) aims to enable research and analyses of modern real-world electricity markets in the context of the physical transmission system and its secure operation. To analyze relevant market and operational processes POMATO provides: Simultaneous zonal market clearing and nodal power flow computation, a scalable and fast SCOPF algorithm, and risk-aware optimal power flow via chance constraints. All optimization features rely on Julia/JuMP, leveraging its accessibility, computational performance, and solver interfaces.

This repository holds these features as a individual julia package for easier maintenance, testing and documentation. However it is eventually embedded in a Python front-end, providing flexible and easily maintainable data processing and user interaction features and is not meant so be used on its own. Please see the corresponding publication and the [(POMATO repository)](https://github.com/richard-weinhold/pomato) for further information.

Installation
------------

The MarketModel can be cloned and added to you julia projects or just used from the repository. The model requires Julia 1.3 and works with the open Clp solver. As said before it is meant to be using in conjunction with the python [(POMATO)](https://github.com/richard-weinhold/pomato) model, which embeds its features and installs it automatically.

Related Publications
--------------------
- (*preprint*) [Weinhold and Mieth (2020), Power Market Tool (POMATO) for the Analysis of Zonal 
   Electricity Markets](https://arxiv.org/abs/2011.11594)
- [Weinhold and Mieth (2020), Fast Security-Constrained Optimal Power Flow through 
   Low-Impact and Redundancy Screening](https://ieeexplore.ieee.org/document/9094021)
- [Schönheit, Weinhold, Dierstein (2020), The impact of different strategies for generation 
   shift keys (GSKs) on  the flow-based market coupling domain: A model-based analysis of Central Western Europe](https://www.sciencedirect.com/science/article/pii/S0306261919317544)

