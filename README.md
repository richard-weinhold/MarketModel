
<img  height="24" src="https://raw.githubusercontent.com/richard-weinhold/pomato/main/docs/_static/graphics/pomato_logo_small.png"> POMATO - Power Market Tool <img  height="24" src="https://raw.githubusercontent.com/richard-weinhold/pomato/main/docs/_static/graphics/pomato_logo_small.png">
=========================================================================================================================================================

Main Branch: ![MarketModel](https://github.com/richard-weinhold/MarketModel/workflows/MarketModel/badge.svg)
![codecov](https://codecov.io/gh/richard-weinhold/MarketModel/branch/master/graph/badge.svg?token=qUK3at8Am2)

Construction Branch: 
![MarketModel](https://github.com/richard-weinhold/MarketModel/actions/workflows/MarketModel_testing.yml/badge.svg?branch=construction)
![codecov](https://codecov.io/gh/richard-weinhold/MarketModel/branch/construction/graph/badge.svg?token=qUK3at8Am2)

Overview
--------

The Power Market Tool [(POMATO)](https://github.com/richard-weinhold/pomato) aims to enable research and analyses of modern real-world electricity markets in the context of the physical transmission system and its secure operation. To analyze relevant market and operational processes POMATO provides: Simultaneous zonal market clearing and nodal power flow computation, a scalable and fast SCOPF algorithm, and risk-aware optimal power flow via chance constraints. All optimization features rely on Julia/JuMP, leveraging its accessibility, computational performance, and solver interfaces.

This repository holds these features as a individual julia package for easier maintenance, testing and documentation. However it is eventually embedded in a Python front-end, providing flexible and easily maintainable data processing and user interaction features and is not meant so be used on its own. Please see the corresponding publication and the [(POMATO repository)](https://github.com/richard-weinhold/pomato) for further information.


Installation
------------

This package is meant to be used in conjunction with the python
[POMATO](https://github.com/richard-weinhold/pomato) model, which embeds its features and installs
it automatically. For stand alone usage see the documentation page. 

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://richard-weinhold.github.io/MarketModel/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://richard-weinhold.github.io/MarketModel/dev)


Related Publications
--------------------
- (*preprint*) [Weinhold and Mieth (2020), Power Market Tool (POMATO) for the Analysis of Zonal 
   Electricity Markets](https://arxiv.org/abs/2011.11594)
- [Weinhold and Mieth (2020), Fast Security-Constrained Optimal Power Flow through 
   Low-Impact and Redundancy Screening](https://ieeexplore.ieee.org/document/9094021)
- [Schönheit, Weinhold, Dierstein (2020), The impact of different strategies for generation 
   shift keys (GSKs) on  the flow-based market coupling domain: A model-based analysis of Central Western Europe](https://www.sciencedirect.com/science/article/pii/S0306261919317544)

