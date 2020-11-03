[![PomatoLogo](https://github.com/richard-weinhold/pomato/blob/master/docs/pomato_logo_small.png "Pomato Soup")](#) MarketModel package for use in POMATO, a Power Market Tool For The Analysis Of Modern Electricity Markets
=====================================================================================================================================
[![Build Status](https://travis-ci.org/richard-weinhold/MarketModel.svg?branch=master)](https://travis-ci.org/richard-weinhold/MarketModel)

Overview
--------

The Power Market Tool [(POMATO)](https://github.com/richard-weinhold/pomato) aims to enable research and analyses of modern real-world electricity markets in the context of the physical transmission system and its secure operation. To analyze relevant market and operational processes POMATO provides: Simultaneous zonal market clearing and nodal power flow computation, a scalable and fast SCOPF algorithm, and risk-aware optimal power flow via chance constraints. All optimization features rely on Julia/JuMP, leveraging its accessibility, computational performance, and solver interfaces.

This repository holds these features as a individual julia package for easier maintenance, testing and documentation. However it is eventually embedded in a Python front-end, providing flexible and easily maintainable data processing and user interaction features and is not meant so be used on its own. Please see the corresponding publication and the [(POAMTO repository)](https://github.com/richard-weinhold/pomato) for further information.

Installation
------------

The MarketModel can be cloned and added to you julia projects or just used from the repository. The model requires Julia 1.3 and works with the open Clp solver. As said before it is meant tobe using in conjunction with the python [(POAMTO)](https://github.com/richard-weinhold/pomato) model, which embeds its features and installs it automatically.

Related Publications
--------------------

- [Weinhold and Mieth (2020), Fast Security-Constrained Optimal Power Flow through
   Low-Impact and Redundancy Screening](https://ieeexplore.ieee.org/document/9094021)
- [Schönheit, Weinhold, Dierstein (2020), The impact of different strategies for generation shift keys (GSKs) on the flow-based market coupling domain: A model-based analysis of Central Western Europe](https://www.sciencedirect.com/science/article/pii/S0306261919317544)
